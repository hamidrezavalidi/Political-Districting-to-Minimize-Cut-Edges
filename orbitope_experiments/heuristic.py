###########################
# Run options
###########################  

levels = { 'county', 'tract' }
iteration_options = { 100, 1000, 10000 }

###########################
# Imports
###########################  

import time
import json
import os

from gerrychain import (GeographicPartition, Graph, MarkovChain, updaters, constraints, accept)
from gerrychain.tree import recursive_tree_part
from gerrychain.proposals import recom
from functools import partial

import geopandas as gpd
import matplotlib.pyplot as plt

###########################
# Hard-coded inputs
###########################  

state_codes = {
    'WA': '53', 'DE': '10', 'WI': '55', 'WV': '54', 'HI': '15',
    'FL': '12', 'WY': '56', 'NJ': '34', 'NM': '35', 'TX': '48',
    'LA': '22', 'NC': '37', 'ND': '38', 'NE': '31', 'TN': '47', 'NY': '36',
    'PA': '42', 'AK': '02', 'NV': '32', 'NH': '33', 'VA': '51', 'CO': '08',
    'CA': '06', 'AL': '01', 'AR': '05', 'VT': '50', 'IL': '17', 'GA': '13',
    'IN': '18', 'IA': '19', 'MA': '25', 'AZ': '04', 'ID': '16', 'CT': '09',
    'ME': '23', 'MD': '24', 'OK': '40', 'OH': '39', 'UT': '49', 'MO': '29',
    'MN': '27', 'MI': '26', 'RI': '44', 'KS': '20', 'MT': '30', 'MS': '28',
    'SC': '45', 'KY': '21', 'OR': '41', 'SD': '46'
}

congressional_districts = {
    'WA': 10, 'DE': 1, 'WI': 8, 'WV': 3, 'HI': 2,
    'FL': 27, 'WY': 1, 'NJ': 12, 'NM': 3, 'TX': 36,
    'LA': 6, 'NC': 13, 'ND': 1, 'NE': 3, 'TN': 9, 'NY': 27,
    'PA': 18, 'AK': 1, 'NV': 4, 'NH': 2, 'VA': 11, 'CO': 7,
    'CA': 53, 'AL': 7, 'AR': 4, 'VT': 1, 'IL': 18, 'GA': 14,
    'IN': 9, 'IA': 4, 'MA': 9, 'AZ': 9, 'ID': 2, 'CT': 5,
    'ME': 2, 'MD': 8, 'OK': 5, 'OH': 16, 'UT': 4, 'MO': 8,
    'MN': 8, 'MI': 14, 'RI': 2, 'KS': 4, 'MT': 1, 'MS': 4,
    'SC': 7, 'KY': 6, 'OR': 5, 'SD': 1
}

skips = {
    ('WA','tract'), ('WA','county'), ('DE','tract'), ('DE','county'), ('WI','tract'), 
    ('WI','county'), ('HI','tract'), ('HI','county'), ('FL','tract'), ('FL','county'), 
    ('WY','tract'), ('WY','county'), ('NJ','tract'), ('NJ','county'), ('TX','tract'), 
    ('TX','county'), ('LA','tract'), ('LA','county'), ('NC','tract'), ('NC','county'), 
    ('ND','tract'), ('ND','county'), ('TN','tract'), ('TN','county'), ('NY','tract'), 
    ('NY','county'), ('PA','tract'), ('PA','county'), ('AK','tract'), ('AK','county'), 
    ('NV','county'), ('NH','county'), ('VA','tract'), ('VA','county'), ('CO','tract'), 
    ('CO','county'), ('CA','tract'), ('CA','county'), ('AL','tract'), ('VT','tract'), 
    ('VT','county'), ('IL','tract'), ('IL','county'), ('GA','tract'), ('GA','county'),
    ('IN','tract'), ('IN','county'), ('IA','tract'), ('MA','tract'), ('MA','county'), 
    ('AZ','tract'), ('AZ','county'), ('CT','tract'), ('CT','county'), ('ME','county'), 
    ('MD','tract'), ('MD','county'), ('OK','tract'), ('OH','tract'), ('OH','county'), 
    ('UT','county'), ('MO','tract'), ('MO','county'), ('MN','tract'), ('MN','county'), 
    ('MI','tract'), ('MI','county'), ('RI','tract'), ('RI','county'), ('KS','tract'), 
    ('MT','tract'), ('MT','county'), ('SC','tract'), ('SC','county'), ('KY','tract'), 
    ('KY','county'), ('OR','tract'), ('OR','county'), ('SD','tract'), ('SD','county')
}

################################################
# Draws districts and saves to png file
################################################ 

def export_to_png(G, df, districts, filename):
    
    assignment = [ -1 for u in G.nodes ]
    
    for j in range(len(districts)):
        for i in districts[j]:
            geoID = G.nodes[i]["GEOID10"]
            for u in G.nodes:
                if geoID == df['GEOID10'][u]:
                    assignment[u] = j
    
    if min(assignment[v] for v in G.nodes) < 0:
        print("Error: did not assign all nodes in district map png.")
    else:
        df['assignment'] = assignment
        my_fig = df.plot(column='assignment').get_figure()
        RESIZE_FACTOR = 3
        my_fig.set_size_inches(my_fig.get_size_inches()*RESIZE_FACTOR)
        plt.axis('off')
        my_fig.savefig(filename)
        

####################################
# Function for GerryChain call
####################################                       

def run_GerryChain_heuristic(G,population_deviation,k,iterations):
    
    my_updaters = {"population": updaters.Tally("TOTPOP", alias="population")}
    start = recursive_tree_part(G,range(k),sum(G.nodes[i]["TOTPOP"] for i in G.nodes())/k,"TOTPOP", population_deviation/2,1)
    initial_partition = GeographicPartition(G, start, updaters = my_updaters)
    
    proposal = partial(recom,
                       pop_col="TOTPOP",
                       pop_target=sum(G.nodes[i]["TOTPOP"] for i in G.nodes())/k,
                       epsilon=population_deviation/2,
                       node_repeats=2
                      )
    
    compactness_bound = constraints.UpperBound(
        lambda p: len(p["cut_edges"]),
        1.5*len(initial_partition["cut_edges"])
    )
    
    pop_constraint = constraints.within_percent_of_ideal_population(initial_partition, population_deviation/2)
    
    my_chain = MarkovChain(
        proposal=proposal,
        constraints=[
            pop_constraint,
            compactness_bound
        ],
        accept=accept.always_accept,
        initial_state=initial_partition,
        total_steps=iterations
    )
    
    min_cut_edges = sum(G[i][j]['edge_length'] for i,j in G.edges)
    print("In GerryChain heuristic, current # of cut edges: ",end='')
    print(min_cut_edges,",",sep='',end=' ')
    for partition in my_chain:
        current_cut_edges = sum(G[i][j]['edge_length'] for i,j in partition["cut_edges"])
        print(current_cut_edges,",",sep='',end=' ')
        if current_cut_edges < min_cut_edges:
            best_partition = partition
            min_cut_edges = current_cut_edges
    
    print("Best heuristic solution has # cut edges =",min_cut_edges)
    return ([[i for i in G.nodes if best_partition.assignment[i]==j] for j in range(k)],min_cut_edges)


###########################
# Main part of the code
###########################  

# create directories for results
os.mkdir("../heuristic-results")
for iterations in iteration_options:
    os.mkdir("../heuristic-results/"+str(iterations)+"-iterations") 

# run all settings
for state in state_codes.keys():
    
    # parameters            
    k = congressional_districts[state]
    deviation = 0.01
    code = state_codes[state]
    
    for level in levels:
        
        # skip certain (state,level) pairs that we know:
        #   1. are infeasible, 
        #   2. are outside our scope (because of size), or
        #   3. gerrychain gets stuck on (infinite loop).
        
        if (state,level) in skips:
            continue
        
        # read input graph and shapefile df
        G = Graph.from_json("../data/"+level+"/dual_graphs/"+level+code+".json")
        df = gpd.read_file("../data/"+level+"/shape_files/"+state+"_"+level+".shp")
        
        # give each edge a "length" of one
        for i,j in G.edges:
            G[i][j]['edge_length'] = 1
        
        for iterations in iteration_options:
        
            # run GerryChain 
            start = time.time()
            (districts,heur_obj) = run_GerryChain_heuristic(G,deviation,k,iterations)
            stop = time.time()
            
            # filename for outputs
            fn = "../heuristic-results/"+str(iterations)+"-iterations/heur_"+state+"_"+level
            
            # draw the solution on a map
            png_fn = fn + ".png"
            export_to_png(G, df, districts, png_fn)
            
            # dump the solution info to json file
            json_fn = fn + ".json"
            with open(json_fn, 'w') as outfile:
                data = {}
                data['obj'] = heur_obj
                data['time'] = '{0:.2f}'.format(stop-start)
                data['iterations'] = iterations
                data['nodes'] = list()
        
                for j in range(k):
                    for i in districts[j]:
                        data['nodes'].append({
                                'name': G.nodes[i]["NAME10"],
                                'index': i,
                                'district': j
                                })
                json.dump(data, outfile)
                
