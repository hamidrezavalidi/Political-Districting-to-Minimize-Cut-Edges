#import networkx as nx
import time
from gerrychain import (GeographicPartition, Partition, Graph, MarkovChain, proposals, updaters, constraints, accept, Election)
from gerrychain.proposals import recom
from functools import partial

import networkx as nx

#from gerrychain.metrics import efficiency_gap, mean_median
#from gerrychain.proposals import propose_random_flip
#from gerrychain.updaters import cut_edges
from gerrychain.tree import recursive_tree_part

import geopandas as gpd
import matplotlib.pyplot as plt

from matplotlib.colors import LinearSegmentedColormap

import json


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

def set_edge_lengths(G, weighted):
    if weighted:
        for i,j in G.edges:
            #print("Distance between county", G.node[i]['NAME10'], "and", G.node[j]['NAME10'],"is",G[i][j]['shared_perim'])
            G[i][j]['edge_length'] = round(G[i][j]['shared_perim']*100)
    else:
        for i,j in G.edges:
            G[i][j]['edge_length'] = 1
            
            
land_parcel = 'tract'            
            
#weighted = False     # weight the edges based on border lengths

         

def run_GerryChain_heuristic(G,population_deviation,k,iterations):
    
    my_updaters = {"population": updaters.Tally("TOTPOP", alias="population")}
    start = recursive_tree_part(G,range(k),sum(G.node[i]["TOTPOP"] for i in G.nodes())/k,"TOTPOP", population_deviation/2,1)
    initial_partition = GeographicPartition(G, start, updaters = my_updaters)
    
    proposal = partial(recom,
                       pop_col="TOTPOP",
                       pop_target=sum(G.node[i]["TOTPOP"] for i in G.nodes())/k,
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
    
    print("Best heuristic solution has (possibly weighted) # cut edges =",min_cut_edges)
    return ([[i for i in G.nodes if best_partition.assignment[i]==j] for j in range(k)],min_cut_edges)


if land_parcel == 'county':
    heuristic_iterations = 10000
elif land_parcel == 'tract':
    heuristic_iterations = 100
else: print("Error: land parcel is not valid!")    

weighted = False

for state in state_codes.keys():
    if state == 'UT' or state == 'MS' or state == 'AR' or state == 'NV':
    #if state == 'MS':
    #if state == 'OK' or state == 'AL' or state == 'NE' or state == 'AR' or state == 'ID' or state == 'WV' or state == 'NM' or state == 'MS':
    #if state == 'OK' or state == 'AL' or state == 'NE' or state == 'AR' or state == 'KS' or state == 'IA' or state == 'ID' or state == 'MS':
    #if state == 'OK' or state == 'AL' or state == 'AR' or state == 'KS' or state == 'IA' or state == 'ID' or state == 'MS':  
    #if state == 'OK' or state == 'KS' or state == 'MS':
    #if state == 'NH' or state == 'ID' or state == 'ME' or state == 'WV' or state == 'NM' or state == 'NE':
            fn = "heur_"+state+"_"+land_parcel 
            if weighted:
                 fn += "_weighted"
            fn += ".json"
                        
            start = time.time()
                        
            k = congressional_districts[state]
            #k = 17
            #population_deviation = 0.1
            population_deviation = 0.01
            state_code = state_codes[state]
            
            G = Graph.from_json("C:/data/Your-State/"+land_parcel+"/dual_graphs/"+land_parcel+state_code+".json")

            #G = nx.convert_node_labels_to_integers(F, first_label=0, ordering='default', label_attribute=None)
            
            set_edge_lengths(G, weighted)
            
            (heur_districts,heur_obj) = run_GerryChain_heuristic(G,population_deviation,k,heuristic_iterations)
            stop = time.time()
            if land_parcel == 'county':
                df = gpd.read_file("C:/data/districting/"+state+"/counties/maps/"+state+"_counties.shp") 
            elif land_parcel == 'tract':
                df = gpd.read_file("C:/data/districting/"+state+"/tracts/maps/"+state+"_tracts.shp")
            else:
                print("Please enter a valid land_parcel!")
            
            # heuristic maps
            # create a new column for heuristic district labels, and fill it in
            df['heuristic']= -1
            for j in range(k):
                for i in heur_districts[j]:
                    geoID = G.node[i]["GEOID10"]
                    for u in G.nodes:
                        if geoID == df['GEOID10'][u]:
                            df['heuristic'][u] = j
            
            '''               
            cmap = LinearSegmentedColormap.from_list('heuristic', [(0, 'blue'), (0.5, 'yellow'), (1, 'green')])    
            #colors = ['grey', 'white']  
            
            splot = df.plot(cmap=cmap, column='heuristic',figsize=(10, 10), linewidth=1, edgecolor='0.25').get_figure()  # display the S map
            plt.axis('off')
            splot.savefig("heur_"+state+"_"+land_parcel+".png")
            '''
            
            hplot = df.plot(column='heuristic',figsize=(10, 10)).get_figure()
            plt.axis('off')
            if weighted:
                hplot.savefig("heur_"+state+"_"+land_parcel+"_weighted.png")
            else:
                hplot.savefig("heur_"+state+"_"+land_parcel+".png")
            
            with open(fn, 'w') as outfile:
                #open(fn, 'w', newline='')
                data = {}
                data['obj'] = heur_obj
                data['time'] = round(stop-start,2)
                data['iterations'] = heuristic_iterations
                data['nodes'] = []
                #textFile.write(str(heur_obj)+" "+str(round(stop-start,2))+"\n")
                for j in range(k):
                    #print("district: ", j)
                    for i in heur_districts[j]:
                        #print("county: ", G.node[i]["NAME10"])
                        data['nodes'].append({
                                'name': G.node[i]["NAME10"],
                                'index': i,
                                'district': j
                                })
                        #L = [str(i)," ",str(j),"\n"]
                        #textFile.writelines(L)
                json.dump(data, outfile)