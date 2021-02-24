import gurobipy as gp
from gurobipy import GRB 

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

import math
import networkx as nx
import csv
import time
import json
#from time import process_time

import main_hess
import main_labeling
import ordering
#import heuristic

from gerrychain import Graph
import geopandas as gpd


###########################
# Input codes
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

epsg_codes = {
    'AL': '26730', 'AK': '3471', 'AZ': '3478', 'AR': '3484', 'CA': '3493',
    'CO': '3501', 'CT': '3507', 'DE': '3509', 'FL': '3514', 'GA': '3518',
    'HI': '2784', 'ID': '3524', 'IL': '3528', 'IN': '3532', 'IA': '3536',
    'KS': '3540', 'KY': '3544', 'LA': '3550', 'ME': '3557', 'MD': '3559',
    'MA': '3585', 'MI': '3587', 'MN': '3594', 'MS': '3597', 'MO': '3602',
    'MT': '3604', 'NE': '3606', 'NV': '3607', 'NH': '3613', 'NJ': '3615',
    'NM': '3617', 'NY': '3623', 'NC': '3631', 'ND': '3633', 'OH': '3637',
    'OK': '3639', 'OR': '3645', 'PA': '3649', 'RI': '3653', 'SC': '3655',
    'SD': '3657', 'TN': '3661', 'TX': '3669', 'UT': '3675', 'VT': '3684',
    'VA': '3685', 'WA': '3689', 'WV': '3693', 'WI': '3695', 'WY': '3703',
}

land_parcel = 'tract'

hess_models = {'lcut'}

#hess_models = {'lcut', 'scf', 'shir'}

#labeling_models = {'mcf', 'shir', 'cut', 'lcut', 'austin_cut', 'austin_lcut'}

labeling_models = {'lcut'}

#labeling_models = {'lcut', 'scf', 'shir', 'labeling'}


def set_edge_lengths(G, weighted):
    if weighted:
        for i,j in G.edges:
            G[i][j]['edge_length'] = round(G[i][j]['shared_perim']*100)
    else:
        for i,j in G.edges:
            G[i][j]['edge_length'] = 1
            

base = "labeling"
# weight of objective coefficients
weighted = False     # weight the edges based on border lengths

# extended formulation
extended = True      # use stronger (but larger) extended model to capture objective function



# order-choice for finding an ordering
OrderChoice = "S_decreasing"  # How to order the vertices? 
    # none: no particular ordering; use whatever order is given when graph is read
    # S : at the back of the ordering, place a maximum cardinality subset S of vertices such that every component of G[S] has population less than L
    # increasing_population : order the vertices by increasing population
    # decreasing_population : order the vertices by decreasing population. Use it for base models without contiguity constraints 
    # S_decreasing : Same as S, but the V\S portion is sorted by decreasing population
    # S_indep_decreasing : Solves max independent set for finding S. V\S portion is sorted by decreasing population
    # heur_ordering : do LFix heuristically. You must have a heuristic solution

# optimal solutions
all_optima = False   # find all optimal solutions?
optima_limit = 100   # how many optimal solutions (at most) to obtain?

#symmetry parameter? (recommended for labeling-base formulations)
symmetry = 2   # options: -1 (automatic), 0 (no symmetry), 1 (conservative), 2 (aggressive)



def build_cut_edges(state,base,model,weighted,extended,OrderChoice,all_optima,optima_limit,symmetry):
    use_heuristic = True # run GerryChain to find a heuristic solution?
    #total_start = time.time() 
    m = gp.Model()
    ###########################
    # READ MORE INPUTS
    ###########################
    m._state = state
    m._model = model
    m._extended = extended
    m._weighted = weighted
    #m._use_heuristic = use_heuristic
    # extended formulation
    m._extended = extended      # use stronger (but larger) extended model to capture objective function
    # contiguity model
    #m._model = "austin_lcut"        # options: hess, shir, mcf, cut, lcut
    #symmetry parameter? (recommended for labeling-base formulations)
    m.Params.symmetry = symmetry   # options: -1 (automatic), 0 (no symmetry), 1 (conservative), 2 (aggressive)
    # optimal solutions
    m._all_optima = all_optima   # find all optimal solutions?
    
        
    # fix me?        
    #if weighted==True:
     #   print("ERROR: weighted edges not yet supported.")
        
    k = congressional_districts[state]
    #k = 17
    #population_deviation = 0.1
    population_deviation = 0.01
    state_code = state_codes[state]
    '''
    if land_parcel == 'county':
        land = 'counties'
    else:
        land = 'tracts'
    '''    
    G = Graph.from_json("C:/data/Your-State/"+land_parcel+"/dual_graphs/"+land_parcel+state_code+".json")
    #G = nx.read_edgelist("C:/data/districting/"+state+"/"+land+"/graph/"+state+".txt",nodetype=int)
    #population = read_population("C:/data/districting/"+state+"/"+land+"/graph/"+state+".population")
    if land_parcel == 'county':
        m._df = gpd.read_file("C:/data/districting/"+state+"/counties/maps/"+state+"_counties.shp") 
        if weighted and use_heuristic:
            heur_data = open("C:/data/Heuristic/JSON/"+land_parcel+"/Weighted/"+"heur_"+state+"_"+land_parcel+"_weighted.json",)
        elif not weighted and use_heuristic:
            heur_data = open("C:/data/Heuristic/JSON/"+land_parcel+"/Unweighted/"+"heur_"+state+"_"+land_parcel+".json",)
    elif land_parcel == 'tract' and use_heuristic:
        m._df = gpd.read_file("C:/data/districting/"+state+"/tracts/maps/"+state+"_tracts.shp")
        if weighted and use_heuristic:
            heur_data = open("C:/data/Heuristic/JSON/"+land_parcel+"/Weighted/"+"heur_"+state+"_"+land_parcel+"_weighted.json",)
        elif not weighted and use_heuristic:
            heur_data = open("C:/data/Heuristic/JSON/"+land_parcel+"/Unweighted/"+"heur_"+state+"_"+land_parcel+".json",)
    elif land_parcel == 'tract' and not use_heuristic:
        m._df = gpd.read_file("C:/data/districting/"+state+"/tracts/maps/"+state+"_tracts.shp")
    else:
        print("Please enter a valid land_parcel!")
    if use_heuristic:    
        heur = json.load(heur_data) 
    else:
        heur = None
    #print("EPSG:"+"epsg_codes[state]")
    #m._df = m._df.to_crs("EPSG:"+epsg_codes[state])
    #print("data frame: ", m._df)
    ###########################
    # MODEL PARAMETERS
    ###########################
    # read heuristic files
    # Opening JSON file 
    
    #f = open('data.json',) 
  
    # returns JSON object as  
    # a dictionary 
    #data = json.load(f) 
    
    m._use_heuristic = use_heuristic  # run GerryChain to find a heuristic solution?
    #if m._state == 'ME':
     #   heuristic_iterations = 1
    #else:    
     #   heuristic_iterations = 5000 # of iterations of GerryChain   
    population = [G.node[i]["TOTPOP"] for i in G.nodes()]
    L = math.ceil((1-population_deviation/2)*sum(population)/k)
    U = math.floor((1+population_deviation/2)*sum(population)/k)
    print("L =",L,"U =",U,"k =",k)
    m._row = [state, G.number_of_nodes(),G.number_of_edges(), k, L, U]
    maxp = max(population[i] for i in G.nodes)
    if maxp>U and base == 'hess':
        print("ERROR: infeasible instance. County with population",maxp,"which is larger than U =",U)
        m._row += ['-','-','-','-','-','-','-','-','-','-','-','-','-']
        return m._row
    elif k==1 and base == 'hess':
        print("ERROR: k=1. Nothing to solve for.")
        m._row += ['-','-','-','-','-','-','-','-','-','-','-','-','-']
        return m._row
    elif maxp>U and base == 'labeling':
        print("ERROR: infeasible instance. County with population",maxp,"which is larger than U =",U)
        m._row += ['-','-','-','-','-','-','-','-']
        return m._row
    elif k==1 and base == 'labeling':
        print("ERROR: k=1. Nothing to solve for.")
        m._row += ['-','-','-','-','-','-','-','-']
        return m._row

    ###########################
    # BUILD AND SOLVE MIP 
    ###########################  
    
    
    DG = nx.DiGraph(G) # bidirected version of G
    
    set_edge_lengths(G, weighted)
    
    
    #################
    # Apply heuristic
    #################
    if use_heuristic:
        heur_obj = heur["obj"]
        heur_districts = []
        for j in range(k):
            district =[]
            for node in heur["nodes"]:
                if node["district"] == j:
                    district.append(node["index"])
            heur_districts.append(district)        
        m._row.append(heur_obj)
        #heur_stop = time.time()
        m._row.append(heur["time"])
    else:
        heur_districts = []
        m._row.append("None")
        m._row.append("None")
    
    ###########################
    # Ordering
    ###########################    
    if OrderChoice!= "None":       
        order_start = time.time()  
        (ordered_vertices,S) = ordering.find_ordering(G, population, L, heur, k, OrderChoice, state, land_parcel)
        position = ordering.construct_position(ordered_vertices)
        order_stop = time.time()    
        
        m._S = S
        
        m._row.append(len(S))
        m._row.append(round(order_stop-order_start,2))
        #m._df['color'] = -1
        #m._df['label'] = -1
        #colors = []
        if S!=[]:
            #print("S=",S)
            m._df['S']= -1  # create a new column for S
            #for i in G.nodes:
             #   m._df['S'][i] = -1
            for i in G.nodes:
                #m._df['label'][i] = i
                if i in S:
                    geoID = G.node[i]["GEOID10"]
                    for u in G.nodes:
                        if geoID == m._df['GEOID10'][u]:
                            m._df['S'][u] = 0
                    #m._df['S'][i] = 0
                    #colors.append('white')
                else:
                    #print("non_S: ", i)
                    geoID = G.node[i]["GEOID10"]
                    for u in G.nodes:
                        if geoID == m._df['GEOID10'][u]:
                            m._df['S'][u] = 1
                    #m._df['S'][i] = 1
                    #colors.append('grey')
            cmap = LinearSegmentedColormap.from_list('S', [(0, 'white'), (1, 'lightgray')])    
            #colors = ['grey', 'white']  
            
            splot = m._df.plot(cmap=cmap, column='S',figsize=(10, 10), linewidth=1, edgecolor='0.25').get_figure()  # display the S map
            plt.axis('off')
            if land_parcel == "county":
                splot.savefig(state+"_S_county.png")
            elif land_parcel == "tract":
                splot.savefig(state+"_S_tract.png")
      
    else:
        m._row.append("None")
        m._row.append("None")
        ordered_vertices = None
        
        '''
        new_ordering = []
        
        for i in ordered_vertices:
            new_ordering.append(i+1)
        '''
            
        #print("ordered_vertices: ", new_ordering)
        
        #print("L:", L, ", U: ", U, ", Population of Bernalillo: ", G.node[18]["TOTPOP"])
        
        #for i in ordered_vertices:
         #   print("County: ", i+1, " with population: ", G.node[i]["TOTPOP"])
            
    
    ####################
    # All optima
    ####################
    # collect alternative optima?
    if m._all_optima==True:
        m.Params.PoolSearchMode=2 # find n best solutions
        m.Params.PoolSolutions=optima_limit # n = 1,000
        m.Params.PoolGap=0.1 # only find solutions within 0.1% of optimal (if opt<1000, this returns only optimal solutions)
    
    
    ############################
    # Base model
    ############################    
    if base == "hess":
        # X[i,j]=1 if vertex i is assigned to (district centered at) vertex j
        m._X = m.addVars(DG.nodes, DG.nodes, vtype=GRB.BINARY)
        row = main_hess.build_base_hess_model(m, G, population, L, U, k, heur_districts, ordered_vertices, position)
    elif base == "labeling":        
        m._X = m.addVars(DG.nodes, range(k), vtype=GRB.BINARY)
        if m._model != "cut" or m._model != "lcut":
            m._R = m.addVars(DG.nodes, range(k), vtype=GRB.BINARY)
        #add labeling constraints
        row = main_labeling.build_base_labeling_model(m, G, population, L, U, k, heur_districts, ordered_vertices, position)
    else:
        print("Error: please enter a valid base (i.e., either hess or labeling)")
   
    #############################
    # RETRIEVE AND DRAW SOLUTION(S)
    #############################
 
    return row
#############################
# RUN ALL STATES AND SAVE TO CSV/PNG
#############################
# write results to csv file. The csv filename includes the parameters used
fn = base+ "_cut_edges_" + OrderChoice #+ model + "_" + OrderChoice 
if weighted:
    fn += "_weighted"
if extended==True:
    fn += "_extended"
fn += ".csv"

with open(fn, 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    if base == 'hess':
        fields = ['State', 'n', 'm', 'k', 'L', 'U', 'heur', 'heur time', '|S|', 'Order time', 'DiagFixed', 'LFixed', 'UFixed', 'XFixed', 'Total X', 'ZFixed', 'Total Z', '#B&B', 'bound', 'obj', 'MIP time', 'Model']
    elif base == 'labeling':
        fields = ['State', 'n', 'm', 'k', 'L', 'U', 'heur', 'heur time', '|S|', 'Order time', 'LFixed', 'UFixed', 'X&R Fixed', 'Total X&R', 'ZFixed', 'Total Z', '#B&B', 'bound', 'obj', 'MIP time', 'Model']
    else:
        print("No correct base! Please enter a valid base (i.e., 'hess' or 'labeling')")
    csvwriter.writerow(fields)
    
    for state in state_codes.keys():
        #if state == 'UT' or state == 'MS' or state == 'AR' or state == 'NV':
        #if state == 'AL':
        #if state == 'ME' or state == 'LA':
        if state == 'NH' or state == 'ID' or state == 'ME' or state == 'WV' or state == 'NM' or state == 'NE':
        #if state == 'OK' or state == 'AL' or state == 'NE' or state == 'AR' or state == 'KS' or state == 'IA' or state == 'ID' or state == 'WV' or state == 'NM' or state == 'MS':
            if base == 'hess':
                for model in hess_models:
                    #row = build_cut_edges(state,base,model,weighted,extended,OrderChoice,all_optima,optima_limit,symmetry) + [model]
                    row = build_cut_edges(state,base,model,weighted,extended,OrderChoice,all_optima,optima_limit,symmetry) + [model]
                    print(row)
                    csvwriter.writerow(row)
            elif base == 'labeling':
                for model in labeling_models:
                    row = build_cut_edges(state,base,model,weighted,extended,OrderChoice,all_optima,optima_limit,symmetry) + [model]
                    print(row)
                    csvwriter.writerow(row)
                    #for i in G.nodes:
                     #   print()
            else:
                print("No correct base! Please enter a valid base (i.e., 'hess' or 'labeling')")
