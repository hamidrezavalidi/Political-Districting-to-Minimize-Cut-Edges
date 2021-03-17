###########################
# Imports
###########################  

import gurobipy as gp
from gurobipy import GRB 

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

import math
import networkx as nx
import csv
import time
import json
import sys

import hess
import labeling
import ordering
import fixing
import separation

from gerrychain import Graph
import geopandas as gpd

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

number_of_congressional_districts = {
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

default_config = {
    'state' : 'OK',
    'level' : 'county',
    'base' : 'labeling',
    'contiguity' : 'scf',
    'symmetry' : 'orbitope',
    'extended' : True,
    'order' : 'B_decreasing',
    'heuristic' : True
    }

available_config = {
    'state' : { key for key in state_codes.keys() },
    'level' : {'county', 'tract'},
    'base' : {'hess', 'labeling'},
    'contiguity' : {'none', 'lcut', 'scf', 'shir'},
    'symmetry' : {'default', 'aggressive', 'orbitope'},  # orbitope only for labeling
    'extended' : {True, False},
    'order' : {'none', 'decreasing', 'B_decreasing'},
    'heuristic' : {True, False}
    }


###########################
# Read configs/inputs and set parameters
###########################  

# read configs file and load into a Python dictionary
configs_file = open('config.json','r')
batch_configs = json.load(configs_file)
configs_file.close()


# run experiment for each config in batch_config file
for key in batch_configs.keys():
      
    # get config and check for errors
    config = batch_configs[key]
    print("In run",key,"using config",config,end='.')
    for ckey in config.keys():
        if config[ckey] not in available_config[ckey]:
            errormessage = "Error: the config option"+ckey+":"+config[ckey]+"is not known."
            sys.exit(errormessage)
    print("")
                   
    # read input data
    state = config['state']
    code = state_codes[state]
    level = config['level']
    G = Graph.from_json("data/"+level+"/dual_graphs/"+level+code+".json")
    DG = nx.DiGraph(G) # bidirected version of G
    df = gpd.read_file("data/"+level+"/shape_files/"+state+"_"+level+".shp")      

    # set parameters
    k = number_of_congressional_districts[state]    
    
    population = [G.nodes[i]['TOTPOP'] for i in G.nodes()]
    deviation = 0.01
    L = math.ceil((1-deviation/2)*sum(population)/k)
    U = math.floor((1+deviation/2)*sum(population)/k)
    print("L =",L,", U =",U,", k =",k)
    
    # abort early for trivial or overtly infeasible instances
    maxp = max(population[i] for i in G.nodes)
    if k==1 or maxp>U:
        print("k=",k,", max{ p_v | v in V } =",maxp,", U =",U,end='.')
        sys.exit("Aborting early, either due to trivial instance or overtly infeasible instance.")
           
    # read heuristic solution from external file (?)
    heuristic = config['heuristic']
    if heuristic:
        heuristic_file = open('data/'+level+"/heuristic/heur_"+state+"_"+level+".json",'r')
        heuristic_dict = json.load(heuristic_file)       
        heuristic_districts = [ [node['index'] for node in heuristic_dict['nodes'] if node['district']==j ] for j in range(k) ]
        heuristic_obj = heuristic_dict['obj']
        heuristic_time = heuristic_dict['time']
        heuristic_iterations = heuristic_dict['iterations']
    else:
        heuristic_districts = None
        heuristic_obj = None
        heuristic_time = None
        heuristic_iterations = None
        
           
    ############################
    # Build base model
    ############################   
    
    m = gp.Model()
    m._DG = DG
    base = config['base']
    
    if base == 'hess':
        # X[i,j]=1 if vertex i is assigned to (district centered at) vertex j
        m._X = m.addVars(DG.nodes, DG.nodes, vtype=GRB.BINARY)
        hess.add_base_constraints(m, population, L, U, k)
    
    if base == 'labeling':        
        # X[i,j]=1 if vertex i is assigned to district j in {0,1,2,...,k-1}
        m._X = m.addVars(DG.nodes, range(k), vtype=GRB.BINARY)
        if config['symmetry']=='orbitope' or config['contiguity'] in {'scf', 'shir'}:
            m._R = m.addVars(DG.nodes, range(k), vtype=GRB.BINARY)
        labeling.add_base_constraints(m, population, L, U, k)

                
    ############################################      
    # Add (extended?) objective 
    ############################################         
    
    extended = config['extended']
    
    if base == 'hess':
        if extended:
            hess.add_extended_objective(m, G)
        else:
            hess.add_objective(m, G)
               
    if base == 'labeling':
        if extended:
            labeling.add_extended_objective(m, G, k)
        else:
            labeling.add_objective(m, G, k)
            
    
    ####################################   
    # Contiguity constraints
    ####################################      
            
    contiguity = config['contiguity']
    m._callback = None
    m._population = population
    m._U = U
    m._k = k
    m._base = base
    
    if base == 'hess':
        if contiguity == 'shir':
            hess.add_shir_constraints(m, population, U)
        elif contiguity == 'scf':
            hess.add_scf_constraints(m, G, extended)
        elif contiguity == 'lcut':
            m.Params.lazyConstraints = 1
            m._callback = separation.lcut_separation_generic
         
    # FIXME -- NEED TO CHECK THESE FUNCTIONS:               
    if base == 'labeling':
        if contiguity == 'shir':
            labeling.add_shir_constraints(m, k, config['symmetry'])
        elif contiguity == 'scf':
            labeling.add_scf_constraints(m, G, k, extended, config['symmetry'])
        elif contiguity == 'lcut':
            m.Params.lazyConstraints = 1
            m._callback = separation.lcut_separation_generic 
         
    m.update()
    
    
    ############################################
    # Vertex ordering and max B problem 
    ############################################  
        
    order = config['order']
    
    if order == 'B_decreasing':
        B = ordering.solve_maxB_problem(DG, population, L, k, heuristic_districts)
    else:
        B = []
    
    my_ordering = ordering.find_ordering(order, B, DG, population)
    position = ordering.construct_position(my_ordering)
        

    ####################################   
    # Symmetry handling
    ####################################    
    
    symmetry = config['symmetry']
    
    if symmetry == 'aggressive':
        m.Params.symmetry = 2
    elif symmetry == 'orbitope':
        if base == 'labeling':
            labeling.add_orbitope_extended_formulation(m, G, k, my_ordering)
        else:
            sys.exit("Error: orbitope only available for labeling base model.")     
            
            
    ####################################   
    # Variable fixing
    ####################################    
    
    if base == 'hess':
        DFixings = fixing.do_Hess_DFixing(m, G, position)
        
        if contiguity == 'none':
            LFixings = fixing.do_Hess_LFixing_without_Contiguity(m, G, population, L, my_ordering)
            UFixings = fixing.do_Hess_UFixing_without_Contiguity(m, G, population, U)
        else:
            LFixings = fixing.do_Hess_LFixing(m, G, population, L, my_ordering)
            UFixings = fixing.do_Hess_UFixing(m, DG, population, U, my_ordering)         
        
        if extended:
            ZFixings = fixing.do_Hess_ZFixing(m, G)
        else:
            ZFixings = 0
    
    
    if base == 'labeling':
        DFixings = fixing.do_Labeling_DFixing(m, G, my_ordering, k)
        
        if contiguity == 'none':
            if symmetry == 'orbitope':
                LFixings = fixing.do_Labeling_LFixing_without_Contiguity(m, G, population, L, my_ordering, k)
            else:
                LFixings = 0
            (UFixings_X, UFixings_R) = fixing.do_labeling_UFixing_without_Contiguity()
        else:
            LFixings = fixing.do_Labeling_LFixing(m, G, population, L, my_ordering, k)
            (UFixings_X, UFixings_R) = fixing.do_Labeling_UFixing(m, DG, population, U, my_ordering, k)
        
        if extended:
            ZFixings = fixing.do_Labeling_ZFixing(m, G)
        else:
            ZFixings = 0
    
    
    ####################################   
    # Inject heuristic warm start
    ####################################    
    
    if heuristic and base == 'hess':
        for district in heuristic_districts:    
            p = min([position[v] for v in district])
            j = my_ordering[p]
            for i in district:
                m._X[i,j].start = 1
                    
    if heuristic and base == 'labeling':
        center_positions = [ min( position[v] for v in heuristic_districts[j] ) for j in range(k) ] 
        cplabel = { center_positions[j] : j for j in range(k) }
    
        # what node r will root the new district j? The one with earliest position.
        for j in range(k):
            min_cp = min(center_positions)
            r = my_ordering[min_cp]
            old_j = cplabel[min_cp]
            
            for i in heuristic_districts[old_j]:
                m._X[i,j].start = 1
                
            center_positions.remove(min_cp)
                
    
    ####################################   
    # Solve MIP
    ####################################  
    
    m.optimize(m._callback)
    
    
    # FIXME : add map pictures and max B pictures
    
    # FIXME : add timers
    
    # FIXME : report runtime results to file
    
    # FIXME : print districting solution to file
    
    
        
    # ###########################
    # # Ordering
    # ###########################    
    # if OrderChoice!= "None":       
    #     order_start = time.time()  
    #     (ordered_vertices,S) = ordering.find_ordering(G, population, L, heur, k, OrderChoice, state, land_parcel)
    #     position = ordering.construct_position(ordered_vertices)
    #     order_stop = time.time()    
        
    #     m._S = S
        
    #     m._row.append(len(S))
    #     m._row.append(round(order_stop-order_start,2))
    #     #m._df['color'] = -1
    #     #m._df['label'] = -1
    #     #colors = []
    #     if S!=[]:
    #         #print("S=",S)
    #         m._df['S']= -1  # create a new column for S
    #         #for i in G.nodes:
    #          #   m._df['S'][i] = -1
    #         for i in G.nodes:
    #             #m._df['label'][i] = i
    #             if i in S:
    #                 geoID = G.node[i]["GEOID10"]
    #                 for u in G.nodes:
    #                     if geoID == m._df['GEOID10'][u]:
    #                         m._df['S'][u] = 0
    #                 #m._df['S'][i] = 0
    #                 #colors.append('white')
    #             else:
    #                 #print("non_S: ", i)
    #                 geoID = G.node[i]["GEOID10"]
    #                 for u in G.nodes:
    #                     if geoID == m._df['GEOID10'][u]:
    #                         m._df['S'][u] = 1
    #                 #m._df['S'][i] = 1
    #                 #colors.append('grey')
    #         cmap = LinearSegmentedColormap.from_list('S', [(0, 'white'), (1, 'lightgray')])    
    #         #colors = ['grey', 'white']  
            
    #         splot = m._df.plot(cmap=cmap, column='S',figsize=(10, 10), linewidth=1, edgecolor='0.25').get_figure()  # display the S map
    #         plt.axis('off')
    #         if land_parcel == "county":
    #             splot.savefig(state+"_S_county.png")
    #         elif land_parcel == "tract":
    #             splot.savefig(state+"_S_tract.png")
      

