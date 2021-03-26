###########################
# Imports
###########################  

import gurobipy as gp
from gurobipy import GRB 

import matplotlib.pyplot as plt

from datetime import date
import math
import networkx as nx
import csv
import time
import json
import sys
import os

import hess
import labeling
import ordering
import fixing
import separation

from gerrychain import Graph
import geopandas as gpd


################################################
# Summarize computational results to csv file
################################################ 

from csv import DictWriter
def append_dict_as_row(file_name, dict_of_elem, field_names):
    # Open file in append mode
    with open(file_name, 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        dict_writer = DictWriter(write_obj, fieldnames=field_names)
        # Add dictionary as wor in the csv
        dict_writer.writerow(dict_of_elem)
        
        
################################################
# Writes districting solution to json file
################################################ 

def export_to_json(G, districts, filename):
    with open(filename, 'w') as outfile:
        soln = {}
        soln['nodes'] = []
        for j in range(len(districts)):
            for i in districts[j]:
                soln['nodes'].append({
                        'name': G.nodes[i]["NAME10"],
                        'index': i,
                        'district': j
                        })
        json.dump(soln, outfile)
        
        
################################################
# Draws districts and saves to png file
################################################ 

def export_to_png(G, df, districts, filename):
    df['assignment'] = -1
    for j in range(len(districts)):
        for i in districts[j]:
            geoID = G.nodes[i]["GEOID10"]
            for u in G.nodes:
                if geoID == df['GEOID10'][u]:
                    df['assignment'][u] = j
                    
    my_fig = df.plot(column='assignment').get_figure()
    RESIZE_FACTOR = 3
    my_fig.set_size_inches(my_fig.get_size_inches()*RESIZE_FACTOR)
    plt.axis('off')
    my_fig.savefig(filename)


################################################
# Draws max B set and saves to png file
################################################ 

def export_B_to_png(G, df, B, filename):
    
    df['B'] = -1
    for i in G.nodes:
        df['B'][i] = 0
        
    for i in B:
        geoID = G.nodes[i]["GEOID10"]
        for u in G.nodes:
            if geoID == df['GEOID10'][u]:
                df['B'][u] = 1
        
    my_fig = df.plot(column='B').get_figure()
    RESIZE_FACTOR = 3
    my_fig.set_size_inches(my_fig.get_size_inches()*RESIZE_FACTOR)
    plt.axis('off')
    my_fig.savefig(filename)
    

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
    'fixing' : True,
    'contiguity' : 'scf',
    'symmetry' : 'orbitope',
    'extended' : True,
    'order' : 'B_decreasing',
    'heuristic' : True,
    'lp': True
}

available_config = {
    'state' : { key for key in state_codes.keys() },
    'level' : {'county', 'tract'},
    'base' : {'hess', 'labeling'},
    'fixing' : {True, False},
    'contiguity' : {'none', 'lcut', 'scf', 'shir'},
    'symmetry' : {'default', 'aggressive', 'orbitope'},  # orbitope only for labeling
    'extended' : {True, False},
    'order' : {'none', 'decreasing', 'B_decreasing'},
    'heuristic' : {True, False},
    'lp' : {True, False} # solve and report root LP bound? (in addition to MIP)
}


###############################################
# Read configs/inputs and set parameters
############################################### 

# read configs file and load into a Python dictionary

if len(sys.argv)>1:
    # name your own config file in command line, like this: 
    #       python main.py usethisconfig.json
    # to keep logs of the experiments, redirect to file, like this:
    #       python main.py usethisconfig.json 1>>log_file.txt 2>>error_file.txt
    config_filename = sys.argv[1] 
else:
    config_filename = 'config.json' # default
    
print("Reading config from",config_filename)    
config_filename_wo_extension = config_filename.rsplit('.',1)[0]
configs_file = open(config_filename,'r')
batch_configs = json.load(configs_file)
configs_file.close()

# create directory for results
path = os.path.join("..", "results_for_" + config_filename_wo_extension) 
os.mkdir(path) 

# print results to csv file
today = date.today()
today_string = today.strftime("%Y_%b_%d") # Year_Month_Day, like 2019_Sept_16
results_filename = "../results_for_" + config_filename_wo_extension + "/results_" + config_filename_wo_extension + "_" + today_string + ".csv" 

# prepare csv file by writing column headers
with open(results_filename,'w',newline='') as csvfile:   
    my_fieldnames = ['run','state','level','base','fixing','contiguity','symmetry','extended','order','heuristic','lp'] # configs
    my_fieldnames += ['k','L','U','n','m'] # params
    my_fieldnames += ['heur_obj', 'heur_time', 'heur_iter'] # heuristic info
    my_fieldnames += ['B_q', 'B_size', 'B_time', 'B_timelimit'] # max B info
    my_fieldnames += ['DFixings', 'LFixings', 'UFixings_X', 'UFixings_R', 'ZFixings'] # fixing info
    my_fieldnames += ['LP_obj', 'LP_time'] # root LP info
    my_fieldnames += ['MIP_obj','MIP_bound','MIP_time', 'MIP_timelimit', 'MIP_status', 'MIP_nodes', 'connected'] # MIP info
    writer = csv.DictWriter(csvfile, fieldnames = my_fieldnames)
    writer.writeheader()
    
############################################################
# Run experiments for each config in batch_config file
############################################################

for key in batch_configs.keys(): 
      
    # get config and check for errors
    config = batch_configs[key]
    print("In run",key,"using config:",config,end='.')
    for ckey in config.keys():
        if config[ckey] not in available_config[ckey]:
            errormessage = "Error: the config option"+ckey+":"+config[ckey]+"is not known."
            sys.exit(errormessage)
    print("")
    
    # fill-in unspecified configs using default values
    for ckey in available_config.keys():
        if ckey not in config.keys():
            print("Using default value",ckey,"=",default_config[ckey],"since no option was selected.")
            config[ckey] = default_config[ckey]
        
    # initialize dictionary to store this run's results
    result = config
    result['run'] = key            
                   
    # read input data
    state = config['state']
    code = state_codes[state]
    level = config['level']
    G = Graph.from_json("../data/"+level+"/dual_graphs/"+level+code+".json")
    DG = nx.DiGraph(G) # bidirected version of G
    df = gpd.read_file("../data/"+level+"/shape_files/"+state+"_"+level+".shp")      

    # set parameters
    k = number_of_congressional_districts[state]        
    population = [G.nodes[i]['TOTPOP'] for i in G.nodes()]    
    deviation = 0.01
    L = math.ceil((1-deviation/2)*sum(population)/k)
    U = math.floor((1+deviation/2)*sum(population)/k)
    print("L =",L,", U =",U,", k =",k)
    result['k'] = k
    result['L'] = L
    result['U'] = U
    result['n'] = G.number_of_nodes()
    result['m'] = G.number_of_edges()
    
    # abort early for trivial or overtly infeasible instances
    maxp = max(population[i] for i in G.nodes)
    if k==1 or maxp>U:
        print("k=",k,", max{ p_v | v in V } =",maxp,", U =",U,end='.')
        sys.exit("Aborting early, either due to trivial instance or overtly infeasible instance.")
           
    # read heuristic solution from external file (?)
    heuristic = config['heuristic']
    if heuristic:
        heuristic_file = open('../data/'+level+"/heuristic/heur_"+state+"_"+level+".json",'r')
        heuristic_dict = json.load(heuristic_file)       
        heuristic_districts = [ [node['index'] for node in heuristic_dict['nodes'] if node['district']==j ] for j in range(k) ]
        result['heur_obj'] = heuristic_dict['obj']
        result['heur_time'] = heuristic_dict['time']
        result['heur_iter'] = heuristic_dict['iterations']
    else:
        heuristic_districts = None
        result['heur_obj'] = 'n/a'
        result['heur_time'] = 'n/a'
        result['heur_iter'] = 'n/a'
        
           
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
        if config['symmetry']=='orbitope' or config['contiguity'] in {'scf', 'shir', 'lcut'}:
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
            hess.add_shir_constraints(m)
        elif contiguity == 'scf':
            hess.add_scf_constraints(m, G, extended)
        elif contiguity == 'lcut':
            m.Params.lazyConstraints = 1
            m._callback = separation.lcut_separation_generic
                    
    if base == 'labeling':
        if contiguity == 'shir':
            labeling.add_shir_constraints(m, config['symmetry'])
        elif contiguity == 'scf':
            labeling.add_scf_constraints(m, G, extended, config['symmetry'])
        elif contiguity == 'lcut':
            m.Params.lazyConstraints = 1
            m._callback = separation.lcut_separation_generic 
         
    m.update()
    
    
    ############################################
    # Vertex ordering and max B problem 
    ############################################  
        
    order = config['order']
    
    if order == 'B_decreasing':
        (B, result['B_q'], result['B_time'], result['B_timelimit']) = ordering.solve_maxB_problem(DG, population, L, k, heuristic_districts)
        
        # draw set B on map and save
        fn_B = "../" + "results_for_" + config_filename_wo_extension + "/" + result['state'] + "-" + result['level'] + "-maxB.png"       
        export_B_to_png(G, df, B, fn_B)
    else:
        (B, result['B_q'], result['B_time'], result['B_timelimit']) = (list(),'n/a','n/a', 'n/a')
        
    result['B_size'] = len(B)
    
    vertex_ordering = ordering.find_ordering(order, B, DG, population)
    position = ordering.construct_position(vertex_ordering)
    
    print("Vertex ordering =", vertex_ordering)  
    print("Position vector =", position)
    print("Set B =", B)

    ####################################   
    # Symmetry handling
    ####################################    
    
    symmetry = config['symmetry']
    
    if symmetry == 'aggressive':
        m.Params.symmetry = 2
    elif symmetry == 'orbitope':
        if base == 'labeling':
            labeling.add_orbitope_extended_formulation(m, G, k, vertex_ordering)
        else:
            sys.exit("Error: orbitope only available for labeling base model.")     
            
            
    ####################################   
    # Variable fixing
    ####################################    
    
    do_fixing = config['fixing']
    
    if do_fixing and base == 'hess':
        result['DFixings'] = fixing.do_Hess_DFixing(m, G, position)
        result['UFixings_R'] = 'n/a'
        
        if contiguity == 'none':
            result['LFixings'] = fixing.do_Hess_LFixing_without_Contiguity(m, G, population, L, vertex_ordering)
            result['UFixings_X'] = fixing.do_Hess_UFixing_without_Contiguity(m, G, population, U)
        else:
            result['LFixings'] = fixing.do_Hess_LFixing(m, G, population, L, vertex_ordering)
            result['UFixings_X'] = fixing.do_Hess_UFixing(m, DG, population, U, vertex_ordering)         
        
        if extended:
            result['ZFixings'] = fixing.do_Hess_ZFixing(m, G)
        else:
            result['ZFixings'] = 0
                
    
    if do_fixing and base == 'labeling':
        result['DFixings'] = fixing.do_Labeling_DFixing(m, G, vertex_ordering, k)
        
        if contiguity == 'none':
            if symmetry == 'orbitope':
                result['LFixings'] = fixing.do_Labeling_LFixing_without_Contiguity(m, G, population, L, vertex_ordering, k)
            else:
                result['LFixings'] = 0
            (result['UFixings_X'], result['UFixings_R']) = fixing.do_labeling_UFixing_without_Contiguity()
        else:
            result['LFixings'] = fixing.do_Labeling_LFixing(m, G, population, L, vertex_ordering, k)
            (result['UFixings_X'], result['UFixings_R']) = fixing.do_Labeling_UFixing(m, DG, population, U, vertex_ordering, k)
        
        if extended:
            result['ZFixings'] = fixing.do_Labeling_ZFixing(m, G, k)
        else:
            result['ZFixings'] = 0
            
    if not do_fixing:
        result['DFixings'] = 0
        result['UFixings_R'] = 0
        result['LFixings'] = 0
        result['UFixings_X'] = 0
        result['ZFixings'] = 0
            
    
    ######################################################################################
    # Solve root LP? Used only for reporting purposes. Not used for MIP solve.
    ######################################################################################  
    
    if config['lp']:
        r = m.relax() # LP relaxation of MIP model m
        #r.Params.LogToConsole = 0 # keep log to a minimum
        r.Params.Method = 3 # use concurrent LP solver
        r.Params.TimeLimit = 3600 # one-hour time limit for solving LP
        print("To get the root LP bound, now solving a (separate) LP model.")
        
        lp_start = time.time()
        r.optimize()
        lp_end = time.time()
        
        if r.status == GRB.OPTIMAL:
            result['LP_obj'] = '{0:.2f}'.format(r.objVal)
        elif r.status == GRB.TIME_LIMIT:
            result['LP_obj'] = 'TL'
        else:
            result['LP_obj'] = '?'
        result['LP_time'] = '{0:.2f}'.format(lp_end - lp_start)
        
    else:
        result['LP_obj'] = 'n/a'
        result['LP_time'] = 'n/a'
        
    
    ####################################   
    # Inject heuristic warm start
    ####################################    
    
    if heuristic and base == 'hess':
        for district in heuristic_districts:    
            p = min([position[v] for v in district])
            j = vertex_ordering[p]
            for i in district:
                m._X[i,j].start = 1
                    
    if heuristic and base == 'labeling':
        center_positions = [ min( position[v] for v in heuristic_districts[j] ) for j in range(k) ] 
        cplabel = { center_positions[j] : j for j in range(k) }
    
        # what node r will root the new district j? The one with earliest position.
        for j in range(k):
            min_cp = min(center_positions)
            r = vertex_ordering[min_cp]
            old_j = cplabel[min_cp]
            
            for i in heuristic_districts[old_j]:
                m._X[i,j].start = 1
                
            center_positions.remove(min_cp)
                
    
    ####################################   
    # Solve MIP
    ####################################  
    
    result['MIP_timelimit'] = 3600 # set a one hour time limit
    m.Params.TimeLimit = result['MIP_timelimit']
    m.Params.Method = 3 # use concurrent method for root LP. Useful for degenerate models
    
    start = time.time()
    m.optimize(m._callback)
    end = time.time()
    result['MIP_time'] = '{0:.2f}'.format(end-start)
    
    result['MIP_status'] = int(m.status)
    result['MIP_nodes'] = int(m.NodeCount)
    result['MIP_bound'] = m.objBound
    
    # report best solution found
    if m.SolCount > 0:
        result['MIP_obj'] = int(m.objVal)

        if base == 'hess':
            labels = [ j for j in DG.nodes if m._X[j,j].x > 0.5 ]
        else: # base == 'labeling'
            labels = [ j for j in range(k) ]
            
        districts = [ [ i for i in DG.nodes if m._X[i,j].x > 0.5 ] for j in labels]
        print("best solution (found) =",districts)
        
        fn = "../" + "results_for_" + config_filename_wo_extension + "/" + result['state'] + "-" + result['level'] + "-" + result['base'] + "-" + result['contiguity']
        
        # export solution to .json file
        json_fn = fn + ".json"
        export_to_json(G, districts, json_fn)
        
        # export solution to .png file (districting map)
        png_fn = fn + ".png"
        export_to_png(G, df, districts, png_fn)
        
        # is solution connected?
        connected = True
        for district in districts:
            if not nx.is_connected(G.subgraph(district)):
                connected = False
        result['connected'] = connected
        
    else:
        result['MIP_obj'] = 'no_solution_found'
        result['connected'] = 'n/a'
        
            
    ####################################   
    # Summarize results of this run to csv file
    ####################################  
    
    append_dict_as_row(results_filename,result,my_fieldnames)
    
