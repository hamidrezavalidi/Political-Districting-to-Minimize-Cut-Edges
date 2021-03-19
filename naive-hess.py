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
    'contiguity' : True
}

available_config = {
    'state' : { key for key in state_codes.keys() },
    'level' : {'county', 'tract'},
    'contiguity' : {True, False},
}


###############################################
# Read configs/inputs and set parameters
############################################### 

# read configs file and load into a Python dictionary
configs_file = open('naive-hess-config.json','r')
batch_configs = json.load(configs_file)
configs_file.close()

# print results to csv file
today = date.today()
today_string = today.strftime("%Y_%b_%d") # Year_Month_Day, like 2019_Sept_16
results_filename = "naive-hess-results_" + today_string + ".csv" 

# prepare csv file by writing column headers
with open(results_filename,'w',newline='') as csvfile:   
    my_fieldnames = ['run','state','level','contiguity'] # configs
    my_fieldnames += ['k','L','U','n','m'] # params
    my_fieldnames += ['MIP_obj','MIP_bound','MIP_time', 'MIP_timelimit', 'MIP_status', 'MIP_nodes'] # MIP info
    writer = csv.DictWriter(csvfile, fieldnames = my_fieldnames)
    writer.writeheader()
    
############################################################
# Run experiments for each config in batch_config file
############################################################

for key in batch_configs.keys(): 
      
    # get config and check for errors
    config = batch_configs[key]
    print("In run",key,"using config",config,end='.')
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
        
           
    ############################
    # Build base MIP
    ############################   
    
    m = gp.Model()
         
    # X[i,j]=1 if vertex i is assigned to district j in {0,1,2,...,k-1}
    X = m.addVars(DG.nodes, DG.nodes, vtype=GRB.BINARY)
    
    # Y[u,v]=1 if edge {u,v} is cut
    Y = m.addVars(G.edges, vtype=GRB.BINARY)
    
    # objective is to minimize the number of cut edges
    m.setObjective( gp.quicksum(Y), GRB.MINIMIZE )
        
    # If i is assigned to vertex v, but its neighbor j is not, then edge ij is cut
    m.addConstrs( X[i,v] - X[j,v] <= Y[i,j] for i,j in G.edges for v in G.nodes)
    
    # Each vertex i assigned to one district
    m.addConstrs(gp.quicksum(X[i,j] for j in DG.nodes) == 1 for i in DG.nodes)
     
    # Pick k centers
    m.addConstr(gp.quicksum(X[j,j] for j in DG.nodes) == k)
    
    # Population balance: population assigned to vertex j should be in [L,U], if j is a center
    m.addConstrs(gp.quicksum(population[i] * X[i,j] for i in DG.nodes) <= U * X[j,j] for j in DG.nodes)
    m.addConstrs(gp.quicksum(population[i] * X[i,j] for i in DG.nodes) >= L * X[j,j] for j in DG.nodes)
    
    # Add coupling inequalities for added model strength
    m.addConstrs(X[i,j] <= X[j,j] for i in DG.nodes for j in DG.nodes)
        
    # Set branch priority on center vars
    # for j in DG.nodes:
    #     X[j,j].BranchPriority=1      
    
    # Do population-based diagonal fixing
    # for i in DG.nodes:
    #     for j in DG.nodes:
    #         if population[i] > population[j]:
    #             X[i,j].UB = 0
    
    ###################################################################
    # Add contiguity constraints? 
    # Based on flow model of Shirabe, see 
    #   http://www.optimization-online.org/DB_HTML/2020/01/7582.html
    ###################################################################

    if config['contiguity']:
        
        # F[j,u,v] tells how much flow (from source j) is sent across arc (u,v)
        F = m.addVars( DG.nodes, DG.edges, vtype=GRB.CONTINUOUS)
        
        # big-M
        M = G.number_of_nodes()
        
        # Flow of type j cannot enter vertex j. (It can only leave.)
        m.addConstrs( gp.quicksum(F[j,u,j] for u in DG.neighbors(j)) == 0 for j in DG.nodes)
        
        # If vertex i is assigned to vertex j, then i should receive one unit of type j flow
        m.addConstrs( gp.quicksum( F[j,u,i]-F[j,i,u] for u in DG.neighbors(i) ) == X[i,j] for i in DG.nodes for j in DG.nodes if i!=j)
        
        # If vertex i is not assigned to vertex j, then no type j flow can enter it
        m.addConstrs( gp.quicksum( F[j,u,i] for u in DG.neighbors(i) ) <= M*X[i,j] for i in DG.nodes for j in DG.nodes if i!=j)
                

    ####################################   
    # Solve MIP
    ####################################  
    
    result['MIP_timelimit'] = 3600 # set a one hour time limit
    m.params.TimeLimit = result['MIP_timelimit']
    
    # m.Params.Method = 3 # use concurrent method for root LP. Useful for degenerate models
    # m.Params.Symmetry = 2 # use aggressive symmetry handling.
    
    start = time.time()
    m.optimize()
    end = time.time()
    result['MIP_time'] = '{0:.2f}'.format(end-start)
    
    result['MIP_status'] = int(m.status)
    result['MIP_nodes'] = int(m.NodeCount)
    result['MIP_bound'] = int(m.objBound)
    
    # report best solution found
    if m.SolCount > 0:
        result['MIP_obj'] = int(m.objVal)

        centers = [ j for j in DG.nodes if X[j,j].x > 0.5 ]
        districts = [ [ i for i in DG.nodes if X[i,j].x > 0.5 ] for j in centers]
        print("best solution (found) =",districts)
        
        if result['contiguity']:
            contiguity_status = 'contiguity_imposed'
        else:
            contiguity_status = 'contiguity_not_imposed'
            
        fn = result['state'] + "-" + result['level'] + "-naivehess-" + contiguity_status
        
        # export solution to .json file
        json_fn = fn + ".json"
        export_to_json(G, districts, json_fn)
        
        # export solution to .png file (districting map)
        png_fn = fn + ".png"
        export_to_png(G, df, districts, png_fn)
    else:
        result['MIP_obj'] = 'no_solution_found'
        
            
    ####################################   
    # Summarize results of this run to csv file
    ####################################  
    
    append_dict_as_row(results_filename,result,my_fieldnames)
    
