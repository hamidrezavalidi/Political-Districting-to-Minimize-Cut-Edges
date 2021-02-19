import gurobipy as gp
from gurobipy import GRB 

#import networkx as nx

def construct_position(ordering):
    position = [-1 for i in range(len(ordering))]
    for p in range(len(ordering)):
        v = ordering[p]
        position[v] = p
    return position

def find_ordering(G, population, L, heur, k, OrderChoice, state, land_parcel):
    S = []
    if OrderChoice=="none":
        return ([v for v in G.nodes],S)
    elif OrderChoice=="increasing_population":
        nodes_with_population = [(i,population[i]) for i in G.nodes]
        nodes_with_population.sort(key=sort_by_second)
        return ([v for (v,p) in nodes_with_population],S)
    elif OrderChoice=="decreasing_population":
        nodes_with_population = [(i,population[i]) for i in G.nodes]
        nodes_with_population.sort(key=sort_by_second,reverse=True)
        decreasing_order = [v for (v,p) in nodes_with_population]
        #print("decreasing order: ", decreasing_order)
        nodes_with_population.sort(key=sort_by_second)
        increasing_order = [v for (v,p) in nodes_with_population]
        #print("increasing order: ", increasing_order)
        under_population = 0
        for i in increasing_order:
            under_population += population[i]
            if under_population <= L:
                S.append(i)
        return (decreasing_order,S)
    elif OrderChoice=="S":
        S = solve_maxS_recursive_model(G,population,L,state)
        V_S = [i for i in G.nodes if i not in S]
        return (V_S + S,S)
    elif OrderChoice=="S_decreasing":
        #warm = solve_heur_LFix(G, population, L, state, heur, k)
        #S = solve_maxS_recursive_model(G,population,L,state)
        S = solve_maxS_recursive_model_tract(G,population,L,k,state,heur,land_parcel)
        #if land_parcel == "county":
            #S = solve_maxS_recursive_model_county(G,population,L,k,state)
         #   S = solve_maxS_recursive_model_tract(G,population,L,k,state,land_parcel)
        #elif land_parcel == "tract":
         #   S = solve_maxS_recursive_model_tract(G,population,L,k,state,land_parcel)
        #else:
         #   print("Please enter a correct land parcel!")
        V_S_with_population = [(i,population[i]) for i in G.nodes if i not in S]
        V_S_with_population.sort(key=sort_by_second,reverse=True)
        return ([v for (v,p) in V_S_with_population] + S,S)
    elif OrderChoice=="heur_ordering":
        S = solve_heur_LFix(G, population, L, state, heur, k)
        V_S_with_population = [(i,population[i]) for i in G.nodes if i not in S]
        V_S_with_population.sort(key=sort_by_second,reverse=True)
        return ([v for (v,p) in V_S_with_population] + S,S)
    elif OrderChoice=="S_indep_decreasing":
        S = solve_max_indep(G,population,L,state)
        V_S_with_population = [(i,population[i]) for i in G.nodes if i not in S]
        V_S_with_population.sort(key=sort_by_second,reverse=True)
        return ([v for (v,p) in V_S_with_population] + S,S)
    else:
        print("ERROR: OrderChoice=",OrderChoice,"not supported.")
        
def sort_by_second(val):
    return val[1]   
        
def solve_maxS_recursive_model_tract(G,population,L,k,state,heur,land_parcel):
    m = gp.Model()
   
    # Y[i,j]=1 if vertices i and j belong to same component of G[S]
    # Y[i,i]=1 if i is selected in S
    q = k
    X = m.addVars(G.nodes, range(q), vtype=GRB.BINARY)
    
    #R = m.addVars(range(q), vtype=GRB.BINARY)
    
    S = m.addVars(G.nodes, vtype=GRB.CONTINUOUS)
    
    for i in G.nodes:
        S[i].ub = 1
        
    #Austin's warm_start
    for node in heur["nodes"]:  
        i = node["index"]
        for j in range(k):
            if j != node["district"]:
                X[i,j].start = 0.0
                
    
    '''    
    nodes_in_cut_edges = []    
    for node1 in heur["nodes"]:
        i = node1["index"]
        for node2 in heur["nodes"]:
            j = node2["index"]
            if (i,j) in G.edges and node1["district"] != node2["district"]:
                S[i].start = 0.0
                S[j].start = 0.0
                nodes_in_cut_edges.append(i)
                nodes_in_cut_edges.append(j)
                
    for node in heur["nodes"]:
        i = node["index"]
        if i not in nodes_in_cut_edges:
            X[i,node["district"]].start = 1.0
    '''        
    '''
    for i in G.nodes:
        for j in G.nodes:
            if i < j:
                X[i,j].ub = 0
    '''
    # assignment constraints            
    m.addConstrs( gp.quicksum(X[i,j] for j in range(q)) == S[i] for i in G.nodes)
                
    # the population reachable from j is less than L
    m.addConstrs( gp.quicksum(population[i] * X[i,j] for i in G.nodes) <= L-1 for j in range(q))
    
    m.addConstrs( X[u,j] + S[v] <= 1 + X[v,j] for u,v in G.edges for j in range(q))
    
    
    # Symmetry breaking tricks
    
    #for j in range(q-1):
     #   m.addConstr( gp.quicksum(population[i] * X[i,j] for i in G.nodes) >= gp.quicksum(population[i] * X[i,j+1] for i in G.nodes))
    
    '''    
    for i in G.nodes:
        for j in range(q):
            if j>i:
                X[i,j].ub = 0
    '''
    # Set objective function and solve MIP
    m.setObjective(gp.quicksum(S[i] for i in G.nodes), GRB.MAXIMIZE)
    #m.Params.Cuts=0 # turn off cuts
    m.Params.MIPFocus=1 # turn on MIPFocus
    #m.Params.Symmetry=2 # set symmetry parameter   
    m.Params.timeLimit=60 # 100-second time limit
    #for i in range(q):
     #   R[i].branchPriority=1
        
    #m.params.OutputFlag = 0

    m.params.LogFile='max_S_'+state+'_'+land_parcel+'.log'
    m.params.LogToConsole = 0    
    
    m.optimize()
    
    if m.status == GRB.OPTIMAL or m.status==GRB.TIME_LIMIT:
        #for i in G.nodes:
         #   if S[i].x > 0.5:
          #      print("vertex ", i, " with S value ", S[i].x)
        return [i for i in G.nodes if S[i].x > 0.5 ]
    else:
        return "Model status != optimal"


def solve_maxS_recursive_model_county(G,population,L,k,state):
    m = gp.Model()
   
    # Y[i,j]=1 if vertices i and j belong to same component of G[S]
    # Y[i,i]=1 if i is selected in S
    X = m.addVars(G.nodes, G.nodes, vtype=GRB.BINARY)
    
    S = m.addVars(G.nodes, vtype=GRB.CONTINUOUS)
    
    for i in G.nodes:
        S[i].ub = 1
    
    for i in G.nodes:
        for j in G.nodes:
            if i < j:
                X[i,j].ub = 0
                
                
    # assignment constraints            
    m.addConstrs( gp.quicksum(X[i,j] for j in G.nodes) == S[i] for i in G.nodes)
    
    
    #m.addConstrs(X[i,j] <= X[j,j] for i in G.nodes for j in G.nodes)
          
              
    # coupling constraints
    couplingConstr = m.addConstrs(X[i,j] <= X[j,j] for i in G.nodes for j in G.nodes)
    
    # Madke them user cuts
    for i in G.nodes:
        for j in G.nodes:
            couplingConstr[i,j].Lazy = -1
    
    #m.addConstrs(Y[i,j] <= Y[j,j] for i in G.nodes for j in G.nodes)
    
    # the population reachable from j is less than L
    m.addConstrs( gp.quicksum(population[i] * X[i,j] for i in G.nodes) <= (L-1)*X[j,j] for j in G.nodes)
    
    # if G[S] contains an i,v-path and (v,j) is an edge and j is selected, 
    #    then G[S] contains an i,j path
    #for v in G.nodes:
    m.addConstrs( X[u,j] + S[v] <= 1 + X[v,j] for u,v in G.edges for j in G.nodes)
    
    # indep set constraints
    m.addConstrs( X[u,u] + X[v,v] <= 1 for u,v in G.edges)
    
    #m.addConstrs( X[j,j] <= S[j] for j in G.nodes)
    # number of centers
    m.addConstr( gp.quicksum(X[j,j] for j in G.nodes) <= k)
    
    '''
    # fix far assignments to zero
    for (i,j) in G.edges:
        G[i][j]['weight'] = population[j]
    
    for j in G.nodes:
        dist = nx.shortest_path_length(G,source=j,weight='weight')
        for i in G.nodes:
            if dist[i]+population[j] >= L and i>j: #position[i]>position[j]:
                #print("UFixed: ", (i+1,j+1))
                X[i,j].UB=0
                #UFixed += 1
    '''
    # impose connectivity
    
    #for v in G.nodes:
     #   m.addConstrs( X[i,v] + S[j] <= 1 + X[j,v] for i in G.nodes for j in G.neighbors(v))
    
    
    #for i in G.nodes:
     #   for j in G.nodes:
      #      if i > j:
       #         m.addConstr( X[i,j] <= gp.quicksum(X[c,j] for c in G.neighbors(j)) )
               
    
    
    # Set objective function and solve MIP
    m.setObjective(gp.quicksum(S[i] for i in G.nodes), GRB.MAXIMIZE)
    m.Params.Cuts=0 # turn off cuts
    m.Params.timeLimit=60 # 100-second time limit
    for i in G.nodes:
        X[i,i].branchPriority=1
        
    #m.params.OutputFlag = 0

    m.params.LogFile='max_S_'+state+'_county.log'
    m.params.LogToConsole = 0    
    
    m.optimize()
    
    
    
    if m.status == GRB.OPTIMAL or m.status==GRB.TIME_LIMIT:
        #for i in G.nodes:
         #   if S[i].x > 0.5:
          #      print("vertex ", i, " with S value ", S[i].x)
        return [i for i in G.nodes if S[i].x > 0.5 ]
    else:
        return "Model status != optimal"

    
def solve_maxS_recursive_model(G,population,L,state):
    m = gp.Model()
   
    # Y[i,j]=1 if vertices i and j belong to same component of G[S]
    # Y[i,i]=1 if i is selected in S
    Y = m.addVars(G.nodes, G.nodes, vtype=GRB.BINARY)
    
    # warm start
    #for i in warm:
     #   Y[i,i].start = 1.0
    
    m.addConstrs(Y[i,j] == Y[j,i] for i in G.nodes for j in G.nodes)
    
    # coupling constraints
    couplingConstr = m.addConstrs(Y[i,j] <= Y[j,j] for i in G.nodes for j in G.nodes)
    
    # Madke them user cuts
    for i in G.nodes:
        for j in G.nodes:
            couplingConstr[i,j].Lazy = -1
    
    #m.addConstrs(Y[i,j] <= Y[j,j] for i in G.nodes for j in G.nodes)
    
    # the population reachable from j is less than L
    m.addConstrs( gp.quicksum(population[i] * Y[i,j] for i in G.nodes) <= (L-1)*Y[j,j] for j in G.nodes)
    
    # if G[S] contains an i,v-path and (v,j) is an edge and j is selected, 
    #    then G[S] contains an i,j path
    for v in G.nodes:
        m.addConstrs( Y[i,v] + Y[j,j] <= 1 + Y[i,j] for i in G.nodes for j in G.neighbors(v))
    
    # Set objective function and solve MIP
    m.setObjective(gp.quicksum(Y[i,i] for i in G.nodes), GRB.MAXIMIZE)
    m.Params.Cuts=0 # turn off cuts
    m.Params.timeLimit=60 # 100-second time limit
    for i in G.nodes:
        Y[i,i].branchPriority=1
        
    #m.params.OutputFlag = 0

    m.params.LogFile='max_S_'+state+'.log'
    m.params.LogToConsole = 0    
    
    m.optimize()
    
    
    
    if m.status == GRB.OPTIMAL or m.status==GRB.TIME_LIMIT:
        return [i for i in G.nodes if Y[i,i].x > 0.5 ]
    else:
        return "Model status != optimal"    
   
def solve_heur_LFix(G, population, L, state, heur, k):
    m = gp.Model()
   
    # Y[i,j]=1 if vertices i and j belong to same component of G[S]
    # Y[i,i]=1 if i is selected in S
    S = m.addVars(G.nodes, vtype=GRB.BINARY)
    
    for j in range(k):
        m.addConstr(gp.quicksum(population[node["index"]] * S[node["index"]] for node in heur["nodes"] if node["district"] == j) <= L-1)
        
    
    
    for node1 in heur["nodes"]:
        i = node1["index"]
        for node2 in heur["nodes"]:
            j = node2["index"]
            if (i,j) in G.edges and node1["district"] != node2["district"]:
                m.addConstr(S[i] + S[j] <= 1)
                
    
    #m.addConstrs(S[node1] + S[node2] <= 1 for (node1,node2) in G.edges if node1["district"] != node2["district"])
    
    # the population reachable from j is less than L
    #m.addConstr( gp.quicksum(population[i] * Y[i] for i in G.nodes) <= L-1 )
    
    # Set objective function and solve MIP
    m.setObjective(gp.quicksum(S[i] for i in G.nodes), GRB.MAXIMIZE)
    m.Params.Cuts=0 # turn off cuts
    m.Params.timeLimit=100 # 100-second time limit

    m.params.LogFile='heur_LFix_'+state+'.log'
    m.params.LogToConsole = 0    
    
    m.optimize()
    
    
    
    if m.status == GRB.OPTIMAL or m.status==GRB.TIME_LIMIT:
        return [i for i in G.nodes if S[i].x > 0.5 ]
    else:
        return "Model status != optimal"     
    
def solve_max_indep(G, population, L, state):
    m = gp.Model()
   
    # Y[i,j]=1 if vertices i and j belong to same component of G[S]
    # Y[i,i]=1 if i is selected in S
    Y = m.addVars(G.nodes, vtype=GRB.BINARY)
    
    m.addConstrs(Y[i] + Y[j] <= 1 for (i,j) in G.edges)
    
    # the population reachable from j is less than L
    #m.addConstr( gp.quicksum(population[i] * Y[i] for i in G.nodes) <= L-1 )
    
    # Set objective function and solve MIP
    m.setObjective(gp.quicksum(Y[i] for i in G.nodes), GRB.MAXIMIZE)
    m.Params.Cuts=0 # turn off cuts
    m.Params.timeLimit=100 # 100-second time limit

    m.params.LogFile='max_indep_'+state+'.log'
    m.params.LogToConsole = 0    
    
    m.optimize()
    
    
    
    if m.status == GRB.OPTIMAL or m.status==GRB.TIME_LIMIT:
        return [i for i in G.nodes if Y[i].x > 0.5 ]
    else:
        return "Model status != optimal"    
