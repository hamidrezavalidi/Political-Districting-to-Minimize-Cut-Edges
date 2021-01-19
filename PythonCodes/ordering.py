import gurobipy as gp
from gurobipy import GRB 

def construct_position(ordering):
    position = [-1 for i in range(len(ordering))]
    for p in range(len(ordering)):
        v = ordering[p]
        position[v] = p
    return position

def find_ordering(G, population, L, heur, k, OrderChoice, state):
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
        S = solve_maxS_recursive_model(G,population,L,state)
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
        
def solve_maxS_recursive_model(G,population,L,state):
    m = gp.Model()
   
    # Y[i,j]=1 if vertices i and j belong to same component of G[S]
    # Y[i,i]=1 if i is selected in S
    Y = m.addVars(G.nodes, G.nodes, vtype=GRB.BINARY)
    
    # warm start
    #for i in warm:
     #   Y[i,i].start = 1.0
    
    m.addConstrs(Y[i,j] == Y[j,i] for i in G.nodes for j in G.nodes)
    m.addConstrs(Y[i,j] <= Y[i,i] for i in G.nodes for j in G.nodes)
    
    # the population reachable from j is less than L
    m.addConstrs( gp.quicksum(population[i] * Y[i,j] for i in G.nodes) <= (L-1)*Y[j,j] for j in G.nodes)
    
    # if G[S] contains an i,v-path and (v,j) is an edge and j is selected, 
    #    then G[S] contains an i,j path
    for v in G.nodes:
        m.addConstrs( Y[i,v] + Y[j,j] <= 1 + Y[i,j] for i in G.nodes for j in G.neighbors(v))
    
    # Set objective function and solve MIP
    m.setObjective(gp.quicksum(Y[i,i] for i in G.nodes), GRB.MAXIMIZE)
    m.Params.Cuts=0 # turn off cuts
    m.Params.timeLimit=100 # 100-second time limit
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
