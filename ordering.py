import gurobipy as gp
from gurobipy import GRB 

def sort_by_second(val):
    return val[1]


def construct_position(ordering):
    position = [-1 for i in range(len(ordering))]
    for p in range(len(ordering)):
        v = ordering[p]
        position[v] = p
    return position


def find_ordering(order, B, DG, population):
    if order == 'decreasing':
        nodes_with_population = [(i,population[i]) for i in DG.nodes]
        nodes_with_population.sort(key=sort_by_second,reverse=True)
        return [v for (v,p) in nodes_with_population]
    elif order == 'B_decreasing':
        V_B_with_population = [(i,population[i]) for i in DG.nodes if i not in B]
        V_B_with_population.sort(key=sort_by_second,reverse=True)
        return [v for (v,p) in V_B_with_population] + B
    else:
        return [v for v in DG.nodes]
    

def solve_maxB_problem(DG, population, L, k, heuristic_districts):
    m = gp.Model()
    m.params.LogToConsole = 0 # keep log to a minimum
    q = k
    
    # X[i,j]=1 if vertex i is assigned to bin j
    X = m.addVars(DG.nodes, range(q), vtype=GRB.BINARY)
    
    # B[i]=1 if vertex i is selected in set B
    B = m.addVars(DG.nodes, vtype=GRB.BINARY)
   
    # assignment constraints            
    m.addConstrs( gp.quicksum(X[i,j] for j in range(q)) == B[i] for i in DG.nodes )
                
    # bin population should be less than L
    m.addConstrs( gp.quicksum(population[i] * X[i,j] for i in DG.nodes) <= L-1 for j in range(q) )
    
    # bins shouldn't touch each other
    m.addConstrs( X[u,j] + B[v] <= 1 + X[v,j] for u,v in DG.edges for j in range(q) )
    
    # objective is to maximize size of set B
    m.setObjective( gp.quicksum( B[i] for i in DG.nodes), GRB.MAXIMIZE )
    
    m.Params.MIPFocus = 1 # turn on MIPFocus
    m.Params.timeLimit = 60 # 60-second time limit
    
    # suggest a (partial) warm start
    for district in heuristic_districts:
        for i in district:
            for t in range(q):
                if t != district:
                    X[i,t].start = 0.0
    
    m.optimize()
    sol = []
    if m.status in { GRB.OPTIMAL, GRB.TIME_LIMIT }:
        sol = [i for i in DG.nodes if B[i].x > 0.5 ]
        print("max B obj val =",m.objVal)
        
    return sol


