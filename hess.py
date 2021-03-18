import gurobipy as gp
from gurobipy import GRB 

import bisect
import networkx as nx

def add_base_constraints(m, population, L, U, k):
    DG = m._DG
    # Each vertex i assigned to one district
    m.addConstrs(gp.quicksum(m._X[i,j] for j in DG.nodes) == 1 for i in DG.nodes)
     
    # Pick k centers
    m.addConstr(gp.quicksum(m._X[j,j] for j in DG.nodes) == k)
    
    # Population balance: population assigned to vertex j should be in [L,U], if j is a center
    m.addConstrs(gp.quicksum(population[i] * m._X[i,j] for i in DG.nodes) <= U * m._X[j,j] for j in DG.nodes)
    m.addConstrs(gp.quicksum(population[i] * m._X[i,j] for i in DG.nodes) >= L * m._X[j,j] for j in DG.nodes)
    
    # Add coupling inequalities for added model strength
    couplingConstrs = m.addConstrs(m._X[i,j] <= m._X[j,j] for i in DG.nodes for j in DG.nodes)
    
    # Make them user cuts
    for i in DG.nodes:
        for j in DG.nodes:
            couplingConstrs[i,j].Lazy = -1
    
    # Set branch priority on center vars
    for j in DG.nodes:
        m._X[j,j].BranchPriority=1         

        
def add_objective(m, G):
    # Y[i,j] = 1 if edge {i,j} is cut
    m._Y = m.addVars(G.edges, vtype=GRB.BINARY)
    m.addConstrs( m._X[i,v]-m._X[j,v] <= m._Y[i,j] for i,j in G.edges for v in G.nodes)
    m.setObjective( gp.quicksum(m._Y), GRB.MINIMIZE )
    

def add_extended_objective(m, G):
    # Z[i,j,v] = 1 if edge (i,j) is cut because i->v but j!->v
    m._Z = m.addVars(G.edges, G.nodes, vtype=GRB.BINARY) 
    m.addConstrs( m._X[i,v]-m._X[j,v] <= m._Z[i,j,v] for i,j in G.edges for v in G.nodes)
    m.setObjective( gp.quicksum(m._Z), GRB.MINIMIZE )
    
    
def compute_bigM_shir(DG, population, U):
    # Compute small values for big-M in shir model that depend on i and j.
    M = [[0 for i in DG.nodes] for j in DG.nodes]
    for (i,j) in DG.edges:
        DG[i][j]['weight']=population[j] # weight of edge (i,j) is population of its head j
    
    cumulative_population = [0 for i in range(len(population))]
    t=0
    for ipopulation in sorted(population):
        if t==0:
            cumulative_population[t] = ipopulation
        else:
            cumulative_population[t] = cumulative_population[t-1]+ipopulation
        t += 1
    
    for j in DG.nodes:
        dist = nx.shortest_path_length(DG,source=j,weight='weight')
        for i in DG.nodes:
            R = int(U - ( dist[i] + population[j] ))
            if R<0: # infeasible
                M[i][j] = 0
            else:   # how many nodes could be reached after getting to i from j? (then +1 for zero-based, and another +1 for i)
                M[i][j] = (bisect.bisect(cumulative_population, R)+1) 

    return M    
    
    
def add_shir_constraints(m, population, U):
    DG = m._DG
    # F[j,u,v] tells how much flow (from source j) is sent across arc (u,v)
    F = m.addVars( DG.nodes, DG.edges, vtype=GRB.CONTINUOUS)
    M = compute_bigM_shir(DG, population, U)   
    m.addConstrs( gp.quicksum(F[j,u,j] for u in DG.neighbors(j)) == 0 for j in DG.nodes)
    m.addConstrs( gp.quicksum( F[j,u,i]-F[j,i,u] for u in DG.neighbors(i) ) == m._X[i,j] for i in DG.nodes for j in DG.nodes if i!=j)
    m.addConstrs( gp.quicksum( F[j,u,i] for u in DG.neighbors(i) ) <= M[i][j]*m._X[i,j] for i in DG.nodes for j in DG.nodes if i!=j)
    m.update()
      
        
def add_scf_constraints(m, G, extended):
    DG = m._DG
    # FIXME: big M
    M = DG.number_of_nodes() - m._k + 1
    # F[u,v] tells how much flow is sent across arc (u,v)
    F = m.addVars( DG.edges, vtype=GRB.CONTINUOUS )
    m.addConstrs( gp.quicksum(m._X[i,j] for i in DG.nodes) == gp.quicksum(F[j,u]-F[u,j] for u in DG.neighbors(j)) + 1 for j in DG.nodes)
    m.addConstrs( gp.quicksum(F[u,j] for u in DG.neighbors(j)) <= M*(1-m._X[j,j]) for j in DG.nodes)
        
    if extended:
        m.addConstrs( F[i,j] + F[j,i] <= M*(1 - gp.quicksum(m._Z[i,j,v] for v in G.nodes)) for i,j in G.edges)
    else:
        m.addConstrs( F[i,j] + F[j,i] <= M*(1 - m._Y[i,j]) for i,j in G.edges)
        
