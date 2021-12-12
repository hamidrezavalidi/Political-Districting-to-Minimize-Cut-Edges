import gurobipy as gp
from gurobipy import GRB 


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
    # Y[i,j] = 1 if edge {i,j} is cut
    m._Y = m.addVars(G.edges, vtype=GRB.BINARY)
    m.addConstrs( m._X[i,v]-m._X[j,v] <= m._Y[i,j] for i,j in G.edges for v in G.nodes)
    m.setObjective( gp.quicksum(m._Y), GRB.MINIMIZE )

    
    
def most_possible_nodes_in_one_district(population, U):
    cumulative_population = 0
    num_nodes = 0
    for ipopulation in sorted(population):
        cumulative_population += ipopulation
        num_nodes += 1
        if cumulative_population > U:
            return num_nodes - 1
   
    
def add_shir_constraints(m):
    DG = m._DG
    
    # F[j,u,v] tells how much flow (from source j) is sent across arc (u,v)
    F = m.addVars( DG.nodes, DG.edges, vtype=GRB.CONTINUOUS)
    
    # compute big-M    
    M = most_possible_nodes_in_one_district(m._population, m._U) - 1
    
    m.addConstrs( gp.quicksum(F[j,u,j] for u in DG.neighbors(j)) == 0 for j in DG.nodes)
    m.addConstrs( gp.quicksum( F[j,u,i]-F[j,i,u] for u in DG.neighbors(i) ) == m._X[i,j] for i in DG.nodes for j in DG.nodes if i!=j)
    m.addConstrs( gp.quicksum( F[j,u,i] for u in DG.neighbors(i) ) <= M * m._X[i,j] for i in DG.nodes for j in DG.nodes if i!=j)
    m.update()
      
        
def add_scf_constraints(m, G, extended):
    DG = m._DG
    
    # F[u,v] tells how much flow is sent across arc (u,v)
    F = m.addVars( DG.edges, vtype=GRB.CONTINUOUS )
    
    # compute big-M
    M = most_possible_nodes_in_one_district(m._population, m._U) - 1
    
    m.addConstrs( gp.quicksum(m._X[i,j] for i in DG.nodes) == gp.quicksum(F[j,u]-F[u,j] for u in DG.neighbors(j)) + 1 for j in DG.nodes)
    m.addConstrs( gp.quicksum(F[u,j] for u in DG.neighbors(j)) <= M * (1-m._X[j,j]) for j in DG.nodes)
        
    if extended:
        m.addConstrs( F[i,j] + F[j,i] <= M * (1 - gp.quicksum(m._Z[i,j,v] for v in G.nodes)) for i,j in G.edges)
    else:
        m.addConstrs( F[i,j] + F[j,i] <= M * (1 - m._Y[i,j]) for i,j in G.edges)
        
