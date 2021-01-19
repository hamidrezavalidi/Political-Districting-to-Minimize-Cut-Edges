import gurobipy as gp
from gurobipy import GRB 

def add_labeling_cut_edges_objective(m, G, k):
    # Y[i,j] = 1 if edge {i,j} is cut
    #print("Hi Hamidddddddddddd!!!!!!!")
    m._Y = m.addVars(G.edges, vtype=GRB.BINARY)
    m.addConstrs( m._X[i,v]-m._X[j,v] <= m._Y[i,j] for i,j in G.edges for v in range(k))
    #m.addConstrs( m._X[j,v]-m._X[i,v] <= m._Y[i,j] for i,j in G.edges for v in range(k))
    m.setObjective( gp.quicksum(G[i][j]['edge_length']*m._Y[i,j] for i,j in G.edges), GRB.MINIMIZE )
    
def add_labeling_cut_edges_extended_objective(m, G, k):
    # Z[i,j,v] = 1 if edge (i,j) is cut because i->v but j!->v
    #print("Hi Hamidddddddddddd!!!!!!!")
    m._Z = m.addVars(G.edges, range(k), vtype=GRB.BINARY)
    m.addConstrs( m._X[i,v]-m._X[j,v] <= m._Z[i,j,v] for i,j in G.edges for v in range(k))
    m.setObjective(gp.quicksum(G[i][j]['edge_length']*m._Z[i,j,v] for i,j in G.edges for v in range(k)), GRB.MINIMIZE)

def build_labeling_model(m, population, L, U, k):
    DG = m._DG # bidirected version of G
    
    # Each vertex i assigned to one district
    m.addConstrs(gp.quicksum(m._X[i,j] for j in range(k)) == 1 for i in DG.nodes)
     
    # Population balance: population assigned to vertex j should be in [L,U]
    m.addConstrs(gp.quicksum(population[i] * m._X[i,j] for i in DG.nodes) <= U for j in range(k))
    m.addConstrs(gp.quicksum(population[i] * m._X[i,j] for i in DG.nodes) >= L for j in range(k)) 
    
def build_Austin_model(m,DG,k,ordered_vertices,orbitope):
    # add local cuts
    m.addConstrs(m._X[a,j] <= gp.quicksum(m._X[c,j] for c in DG.neighbors(a)) + m._R[a,j] for j in range(k) for a in DG.nodes)
        
    # coupling constraints
    m.addConstrs(m._R[v,j] <= m._X[v,j] for j in range(k) for v in DG.nodes)

    
        
def build_shir_model(m,DG,k,ordered_vertices,orbitope):
    g = m.addVars(DG.nodes, range(k), vtype=GRB.CONTINUOUS)
    f = m.addVars(range(k), DG.edges, vtype=GRB.CONTINUOUS)

    m.addConstrs(g[i,j] <= (DG.number_of_nodes()-k+1)*m._R[i,j] for i in DG.nodes for j in range(k))
    m.addConstrs(g[i,j] + gp.quicksum(f[j,u,i] for u in DG.neighbors(i)) == m._X[i,j] + gp.quicksum(f[j,i,u] for u in DG.neighbors(i)) for i in DG.nodes for j in range(k))
    m.addConstrs(gp.quicksum(f[j,u,i] for u in DG.neighbors(i)) <= (DG.number_of_nodes()-k)*(m._X[i,j]-m._R[i,j]) for i in DG.nodes for j in range(k))
    
    # the following constraints are weaker than some in the orbitope EF
    if orbitope==False:
        m.addConstrs(gp.quicksum(m._R[i,j] for i in DG.nodes)==1 for j in range(k))
        m.addConstrs(m._R[i,j] <= m._X[i,j] for i in DG.nodes for j in range(k))
        
def build_scf_model(m,G,DG,k,ordered_vertices,orbitope):
    f = m.addVars(DG.edges, vtype=GRB.CONTINUOUS)
    M = len(G.nodes) - k + 1
    
    m.addConstrs(gp.quicksum(f[u,v] for u in DG.neighbors(v)) - gp.quicksum(f[v,u] for u in DG.neighbors(v)) >= 1 - M*gp.quicksum(m._R[v,j] for j in range(k)) for v in G.nodes)
    
    for (i,j) in G.edges:
        if m._extended==True:
            m.addConstr( f[i,j] + f[j,i] <= M*(1 - gp.quicksum(m._Z[i,j,v] for v in range(k))) )
        else:
            m.addConstr( f[i,j] + f[j,i] <= M*(1 - m._Y[i,j]) )
            
    # the following constraints are weaker than some in the orbitope EF
    if orbitope==False:
        m.addConstrs(gp.quicksum(m._R[i,j] for i in DG.nodes)==1 for j in range(k))
        m.addConstrs(m._R[i,j] <= m._X[i,j] for i in DG.nodes for j in range(k))
    #m.addConstrs(g[i,j] <= (DG.number_of_nodes()-k+1)*m._R[i,j] for i in DG.nodes for j in range(k))
    #m.addConstrs(g[i,j] + gp.quicksum(f[j,u,i] for u in DG.neighbors(i)) == m._X[i,j] + gp.quicksum(f[j,i,u] for u in DG.neighbors(i)) for i in DG.nodes for j in range(k))
    #m.addConstrs(gp.quicksum(f[j,u,i] for u in DG.neighbors(i)) <= (DG.number_of_nodes()-k)*(m._X[i,j]-m._R[i,j]) for i in DG.nodes for j in range(k))
    
    # the following constraints are weaker than some in the orbitope EF
    if orbitope==False:
        m.addConstrs(gp.quicksum(m._R[i,j] for i in DG.nodes)==1 for j in range(k))
        m.addConstrs(m._R[i,j] <= m._X[i,j] for i in DG.nodes for j in range(k))        
       
def build_mcf_model(m,DG,k,population,orbitope):
    g = m.addVars(range(k), DG.nodes, DG.nodes, vtype=GRB.CONTINUOUS)
    f = m.addVars(range(k), DG.nodes, DG.edges, vtype=GRB.CONTINUOUS)
    
    for i in DG.nodes:
        for t in DG.nodes:
            for j in range(k):
                if i==t:
                    m.addConstr( m._R[i,j] + gp.quicksum(f[j,i,u,i] for u in DG.neighbors(i)) == m._X[i,j] )
                    m.addConstr( gp.quicksum(f[j,i,i,u] for u in DG.neighbors(i)) == 0 )
                else:
                    m.addConstr( g[j,t,i] + gp.quicksum(f[j,t,u,i] for u in DG.neighbors(i)) == gp.quicksum(f[j,t,i,u] for u in DG.neighbors(i)) )
                    m.addConstr( g[j,t,i] <= m._R[i,j] )
                    m.addConstr( gp.quicksum(f[j,t,u,i] for u in DG.neighbors(i)) <= m._X[i,j] - m._R[i,j] )
                    
    
    # the following constraints are weaker than some in the orbitope EF
    if orbitope==False:
        m.addConstrs(gp.quicksum(m._R[i,j] for i in DG.nodes)==1 for j in range(k))
        m.addConstrs(m._R[i,j] <= m._X[i,j] for i in DG.nodes for j in range(k))  
               
def add_orbitope_extended_formulation(m, G, k, ordering, position):
    s = m.addVars(G.nodes, range(k), vtype=GRB.CONTINUOUS) 
    u = m.addVars(G.nodes, range(k), vtype=GRB.CONTINUOUS) 
    w = m.addVars(G.nodes, range(k), vtype=GRB.CONTINUOUS) 
    
    m.addConstrs(m._X[i,j] == s[i,j]-s[i,j+1] for i in G.nodes for j in range(k-1))
    m.addConstrs(m._X[i,k-1] == s[i,k-1] for i in G.nodes)
    
    m.addConstrs(m._R[ordering[0],j] == w[ordering[0],j] for j in range(k))
    m.addConstrs(m._R[ordering[i],j] == w[ordering[i],j] - w[ordering[i-1],j] for i in range(1,G.number_of_nodes()) for j in range(k))
    
    m.addConstrs(m._R[i,j] <= m._X[i,j] for i in G.nodes for j in range(k))
    m.addConstrs(s[i,j] <= w[i,j] for i in G.nodes for j in range(k))
    
    m.addConstrs(u[ordering[i],j]+m._R[ordering[i],j] == u[ordering[i+1],j] + m._R[ordering[i+1],j+1] for i in range(0,G.number_of_nodes()-1) for j in range(k-1))
    m.addConstrs(u[ordering[i],k-1]+m._R[ordering[i],k-1] == u[ordering[i+1],k-1] for i in range(0,G.number_of_nodes()-1))
    m.addConstrs(u[ordering[G.number_of_nodes()-1],j]+m._R[ordering[G.number_of_nodes()-1],j] == 0 for j in range(k-1))
    
    m._R[ordering[0],0].LB=1
    m.addConstr( u[ordering[G.number_of_nodes()-1],k-1] + m._R[ordering[G.number_of_nodes()-1],k-1]==1 )  
