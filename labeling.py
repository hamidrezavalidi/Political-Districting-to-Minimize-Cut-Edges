import gurobipy as gp
from gurobipy import GRB 

def add_base_constraints(m, population, L, U, k):
    DG = m._DG # bidirected version of G
    
    # Each vertex i assigned to one district
    m.addConstrs(gp.quicksum(m._X[i,j] for j in range(k)) == 1 for i in DG.nodes)
     
    # Population balance: population assigned to vertex j should be in [L,U]
    m.addConstrs(gp.quicksum(population[i] * m._X[i,j] for i in DG.nodes) <= U for j in range(k))
    m.addConstrs(gp.quicksum(population[i] * m._X[i,j] for i in DG.nodes) >= L for j in range(k)) 
    
 
def add_objective(m, G, k):
    # Y[i,j] = 1 if edge {i,j} is cut
    m._Y = m.addVars(G.edges, vtype=GRB.BINARY)
    m.addConstrs( m._X[i,v]-m._X[j,v] <= m._Y[i,j] for i,j in G.edges for v in range(k))
    m.setObjective( gp.quicksum(m._Y[i,j] for i,j in G.edges), GRB.MINIMIZE )

    
def add_extended_objective(m, G, k):
    # Z[i,j,v] = 1 if edge (i,j) is cut because i->v but j!->v
    m._Z = m.addVars(G.edges, range(k), vtype=GRB.BINARY)
    m.addConstrs( m._X[i,v]-m._X[j,v] <= m._Z[i,j,v] for i,j in G.edges for v in range(k))
    m.setObjective(gp.quicksum(m._Z[i,j,v] for i,j in G.edges for v in range(k)), GRB.MINIMIZE)


def add_orbitope_extended_formulation(m, G, k, ordering):
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
   
        
def add_shir_constraints(m, k, symmetry):
    
    DG = m._DG
    # FIXME: need to make sure this matches up with model in overleaf document
    
    # g[i,j] = amount of flow generated at node i of type j
    g = m.addVars(DG.nodes, range(k), vtype=GRB.CONTINUOUS)
    
    # f[j,u,v] = amount of flow sent across arc uv of type j
    f = m.addVars(range(k), DG.edges, vtype=GRB.CONTINUOUS)

    # big-M
    M = DG.number_of_nodes()-k
    
    # add constraints
    m.addConstrs( g[i,j] <= (M+1)*m._R[i,j] for i in DG.nodes for j in range(k) )
    m.addConstrs( g[i,j] - m._X[i,j] == gp.quicksum(f[j,i,u]-f[j,u,i] for u in DG.neighbors(i)) for i in DG.nodes for j in range(k) )
    m.addConstrs( gp.quicksum(f[j,u,i] for u in DG.neighbors(i)) <= M*(m._X[i,j]-m._R[i,j]) for i in DG.nodes for j in range(k) )
    
    # the following constraints are weaker than some in the orbitope EF
    if symmetry != 'orbitope':
        m.addConstrs( gp.quicksum(m._R[i,j] for i in DG.nodes)==1 for j in range(k) )
        m.addConstrs( m._R[i,j] <= m._X[i,j] for i in DG.nodes for j in range(k) )
        

def add_scf_constraints(m, G, k, extended, symmetry):
    
    DG = m._DG
    # FIXME : make sure this matches overleaf document
    
    # f[u,v] = amount of flow sent across arc uv
    f = m.addVars(DG.edges, vtype=GRB.CONTINUOUS)
    
    # big M
    M = DG.number_of_nodes() - k + 1
    
    # if not a root, consume some flow.
    # if a root, only send out so much flow.
    m.addConstrs( gp.quicksum(f[u,v]-f[v,u] for u in DG.neighbors(v)) >= 1 - M*gp.quicksum(m._R[v,j] for j in range(k)) for v in G.nodes)
    
    # do not sent flow across cut edges
    if extended:
        m.addConstrs( f[i,j] + f[j,i] <= M*(1 - gp.quicksum( m._Z[i,j,v] for v in range(k) )) for (i,j) in G.edges)
    else:
        m.addConstrs( f[i,j] + f[j,i] <= M*(1 - m._Y[i,j]) for (i,j) in G.edges )
            
    # the following constraints are weaker than some in the orbitope EF
    if symmetry != 'orbitope':
        m.addConstrs( gp.quicksum(m._R[i,j] for i in DG.nodes)==1 for j in range(k) )
        m.addConstrs( m._R[i,j] <= m._X[i,j] for i in DG.nodes for j in range(k) )    
  
       
  # AUSTIN IS HERE???
