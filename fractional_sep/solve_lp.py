import gurobipy as gp
from gurobipy import GRB 

def labeling_relaxation(m, G, k):
    violation = 1.00
    m.relax()
    
    while violation > 0.001:
        m.optimize()
        xval = [ [m._X[i,j].x  for i in G.nodes] for j in range(k)]
        
        yval = [ m._Y[u,v].x for (u,v) in G.edges ]
        
        for (u,v) in G.edges:
            lhs_set = []
            lhs_val = 0
            for j in G.nodes:
                if xval[u][j] - xval[v][j] > 0:
                    lhs_val += xval[u,j] - xval[v,j]
                    lhs_set.append(j)
                
            
        
        