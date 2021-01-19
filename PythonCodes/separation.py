from gurobipy import GRB
import gurobipy as gp 
import networkx as nx

def labeling_cut_callback(m, where):
    labeling_callback(m, where, False)
       
def labeling_lcut_callback(m, where):
    labeling_callback(m, where, True)
    
def labeling_Austin_cut_callback(m, where):
    labeling_Austin_callback(m, where, False)    
    
def labeling_Austin_lcut_callback(m, where):
    labeling_Austin_callback(m, where, True)    

def hess_cut_callback(m, where):
    hess_callback(m, where, False)
        
def hess_lcut_callback(m, where):
    hess_callback(m, where, True)

def labeling_Austin_callback(m, where, lcut_minimalization):
    if where == GRB.Callback.MIPSOL:
        xval = m.cbGetSolution(m._X)
        rval = m.cbGetSolution(m._R)
        DG = m._DG
        U = m._U
        k = m._k
        population = m._population
        
        for j in range(k):
            V_j = [v for v in DG.nodes if xval[v,j] > 0.5]
                        
            # find vertex b
            for i in V_j:
                if rval[i,j] > 0.5:
                    b = i
                    
            # find a,b sep
            for component in nx.strongly_connected_components(DG.subgraph(V_j)):
                if b in component:
                    continue
                
                max_pop = max(population[v] for v in component)
                for v in component:
                    if population[v] == max_pop:
                        a = v
                        break
                
                
                # create G'
                boundary_vertices = [v for v in nx.node_boundary(DG, component, None)]
                hidden_edges = [(i,j) for (i,j) in DG.edges if i in boundary_vertices or j in component]                
                DG_prime = nx.restricted_view(DG, [], hidden_edges)
                
                if lcut_minimalization:
                    for (u,v) in DG.edges():
                            DG[u][v]['weight'] = population[v]                            
                    lengths_from_b = nx.single_source_dijkstra_path_length(DG_prime, b, cutoff=U-population[b], weight='weight')
                else:
                    lengths_from_b = nx.single_source_shortest_path_length(DG_prime, b)
                    
                C = [v for v in boundary_vertices if v in lengths_from_b]                           
                m.cbLazy( m._X[a,j] <= gp.quicksum(m._X[c,j] for c in C) + gp.quicksum(m._R[s,j] for s in component))
        

def hess_callback(m, where, lcut_minimalization):
    if where == GRB.Callback.MIPSOL:
        xval = m.cbGetSolution(m._X)
        DG = m._DG
        U = m._U
        population = m._population
        centers = [j for j in DG.nodes if xval[j,j]>0.5]
               
        for b in centers:
            # create component of b
            S = [v for v in DG.nodes if xval[v,b]>0.5]
            
            for component in nx.strongly_connected_components(DG.subgraph(S)):
                if b in component: 
                    continue
                
                # find vertex a
                max_pop = max(population[v] for v in component)
                for v in component:
                    if population[v] == max_pop:
                        a = v
                        break   
                    
                # create G'
                boundary_vertices = [v for v in nx.node_boundary(DG, component, None)]
                hidden_edges = [(i,j) for (i,j) in DG.edges if i in boundary_vertices or j in component]                
                DG_prime = nx.restricted_view(DG, [], hidden_edges)
                
                if lcut_minimalization:
                    for (u,v) in DG.edges():
                            DG[u][v]['weight'] = population[v]                            
                    lengths_from_b = nx.single_source_dijkstra_path_length(DG_prime, b, cutoff=U-population[b], weight='weight')
                else:
                    lengths_from_b = nx.single_source_shortest_path_length(DG_prime, b)
                    
                C = [v for v in boundary_vertices if v in lengths_from_b]                           
                m.cbLazy( m._X[a,b] <= gp.quicksum(m._X[c,b] for c in C) )    
                
                                                                        
def labeling_callback(m, where, lcut_minimalization):
    if where == GRB.Callback.MIPSOL:
        DG = m._DG
        xval = m.cbGetSolution(m._X)       
        population = m._population
        U = m._U
        k = m._k
            
        for j in range(k):
            V_j = [v for v in DG.nodes if xval[v,j] > 0.5]
            b = -1
            for component in sorted(nx.strongly_connected_components(DG.subgraph(V_j)), key = len, reverse = True):
                if b == -1:
                    max_pop = max(population[v] for v in component)
                    for v in component:
                        if population[v] == max_pop:
                            b = v
                            break
                    continue
                    
                # find vertex "a" with maximum population
                max_pop = max(population[v] for v in component)
                for v in component:
                    if population[v] == max_pop:
                        a = v
                        break
                
                # create DG'
                boundary_vertices = [v for v in nx.node_boundary(DG, component, None)]
                hidden_edges = [(i,j) for (i,j) in DG.edges if i in boundary_vertices or j in component]    
                DG_prime = nx.restricted_view(DG, [], hidden_edges)
                
                if lcut_minimalization:
                    for (u,v) in DG.edges():
                        DG[u][v]['weight'] = population[v]                        
                    lengths_from_b = nx.single_source_dijkstra_path_length(DG_prime, b, cutoff=U-population[b], weight='weight')                               
                else:
                    lengths_from_b = nx.single_source_shortest_path_length(DG_prime, b)
                    
                C = [v for v in boundary_vertices if v in lengths_from_b]                                            
                m.cbLazy( m._X[a,j] + m._X[b,j] <= 1 + gp.quicksum(m._X[c,j] for c in C) )    
                    

                        