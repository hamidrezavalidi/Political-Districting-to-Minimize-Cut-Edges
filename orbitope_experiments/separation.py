import gurobipy as gp
from gurobipy import GRB 

import networkx as nx

def find_fischetti_separator(DG, component, b):
    neighbors_component = [False for i in DG.nodes]
    for i in nx.node_boundary(DG, component, None):
        neighbors_component[i] = True
    
    visited = [False for i in DG.nodes]
    child = [b]
    visited[b] = True
    
    while child:
        parent = child
        child = []
        for i in parent:
            if not neighbors_component[i]:
                for j in DG.neighbors(i):
                    if not visited[j]:
                        child.append(j)
                        visited[j] = True
    
    C = [ i for i in DG.nodes if neighbors_component[i] and visited[i] ]
    return C


def lcut_separation_generic(m, where):
    if where == gp.GRB.Callback.MIPNODE and m.cbGet(gp.GRB.Callback.MIPNODE_STATUS) == gp.GRB.OPTIMAL and m._symmetry == 'sci':
            
        m._numOrbitopeCallbacks += 1
        
        #Get current fractional solution
        xval = m.cbGetNodeRel(m._X)
        
        OG = m._orbitope_graph
        
        k = m._k
        
        ordering = m._ordering
        
        n = len(ordering)
        
        epsilon = 0.01
        
        # Get x(bar_ij) for each leader (i,j).
        # See Corollary 11 of Kaibel and Pfetsch (MP, 2008).
        # Or bottom of page 688 of Faenza and Kaibel (MOR, 2009).
        #
        xbar = xval.copy()
        for j in range(k-2,-1,-1):
            for i in range(n):
                xbar[i,j] += xbar[i,j+1]
                
        position = m._position
        
        # Set edge weights. By Faenza and Kaibel (MOR, 2009; pp.688),
        #   we should set weights of vertical arc into (i,j) to be x_ij.
        # 
        for i,j in OG.edges:
            OG.edges[i,j]['edge_weight'] = 0
    
        for p in range(1,n):
    
            i = ordering[p]
            pi = ordering[p-1] # 'previous' i
    
            for j in range(k):
                node = (i,j)
                pnode = (pi,j)
                OG.edges[pnode,node]['edge_weight'] = max(0, xval[i,j]) # vertical arc
                
        # Now solve shortest path problem for each of the first k vertices along 
        #   the diagonal of the orbitope graph OG. The paths \Gamma will give
        #   us the type S vertices along the path S=S(\Gamma). See page 688.
        
        for p in range(k):
        
            v = ordering[p]
    
            path_start = (v,p)
            (distance,path) = nx.single_source_dijkstra(OG,weight='edge_weight',source=path_start)
    
            for node in distance.keys():
    
                if node == 't':
                    continue
    
                (i,j) = node
                my_path = path[node]
                S = list()
    
                for pos in range(len(my_path)):
    
                    path_node = my_path[pos]
                    (path_node_i, path_node_j) = path_node
    
                    if path_node == path_start:
                        S.append(path_node)
                    else:
                        prev_path_node = my_path[pos-1]
                        (prev_i,prev_j) = prev_path_node
                        if path_node_j == prev_j:
                            S.append(path_node)
    
                pos_i = position[i]
    
                if pos_i < n-1 and j < k-1:
    
                    next_i = ordering[pos_i+1]
                    S_weight = sum( xval[myi,myj] for (myi,myj) in S )
                    
                    
    
                    if xbar[next_i,j+1] - S_weight > epsilon:
                        #print("S: ", S)
                        B = [ (next_i,jval) for jval in range(j+1,k) ]
                        m.cbCut(gp.quicksum(m._X[i_val, j_val] for (i_val, j_val) in B) <= gp.quicksum(m._X[u_val, v_val] for (u_val, v_val) in S))
                        m._numOrbitopeCuts += 1
    
    if where == GRB.Callback.MIPSOL:
        m._numCallbacks += 1
        xval = m.cbGetSolution(m._X)
        DG = m._DG
        U = m._U
        base = m._base
        population = m._population
        k = m._k
        
        if base == 'labeling':
            district_labels = [j for j in range(k)]
        elif base == 'hess':
            district_labels = [j for j in DG.nodes if xval[j,j] > 0.5]
               
        for j in district_labels:
            
            # vertices assigned to this district (label j)
            S = [v for v in DG.nodes if xval[v,j] > 0.5]
            
            # what shall we deem as the "root" of this district? call it b
            if base == 'hess':
                b = j
            elif base == 'labeling':
                max_cp = 0
                max_component = []
                
                for component in nx.strongly_connected_components(DG.subgraph(S)):
                    if max_cp < sum(population[v] for v in component):
                        max_cp = sum(population[v] for v in component)
                        max_component = component
                    
                # find some vertex "b" that has largest population in this component
                maxpb = max(population[v] for v in max_component)
                maxpb_vertices = [ v for v in max_component if population[v] == maxpb ]
                b = maxpb_vertices[0]  
            
            for component in nx.strongly_connected_components(DG.subgraph(S)):
                
                if b in component: 
                    continue
                
                # find some vertex "a" that has largest population in this component
                maxp = max(population[v] for v in component)
                maxp_vertices = [ v for v in component if population[v] == maxp ]
                a = maxp_vertices[0]  
                    
                # get minimal a,b-separator
                C = find_fischetti_separator(DG, component, b)
                
                # make it a minimal *length-U* a,b-separator
                for (u,v) in DG.edges():
                    DG[u][v]['lcutweight'] = population[u]   
                    
                # "remove" C from graph
                for c in C:
                    for node in DG.neighbors(c):
                        DG[c][node]['lcutweight'] = U+1
                
                # is C\{c} a length-U a,b-separator still? If so, remove c from C
                drop_from_C = []
                for c in C:
                    
                    # temporarily add c back to graph (i.e., "remove" c from cut C)
                    for node in DG.neighbors(c):
                        DG[c][node]['lcutweight'] = population[c]
                    
                    # what is distance from a to b in G-C ?
                    distance_from_a = nx.single_source_dijkstra_path_length(DG, a, weight='lcutweight')
                    
                    if distance_from_a[b] + population[b] > U:
                        # c was not needed in the cut C. delete c from C
                        drop_from_C.append(c)
                    else:
                        # keep c in C. revert arc weights back to "infinity"
                        for node in DG.neighbors(c):
                            DG[c][node]['lcutweight'] = U+1
                    
                # add lazy cut
                minC = [c for c in C if c not in drop_from_C ]
                if base == 'hess':
                    m.cbLazy( m._X[a,b] <= gp.quicksum(m._X[c,b] for c in minC) )    
                elif base == 'labeling':
                    m.cbLazy( m._X[a,j] + m._X[b,j] <= 1 + gp.quicksum(m._X[c,j] for c in minC) )
                m._numLazyCuts += 1

        