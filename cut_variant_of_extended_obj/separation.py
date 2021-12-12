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


def obj_hess_separation(m, where):
    
    if where != gp.GRB.Callback.MIPNODE:
        return
    
    if m.cbGet(gp.GRB.Callback.MIPNODE_STATUS) != gp.GRB.OPTIMAL:
        return
    
    m._numCallbacks += 1
    
    #Get current fractional solution
    xval = m.cbGetNodeRel(m._X)
    yval = m.cbGetNodeRel(m._Y)
    
    G = m._G
    
    epsilon = 0.01
    
    for (u,v) in G.edges:
        lhs_set = []
        lhs_val = 0
        
        for j in G.nodes:
            if xval[u,j] - xval[v,j] > 0:
                lhs_val += xval[u,j] - xval[v,j]
                lhs_set.append(j)
                
        if lhs_val - yval[u,v] > epsilon:  
            m.cbCut( gp.quicksum(m._X[u,j] - m._X[v,j] for j in lhs_set), GRB.LESS_EQUAL, m._Y[u,v] )
            m._numLazyCuts += 1
            
    nodeCount = m.cbGet(gp.GRB.Callback.MIPNODE_NODCNT)
    
    if nodeCount == 0:
        objValue = 0.0
        for (u,v) in G.edges:
            objValue += yval[u,v]
        print("obj value at node 0: ", objValue)
            
            
def obj_labeling_separation(m, where):
    
    if where != gp.GRB.Callback.MIPNODE:
        return
    
    if m.cbGet(gp.GRB.Callback.MIPNODE_STATUS) != gp.GRB.OPTIMAL:
        return
    
    m._numCallbacks += 1
    
    #Get current fractional solution
    xval = m.cbGetNodeRel(m._X)
    yval = m.cbGetNodeRel(m._Y)
    
    G = m._G
    k = m._k
    
    epsilon = 0.01
    
    threshold = pow(2, k)

    
    for (u,v) in G.edges:
        lhs_set = []
        lhs_val = 0
        
        for j in range(k):
            if xval[u,j] - xval[v,j] > 0:
                lhs_val += xval[u,j] - xval[v,j]
                lhs_set.append(j)
                  
                
        if lhs_val - yval[u,v] > epsilon and m._cutDict[(u, v)] < threshold:  
            m.cbCut( gp.quicksum(m._X[u,j] - m._X[v,j] for j in lhs_set) <= m._Y[u,v] )
            m._cutDict[(u, v)] += 1
            m._numLazyCuts += 1
            
    nodeCount = m.cbGet(gp.GRB.Callback.MIPNODE_NODCNT)
    
    if nodeCount == 0:
        objValue = 0.0
        for (u,v) in G.edges:
            objValue += yval[u,v]
        print("obj value at node 0: ", objValue)


def lcut_separation_generic(m, where):
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

        