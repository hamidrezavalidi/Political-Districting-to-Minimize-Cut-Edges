import networkx as nx

# how many people are reachable from v in G[S]? Uses BFS
def reachable_population(G, population, S, v):
    pr = 0 # population reached
    if not S[v]:
        return 0
    
    visited = [False for i in G.nodes]
    child = [v]
    visited[v] = True
    while child:
        parent = child
        child = list()
        for i in parent:
            pr += population[i]
            for j in G.neighbors(i):
                if S[j] and not visited[j]:
                    child.append(j)
                    visited[j] = True
    return pr   


def do_Hess_DFixing(m, G, position):
    DFixings = 0
    for i in G.nodes:
        for j in G.nodes:
            if position[i] < position[j] and m._X[i,j].UB > 0.5:
                m._X[i,j].UB = 0
                DFixings += 1
    m.update()
    return DFixings


def do_Hess_LFixing(m, G, population, L, ordering):
    LFixings = 0
    S = [True for v in G.nodes]
    for j in ordering:       
        if reachable_population(G, population, S, j) < L:
            for i in G.nodes:
                if m._X[i,j].UB > 0.5:
                    m._X[i,j].UB = 0
                    LFixings += 1
        S[j] = False
        
    m.update()
    return LFixings


def do_Hess_LFixing_without_Contiguity(m, G, population, L, ordering):
    LFixings = 0
    remaining_population = sum(population[v] for v in G.nodes)
    for j in ordering:
        if remaining_population < L:
            for i in G.nodes:
                if m._X[i,j].UB > 0.5:
                    m._X[i,j].UB = 0
                    LFixings += 1
        remaining_population -= population[j]
    m.update()
    return LFixings


def do_Hess_UFixing(m, DG, population, U, ordering): 
    UFixings = 0
    
    for (i,j) in DG.edges:
        DG[i][j]['ufixweight'] = population[j] # weight of edge (i,j) is population of its head j
       
    for j in ordering:
        dist = nx.shortest_path_length(DG,source=j,weight='ufixweight')
        for i in DG.nodes:
            if i != j and dist[i]+population[j] > U and m._X[i,j].UB > 0.5: 
                m._X[i,j].UB = 0
                UFixings += 1
                
        # we should "remove" vertex j from the graph for subsequent distance calculations, so give incoming edges large weights
        for i in DG.neighbors(j):
            DG[i][j]['ufixweight'] = U+1

    m.update()
    return UFixings


def do_Hess_UFixing_without_Contiguity(m, G, population, U): 
    UFixings = 0
    for i in G.nodes:
        for j in G.nodes:
            if i != j and population[i] + population[j] > U and m._X[i,j].UB > 0.5:
                m._X[i,j].UB = 0
                UFixings += 1
    m.update()
    return UFixings
    

def do_Hess_ZFixing(m, G):
    ZFixings = 0
    for u,v in G.edges:
        for j in G.nodes:
            if m._X[u,j].UB < 0.5 or m._X[v,j].LB > 0.5:
                m._Z[u,v,j].UB = 0
                ZFixings += 1
    m.update()
    return ZFixings


def do_Labeling_DFixing(m, G, ordering, k):
    DFixings = 0
    for p in range(G.number_of_nodes()):
        i = ordering[p]
        for j in range(p+1,k):
            if m._X[i,j].UB > 0.5:
                m._X[i,j].UB = 0
                DFixings += 1
    m.update()
    return DFixings


def do_Labeling_ZFixing(m, G, k):
    ZFixings = 0
    for u,v in G.edges:
        for j in range(k):
            if m._X[u,j].UB < 0.5 or m._X[v,j].LB > 0.5:
                m._Z[u,v,j].UB = 0
                ZFixings += 1
            elif m._X[u,j].LB > 0.5 and m._X[v,j].UB < 0.5:
                m._Z[u,v,j].LB = 1
                ZFixings += 1
    m.update()
    return ZFixings


def do_Labeling_LFixing(m, G, population, L, ordering, k):
    LFixings = 0
    
    # find "back" of ordering B = {v_q, v_{q+1}, ..., v_{n-1} }
    n = G.number_of_nodes()
    S = [False for v in G.nodes]
    for p in range(n):
        v_pos = n - p - 1
        v = ordering[v_pos]
        S[v] = True
        pr = reachable_population(G, population, S, v)
        if pr >= L:
            q = v_pos + 1
            break
    
    # none of the vertices at back (in B) can root a district. 
    for p in range(q,n):
        i = ordering[p]
        for j in range(k):
            if m._R[i,j].UB > 0.5:
                m._R[i,j].UB = 0
                LFixings += 1
    
    # vertex v_{q-1} cannot root districts {0, 1, ..., k-2}
    # vertex v_{q-2} cannot root districts {0, 1, ..., k-3}
    # ... 
    # vertex v_{q-t} cannot root districts {0, 1, ..., k-t-1}
    # ...
    # vertex v_{q-(k-1)} cannot root district {0}
    for t in range(1,k):
        i = ordering[q-t]
        for j in range(k-t):
            if m._R[i,j].UB > 0.5:
                m._R[i,j].UB = 0
                LFixings += 1
    
    m.update()
    return LFixings 


def do_Labeling_LFixing_without_Contiguity(m, G, population, L, ordering, k):
    LFixings = 0
    
    # find "back" of ordering B = {v_q, v_{q+1}, ..., v_{n-1} }
    n = G.number_of_nodes()
    cumulative_population = 0
    for p in range(n):
        v = n - p - 1
        cumulative_population += population[v]
        if cumulative_population >= L:
            q = p + 1
            break
    
    # none of the vertices at back (in B) can root a district. 
    for p in range(q,n):
        i = ordering[p]
        for j in range(k):
            if m._R[i,j].UB > 0.5:
                m._R[i,j].UB = 0
                LFixings += 1
    
    # vertex v_{q-1} cannot root districts {0, 1, ..., k-2}
    # vertex v_{q-2} cannot root districts {0, 1, ..., k-3}
    # ... 
    # vertex v_{q-t} cannot root districts {0, 1, ..., k-t-1}
    # ...
    # vertex v_{q-(k-1)} cannot root district {0}
    for t in range(1,k):
        i = ordering[q-t]
        for j in range(k-t):
            if m._R[i,j].UB > 0.5:
                m._R[i,j].UB = 0
                LFixings += 1
    
    m.update()
    return LFixings 


def do_Labeling_UFixing(m, DG, population, U, ordering, k):
    UFixings_X = 0
    UFixings_R = 0
    DG = m._DG
    for (i,j) in DG.edges:
        DG[i][j]['ufixweight'] = population[j] # weight of edge (i,j) is population of its head j
    
    for j in range(k):
        
        v = ordering[j]
        dist = nx.shortest_path_length(DG,source=v,weight='ufixweight')
        
        if j == 0:
            min_dist = U+1
        else:
            min_dist = min( dist[ordering[t]] + population[v] for t in range(j) )
        
        if min_dist <= U:
            break
        
        if m._R[v,j].LB < 0.5:
            m._R[v,j].LB = 1
            UFixings_R += 1
            
        if m._X[v,j].LB < 0.5:
            m._X[v,j].LB = 1
            UFixings_X += 1
            
        for t in range(k):
            if t != j and m._X[v,t].UB > 0.5:
                m._X[v,t].UB = 0
                UFixings_X += 1
        
        for i in DG.nodes:
            if i != v and m._R[i,j].UB > 0.5:
                m._R[i,j].UB = 0
                UFixings_R += 1
        
        for i in DG.nodes:
            if i != v and dist[i] + population[v] > U and m._X[i,j].UB > 0.5:
                m._X[i,j].UB = 0
                UFixings_X += 1
        
    m.update()
    return UFixings_X, UFixings_R    


def do_labeling_UFixing_without_Contiguity():
    # no fixings possible?
    UFixings_X = 0
    UFixings_R = 0
    return UFixings_X, UFixings_R    

