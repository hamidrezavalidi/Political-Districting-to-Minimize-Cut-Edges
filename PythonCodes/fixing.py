import networkx as nx

# how many people are reachable from v in G[S]? Using BFS
def reachable_population(G, population, S, v):
    pr = 0 # population reached
    if S[v]==False:
        return 0
    
    visited = [False for i in G.nodes]
    child = [v]
    while child:
        parent = child
        child=[]
        for i in parent:
            pr += population[i]
            for j in G.neighbors(i):
                if S[j]==True and visited[j]==False:
                    child.append(j)
                    visited[j]=True
    return pr   
'''        
def do_LFixing(m,G,population,L,ordered_vertices,position):
    LFixed = 0
    S = [True for v in G.nodes]
    for v in ordered_vertices:
        if reachable_population(G, population, S, v) < L:
            for i in G.nodes:
                if m._X[i,v].UB>0.5: #position[i]>position[v]:
                    print("LFixed: ",(i,v))
                    m._X[i,v].UB=0
                    LFixed += 1
        S[v]=False
    print("Number of LFixings =",LFixed,"out of",G.number_of_nodes()*G.number_of_nodes())
    return LFixed
'''
def do_LFixing(m,G,population,L,ordered_vertices,position):
    LFixed = 0
    #S = [True for v in G.nodes]
    for v in ordered_vertices:
        if v in m._S:
            for i in G.nodes:
                #print("LFixed: ",(i,v))
                if position[i] >= position[v]:                   
                    m._X[i,v].UB=0
                    LFixed += 1
        #S[v]=False
    print("Number of LFixings =",LFixed,"out of",G.number_of_nodes()*G.number_of_nodes())
    return LFixed

def do_UFixing(m,population,U,ordered_vertices,position):
    DG = m._DG
    UFixed = 0
    for (i,j) in DG.edges:
        DG[i][j]['weight'] = population[j] # weight of edge (i,j) is population of its head j
    
    #for i in ordered_vertices:
     #   print ("vertex ", i, "with position ", position[i])
    
    for j in ordered_vertices:
        dist = nx.shortest_path_length(DG,source=j,weight='weight')
        for i in DG.nodes:
            if dist[i]+population[j]>U and i!=j and m._X[i,j].UB>0.5: #position[i]>position[j]:
                #print("UFixed: ", (i+1,j+1))
                m._X[i,j].UB=0
                UFixed += 1
                
        # we should "remove" vertex j from the graph for subsequent distance calculations, so give incoming edges large weights
        for i in DG.neighbors(j):
            DG[i][j]['weight']=U+1
    print("Number of UFixings =",UFixed,"out of",DG.number_of_nodes()*DG.number_of_nodes())
    return UFixed
    
def do_ZFixing(X,Z,G):
    ZFixed = 0
    for i,j in G.edges:
        for v in G.nodes:
            if X[i,v].UB<0.5:
                Z[i,j,v].UB=0
                ZFixed += 1
    print("Number of ZFixings =",ZFixed,"out of",G.number_of_nodes()*G.number_of_edges())
    return ZFixed

def do_labeling_ZFixing(m,G,k):
    ZFixed = 0
    for i,j in G.edges:
        for v in range(k):
            if m._X[i,v].UB<0.5:
                m._Z[i,j,v].UB=0
                ZFixed += 1
    print("Number of ZFixings =",ZFixed,"out of",G.number_of_nodes()*k)
    return ZFixed
     
def do_DiagFixing(m,G,position):
    DiagFixed = 0
    for i in G.nodes:
        for j in G.nodes:
            if position[i]<position[j]:
                #print("DiagFixing: ", (i,j))
                m._X[i,j].UB=0
                DiagFixed += 1
    print("Number of DiagFixings =",DiagFixed,"out of",G.number_of_nodes()*G.number_of_nodes())
    return DiagFixed

def do_labeling_UFixing(m,G,population,U,ordering,k):
    UFixed = 0
    DG = m._DG
    for (i,j) in DG.edges:
        DG[i][j]['weight'] = population[j] # weight of edge (i,j) is population of its head j
    
    for district in range(k):
        # fix to one
        '''
        if district == 0:
            m._X[ordering[district],district].LB=1
            UFixed += 1
        else:
            counter = 0
            for j in range(district):
                total_population = population[ordering[district]] + population[ordering[j]]
                if total_population > U:
                    counter += 1                   
                else:
                    break
                if counter == district: 
                    m._X[ordering[district],district].LB=1
                    UFixed += 1
        '''            
        
        # fix to zero
        fix = False
        if district != 0:
            counter = 0
            for j in range(district):
                total_population = population[ordering[district]] + population[ordering[j]]
                if total_population > U:
                    counter += 1
                else:
                    break
            if counter == district: fix = True
                
             
        if district == 0 or fix:
            center = ordering[district]
            dist = nx.shortest_path_length(DG,source=center,weight='weight')
            #print("state: ", m._state)
            for i in DG.nodes:
                #print("node: ", i,"district:", district)
                if dist[i]+population[center]>U and i!=center and m._X[i,district].UB > 0.5:
                    m._X[i,district].UB=0
                    UFixed += 1
        # we should "remove" vertex "center" from the graph for subsequent distance calculations, so give incoming edges large weights
        for i in DG.neighbors(center):
            DG[i][center]['weight']=U+1
        
        else:
            break
            
    print("Number of UFixings =",UFixed,"out of",DG.number_of_nodes()*k)
    return UFixed 
    '''   
    counter = 0
    j = ordering[counter]
    dist = nx.shortest_path_length(DG,source=j,weight='weight')
    for i in DG.nodes:
            if dist[i]+population[j]>U and i!=j and m._X[i,j].UB>0.5:
                m._X[i,counter].UB=0
                UFixed += 1
    '''
def do_labeling_LFixing(m,G,population,L,ordering,k):
    count_fixed = 0
    S = [True for v in G.nodes]
    for v in ordering:
        if reachable_population(G, population, S, v) < L:
            #print("This is the root being zero:", v)
            for j in range(k):
                m._R[v,j].UB=0
                count_fixed += 1
        S[v]=False
    print("Number of LFixings =",count_fixed,"out of", G.number_of_nodes()*k)
    return count_fixed 