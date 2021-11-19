import networkx as nx

def construct_position(ordering):
    position = [-1 for i in range(len(ordering))]
    for p in range(len(ordering)):
        v = ordering[p]
        position[v] = p
    return position


def create_orbitope_graph(ordering, k):
    
    OG = nx.DiGraph()
    n = len(ordering)    
    node_pos = dict()
    
    # add nodes s and t
    OG.add_node('s')
    node_pos['s'] = (0,1)

    OG.add_node('t')
    node_pos['t'] = ((k-1)/2,-n)
    
    # add internal nodes
    for p in range(n):
        i = ordering[p]
        for j in range(k):
            node = (i,j)
            OG.add_node(node)
            node_pos[node] = (j,-p)
            
    # edge out of s
    OG.add_edge('s',(ordering[0],0))
    
    # edges into t
    last = ordering[n-1]
    for j in range(k):
        OG.add_edge((last,j),'t')
        
    # internal edges
    for p in range(n-1):
        i = ordering[p]
        ni = ordering[p+1]

        for j in range(k-1):
            OG.add_edge((i,j),(ni,j))
            OG.add_edge((i,j),(ni,j+1))

        OG.add_edge((i,k-1),(ni,k-1))
        
    nx.draw(OG, pos=node_pos, with_labels=True)
    
    return OG


def add_shifted_column_inequalities(OG, xval, ordering, k):
    
    n = len(ordering)
    
    # Get x(bar_ij) for each leader (i,j).
    # See Corollary 11 of Kaibel and Pfetsch (MP, 2008).
    # Or bottom of page 688 of Faenza and Kaibel (MOR, 2009).
    #
    xbar = xval.copy()
    for j in range(k-2,-1,-1):
        for i in range(n):
            xbar[i,j] += xbar[i,j+1]
            
    position = construct_position(ordering)
    
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
            OG.edges[pnode,node]['edge_weight'] = xval[i,j] # vertical arc
            
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

                if xbar[next_i,j+1] > S_weight:
                    #print("distance from",path_start,"to (i,j) =",node,"is",distance[node])
                    #print("distance from",path_start,"to (i,j) =",node,"plus x[path_start] is",distance[node]+xval[path_start])
                    #print("my_path",my_path)
                    
                    B = [ (next_i,jval) for jval in range(j+1,k) ]
                    B_weight = sum( xval[myi,myj] for (myi,myj) in B )
                    
                    print("S_weight =",S_weight)
                    print("B_weight =",B_weight,"=",xbar[next_i,j+1])
                    print("The bar leader is (",next_i,",",j+1,"), meaning that B =",B)

                    print("The inequality to add is x(B) <= x(S) where B =",B,"and S =",S)
                    print("It is violated by this much",xbar[next_i,j+1]-S_weight)
    
    return
    