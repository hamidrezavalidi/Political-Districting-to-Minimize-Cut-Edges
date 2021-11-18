import networkx as nx

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
    