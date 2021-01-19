import gurobipy as gp
from gurobipy import GRB 

import bisect
import networkx as nx

def add_cut_edges_objective(m, G):
    # Y[i,j] = 1 if edge {i,j} is cut
    m._Y = m.addVars(G.edges, vtype=GRB.BINARY)
    for i,j in G.edges:
        m.addConstrs( m._X[i,v]-m._X[j,v] <= m._Y[i,j] for v in G.nodes)
        #m.addConstrs( m._X[j,v]-m._X[i,v] <= m._Y[i,j] for v in G.nodes)
    m.setObjective( gp.quicksum(G[i][j]['edge_length']*m._Y[i,j] for i,j in G.edges), GRB.MINIMIZE )
    

def add_cut_edges_extended_objective(m, G):
    # Z[i,j,v] = 1 if edge (i,j) is cut because i->v but j!->v
    m._Z = m.addVars(G.edges, G.nodes, vtype=GRB.BINARY) #)CONTINUOUS)
    m.addConstrs( m._X[i,v]-m._X[j,v] <= m._Z[i,j,v] for i,j in G.edges for v in G.nodes)
    m.setObjective( gp.quicksum(G[i][j]['edge_length']*m._Z[i,j,v] for i,j in G.edges for v in G.nodes), GRB.MINIMIZE )
        
def most_possible_nodes_in_one_district(DG, population, U):
    cumulative_population = 0
    num_nodes=0
    for ipopulation in sorted(population):
        cumulative_population += ipopulation
        num_nodes+=1
        if cumulative_population>U:
            return num_nodes-1

def compute_bigM_values(DG,population,U,bigM):
    if bigM==0:
        print("Using M=",DG.number_of_nodes()-1)
        return [[DG.number_of_nodes()-1 for i in DG.nodes] for j in DG.nodes]
    elif bigM==1:
        my_bigM = most_possible_nodes_in_one_district(DG,population,U) - 1
        print("Using M=",my_bigM)
        return [[my_bigM for i in DG.nodes] for j in DG.nodes]
    elif bigM!=2:
        print("ERROR: big-M parameter=",bigM,"is not supported.")
        return
    
    # Now, bigM==2.   COMPUTE SMALL VALUES FOR BIG-M, WHICH DEPEND ON i AND j
    M = [[0 for i in DG.nodes] for j in DG.nodes]
    for (i,j) in DG.edges:
        DG[i][j]['weight']=population[j] # weight of edge (i,j) is population of its head j
    
    cumulative_population = [0 for i in range(len(population))]
    t=0
    for ipopulation in sorted(population):
        if t==0:
            cumulative_population[t] = ipopulation
        else:
            cumulative_population[t] = cumulative_population[t-1]+ipopulation
        t += 1
    
    for j in DG.nodes:
        dist = nx.shortest_path_length(DG,source=j,weight='weight')
        for i in DG.nodes:
            R = int(U - ( dist[i] + population[j] ))
            if R<0: # infeasible
                M[i][j] = 0
            else:   # how many nodes could be reached after getting to i from j? (then +1 for zero-based, and another +1 for i)
                M[i][j] = (bisect.bisect(cumulative_population, R)+1) 

    return M

def build_hess_model(m, population, L, U, k):
    DG = m._DG
    # Each vertex i assigned to one district
    m.addConstrs(gp.quicksum(m._X[i,j] for j in DG.nodes) == 1 for i in DG.nodes)
     
    # Pick k centers
    m.addConstr(gp.quicksum(m._X[j,j] for j in DG.nodes) == k)
    
    # Population balance: population assigned to vertex j should be in [L,U]
    m.addConstrs(gp.quicksum(population[i] * m._X[i,j] for i in DG.nodes) <= U * m._X[j,j] for j in DG.nodes)
    m.addConstrs(gp.quicksum(population[i] * m._X[i,j] for i in DG.nodes) >= L * m._X[j,j] for j in DG.nodes)
    
    # Add coupling inequalities for added model strength
    m.addConstrs(m._X[i,j] <= m._X[j,j] for i in DG.nodes for j in DG.nodes)
    #couplingConstr = m.addConstrs(m._X[i,j] <= m._X[j,j] for i in DG.nodes for j in DG.nodes)
    
    # Madke them user cuts
    #for i in DG.nodes:
     #   for j in DG.nodes:
      #      couplingConstr[i,j].Lazy = -1
    
    #m.update()
    # Set branch priority on center vars
    for j in DG.nodes:
        m._X[j,j].BranchPriority=1    
        
        
def build_scf_model(m,G,DG,population):
    # F[j,u,v] tells how much flow (from source j) is sent across arc (u,v)
    F = m.addVars(DG.edges,vtype=GRB.CONTINUOUS)
    #M = compute_bigM_values(DG,population,U,bigM)
    for j in DG.nodes:
        m.addConstr( gp.quicksum(m._X[i,j] for i in DG.nodes) + gp.quicksum(F[u,j] for u in DG.neighbors(j)) == gp.quicksum(F[j,u] for u in DG.neighbors(j)) + 1)
        m.addConstr( gp.quicksum(F[u,j] for u in DG.neighbors(j)) <= (len(G.nodes) - m._k + 1)*(1-m._X[j,j]))
        
    for (i,j) in G.edges:
                if m._extended==True:
                    m.addConstr( F[i,j] + F[j,i] <= (len(G.nodes) - m._k + 1)*(1 - gp.quicksum(m._Z[i,j,v] for v in G.nodes)) )
                else:
                    m.addConstr( F[i,j] + F[j,i] <= (len(G.nodes) - m._k + 1)*(1 - m._Y[i,j]) )

def build_shir_model(m,DG,population,U,bigM):
    # F[j,u,v] tells how much flow (from source j) is sent across arc (u,v)
    F = m.addVars(DG.nodes,DG.edges,vtype=GRB.CONTINUOUS)
    M = compute_bigM_values(DG,population,U,bigM)
    for j in DG.nodes:
        m.addConstr( gp.quicksum(F[j,u,j] for u in DG.neighbors(j)) == 0 )
        for i in DG.nodes:
            if i!=j:
                m.addConstr( gp.quicksum(F[j,u,i] for u in DG.neighbors(i)) - gp.quicksum(F[j,i,u] for u in DG.neighbors(i)) == m._X[i,j] )
                m.addConstr( gp.quicksum(F[j,u,i] for u in DG.neighbors(i)) <= M[i][j]*m._X[i,j] )
    #m.write("shir.lp")
    m.update()
    
def read_dual_graph(file_path, state):   
    #f = open(filename, "r")  
    file = open(file_path+state+"_dual.txt")
    all_the_lines = file.readlines()
    edges = []
    for line in all_the_lines:
        int_list = [int(i) for i in line.split()]
        edges.append(int_list)          
    return edges

def build_williams_model(m,G,population,k):
    file_path = "C:\\data\\DualGraphs\\"
    dual_edges = read_dual_graph(file_path, m._state)
        
    D = nx.read_edgelist(file_path+m._state+"_dual.txt",nodetype=int)
    
    DD = nx.DiGraph(D)
    
    #print(DD.edges)
    
    DG = nx.DiGraph(G)
    
    W = m.addVars(DG.edges,vtype=GRB.CONTINUOUS)
    
    U = m.addVars(DG.edges,vtype=GRB.CONTINUOUS)
    
    #S = m.addVars(DG.edges,vtype=GRB.CONTINUOUS)
    
    #R = m.addVars(DG.nodes,vtype=GRB.CONTINUOUS)
    
    W_star = m.addVars(DD.edges,vtype=GRB.CONTINUOUS)
    for i, j in G.edges:
        r = i
        break

    #print("r",r)        

    r_star = dual_edges[0][0]


    #incomming constraints for primal graph
    for j in DG.nodes:
        if j != r:
            m.addConstr(sum(W[u,j] for u in DG.neighbors(j)) == 1)
        else:
            m.addConstr(sum(W[u,r] for u in DG.neighbors(r)) == 0)    

    #incomming constraints for dual graph
    for j in DD.nodes:
        if j != r_star:
            m.addConstr(sum(W_star[u,j] for u in DD.neighbors(j)) == 1)
        else:
            m.addConstr(sum(W_star[u,r_star] for u in DD.neighbors(r_star)) == 0)


    #crossing edges constraints
    counter = -1
    for i,j in G.edges:
        counter += 1
        u = dual_edges[counter][0]
        v = dual_edges[counter][1]
        m.addConstr(W[i, j] + W[j, i] + W_star[u, v] + W_star[v, u] == 1)    

    #sub_tree edge constraints
    for i,j in G.edges:
        m.addConstr(U[i, j] + U[j, i] <= W[i, j] + W[j, i])
        
    # sum of edges is n-k
    #m.addConstr(gp.quicksum(U[i, j] for i,j in G.edges) == len(G.nodes) - k)

    # coupling with cut-edges
    if m._extended:
        m.addConstrs(U[i, j] + U[j, i] <= 1 - gp.quicksum(m._Z[i,j,v] for v in G.nodes) for (i,j) in G.edges)
    else:    
        m.addConstrs(U[i, j] + U[j, i] <= 1 - m._Y[i,j] for (i,j) in G.edges)


    #either root or incomming
    for j in DG.nodes:
        m.addConstr(gp.quicksum(U[v,j] for v in DG.neighbors(j)) + m._X[j,j] == 1)

    #number of roots
    #m.addConstr(sum(R[j] for j in DG.nodes) == k)
        
    #faces = [(0, 66, 57), (66, 0, 61), (0, 57, 73), (57, 0, 66), (0, 61, 66), (61, 0, 27), (0, 73, 27), (73, 0, 57), (0, 27, 61), (27, 0, 73), (1, 51, 45), (51, 1, 22), (1, 22, 51), (22, 1, 56), (1, 45, 68), (45, 1, 51), (1, 56, 22), (56, 1, 53), (1, 68, 46), (68, 1, 45), (1, 53, 56), (53, 1, 46), (1, 46, 53), (46, 1, 68), (2, 63, 26), (63, 2, 36), (26, 2, 63), (2, 36, 63), (36, 2, 31), (2, 31, 36), (31, 2, 43), (2, 43, 31), (32, 3, 71), (3, 37, 72), (3, 71, 32), (71, 3, 72), (3, 72, 71), (72, 3, 37), (4, 40, 22), (40, 4, 28), (4, 28, 40), (28, 4, 24), (4, 22, 24), (22, 4, 40), (4, 24, 28), (24, 4, 22), (5, 72, 25), (72, 5, 20), (5, 25, 47), (25, 5, 72), (5, 20, 72), (20, 5, 58), (5, 47, 58), (47, 5, 25), (5, 58, 20), (58, 5, 47), (6, 20, 58), (20, 6, 75), (6, 58, 19), (58, 6, 20), (6, 75, 20), (75, 6, 19), (6, 19, 75), (19, 6, 58), (7, 71, 72), (71, 7, 65), (7, 72, 20), (72, 7, 71), (7, 65, 71), (65, 7, 20), (7, 20, 65), (20, 7, 72), (8, 32, 71), (8, 71, 35), (71, 8, 32), (52, 8, 35), (8, 35, 52), (35, 8, 71), (9, 52, 51), (40, 9, 22), (9, 51, 22), (51, 9, 52), (9, 22, 40), (22, 9, 51), (10, 47, 60), (47, 10, 58), (10, 58, 47), (58, 10, 19), (10, 60, 39), (60, 10, 47), (10, 19, 58), (19, 10, 69), (10, 39, 69), (39, 10, 60), (10, 69, 19), (69, 10, 39), (37, 11, 25), (11, 76, 34), (11, 25, 37), (25, 11, 34), (11, 34, 25), (34, 11, 76), (12, 42, 48), (42, 12, 33), (12, 48, 70), (48, 12, 42), (70, 12, 48), (12, 33, 42), (13, 54, 38), (54, 13, 30), (13, 30, 54), (38, 13, 54), (14, 53, 46), (53, 14, 21), (14, 21, 53), (21, 14, 42), (14, 46, 29), (46, 14, 53), (14, 42, 21), (42, 14, 48), (14, 29, 48), (29, 14, 46), (14, 48, 42), (48, 14, 29), (15, 67, 75), (67, 15, 68), (15, 75, 66), (75, 15, 67), (15, 68, 67), (68, 15, 61, 50), (15, 66, 61), (66, 15, 75), (15, 61, 50, 68), (61, 15, 66), (76, 16, 34), (16, 34, 76), (34, 16, 55), (16, 55, 34), (55, 16, 60), (16, 60, 55), (60, 16, 59), (16, 59, 60), (17, 74, 31), (74, 17, 49), (17, 49, 74), (31, 17, 74), (19, 58, 10), (58, 19, 6), (19, 75, 6), (75, 19, 66), (19, 69, 57), (69, 19, 10), (19, 66, 75), (66, 19, 57), (19, 57, 66), (57, 19, 69), (20, 72, 5), (72, 20, 7), (20, 65, 7), (65, 20, 44), (20, 58, 6), (58, 20, 5), (20, 44, 65), (44, 20, 75), (20, 75, 44), (75, 20, 6), (21, 24, 53), (21, 53, 14), (53, 21, 24), (42, 21, 14), (22, 40, 9), (40, 22, 4), (22, 51, 1), (51, 22, 9), (22, 56, 24), (56, 22, 1), (22, 24, 4), (24, 22, 56), (23, 69, 64), (69, 23, 57), (23, 57, 69), (57, 23, 73), (23, 64, 74), (64, 23, 69), (23, 74, 49), (74, 23, 64), (23, 73, 57), (73, 23, 62), (23, 62, 73), (62, 23, 49), (23, 49, 62), (49, 23, 74), (24, 28, 4), (24, 56, 53), (56, 24, 22), (24, 53, 21), (53, 24, 56), (25, 37, 11), (37, 25, 72), (25, 72, 37), (72, 25, 5), (25, 34, 47), (34, 25, 11), (25, 47, 5), (47, 25, 34), (59, 26, 63), (26, 63, 59), (63, 26, 2), (27, 61, 0), (61, 27, 54), (27, 73, 62), (73, 27, 0), (27, 54, 61), (54, 27, 38), (27, 62, 38), (62, 27, 73), (27, 38, 54), (38, 27, 62), (28, 40, 4), (29, 46, 50), (46, 29, 14), (29, 50, 54), (50, 29, 46), (29, 54, 30), (54, 29, 50), (29, 48, 14), (48, 29, 70), (29, 70, 48), (70, 29, 30), (29, 30, 70), (30, 29, 54), (30, 54, 13), (54, 30, 29), (30, 70, 29), (31, 36, 2), (36, 31, 74), (31, 74, 36), (74, 31, 17), (43, 31, 2), (32, 71, 8), (71, 32, 3), (33, 42, 12), (34, 76, 16), (76, 34, 11), (34, 55, 47), (55, 34, 16), (34, 47, 25), (47, 34, 55), (35, 71, 65), (71, 35, 8), (35, 52, 8), (52, 35, 51), (35, 51, 52), (51, 35, 45), (35, 65, 45), (65, 35, 71), (35, 45, 51), (45, 35, 65), (36, 63, 2), (63, 36, 64), (36, 64, 63), (64, 36, 74), (36, 74, 64), (74, 36, 31), (37, 72, 3), (72, 37, 25), (38, 54, 27), (54, 38, 13), (62, 38, 27), (39, 60, 59), (60, 39, 10), (39, 59, 63), (59, 39, 60), (39, 69, 10), (69, 39, 64), (39, 63, 64), (63, 39, 59), (39, 64, 69), (64, 39, 63), (42, 48, 12), (48, 42, 14), (44, 65, 20), (65, 44, 45), (44, 45, 65), (45, 44, 68), (44, 67, 68), (67, 44, 75), (44, 75, 67), (75, 44, 20), (44, 68, 45), (68, 44, 67), (45, 51, 35), (51, 45, 1), (45, 65, 44), (65, 45, 35), (45, 68, 1), (68, 45, 44), (46, 68, 50), (68, 46, 1), (46, 53, 1), (53, 46, 14), (46, 50, 29), (50, 46, 68), (47, 55, 60), (55, 47, 34), (47, 58, 5), (58, 47, 10), (47, 60, 10), (60, 47, 55), (48, 70, 12), (70, 48, 29), (49, 74, 17), (74, 49, 23), (49, 62, 23), (50, 68, 15, 61), (68, 50, 46), (50, 61, 54), (61, 50, 68, 15), (50, 54, 29), (54, 50, 61), (51, 52, 35), (52, 51, 9), (53, 56, 1), (56, 53, 24), (54, 61, 27), (61, 54, 50), (55, 60, 47), (60, 55, 16), (57, 69, 23), (69, 57, 19), (57, 66, 19), (66, 57, 0), (57, 73, 0), (73, 57, 23), (59, 60, 16), (60, 59, 39), (59, 63, 39), (63, 59, 26), (61, 66, 0), (66, 61, 15), (62, 73, 23), (73, 62, 27), (63, 64, 39), (64, 63, 36), (64, 69, 39), (69, 64, 23), (64, 74, 23), (74, 64, 36), (65, 71, 7), (71, 65, 35), (66, 75, 19), (75, 66, 15), (67, 75, 15), (75, 67, 44), (67, 68, 44), (68, 67, 15), (71, 72, 7), (72, 71, 3)]

    #edges_of_faces = []
    
    #print("Edges: ",G.edges)
    '''
    for face in faces:
        new_face=[]
        for i in face:
            for j in face:
                if i<j and (i,j) in G.edges:
                    new_face.append((i,j))
        edges_of_faces.append(new_face)            
    
    for face in edges_of_faces:
        for (i,j) in face:
            for (u,v) in face:
                if (i,j) != (u,v):
                    if (i,j) in G.edges and (u,v) in G.edges:
                        m.addConstr(S[i, j] + S[j, i] + S[u, v] + S[v, u] <= 1 + gp.quicksum(m._Z[i,j,p] for p in G.nodes))
                        m.addConstr(S[i, j] + S[j, i] + S[u, v] + S[v, u] <= 1 + gp.quicksum(m._Z[u,v,p] for p in G.nodes))
    '''                         
                    
                    
    
def build_mcf_model(m,DG,population):    
    commodities = [(a,b) for a in DG.nodes for b in DG.nodes if m._X[a,b].UB>0.5 and a!=b]
    print("# commodities=",len(commodities),"out of",len(DG.nodes)*len(DG.nodes))
    
    # F[a,b,u,v] tells how much flow (from b to a) is sent across arc (u,v)
    F = m.addVars(commodities,DG.edges,vtype=GRB.CONTINUOUS)
    for a,b in commodities:
        m.addConstr( m._X[a,b] == gp.quicksum(F[a,b,b,u] for u in DG.neighbors(b)) - gp.quicksum(F[a,b,u,b] for u in DG.neighbors(b)) )
        m.addConstr( 0 == gp.quicksum(F[a,b,u,b] for u in DG.neighbors(b)) )
        m.addConstrs( 0 == gp.quicksum(F[a,b,i,u] for u in DG.neighbors(i)) - gp.quicksum(F[a,b,u,i] for u in DG.neighbors(i)) for i in DG.nodes if i!=a and i!=b)
        m.addConstrs( gp.quicksum(F[a,b,u,j] for u in DG.neighbors(j)) <= m._X[j,b] for j in DG.nodes if j!=b )
