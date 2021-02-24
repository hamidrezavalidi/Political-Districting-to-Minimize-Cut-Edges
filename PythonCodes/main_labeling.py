from gurobipy import GRB 
import networkx as nx
import time
import fixing
import labeling
import separation
import matplotlib.pyplot as plt

###################
# labeling parameters
###################
# Note 1: For labeling base, set LFixing and UFixing to False





#############
# Apply active procedures
############# 
def build_base_labeling_model(m, G, population, L, U, k, heur_districts, ordered_vertices, position): 
    orbitope = False
    ZFixing = True     # fix Z[i,j,v]=0 when X[i,v]=0 fixed
    LFixing = True     # fix X[j,j]=0 if p(R_j)<L, where R_j is the vertices reachable from j in G[V_j], and V_j is vertices that DiagFixing allows
    UFixing = True     # fix X[i,j]=0 if total population 
    
    
    
    #bigM = 2  # in the SHIR model, what value(s) should be used for big-M?
            # 0: use the trivial value M = n-1
            # 1: use a smarter value M = max {|S| : p(S)<=U } - 1
            # 2: use smart values M_ij that depend on i and j.
            # FIXME the options 0,1,2 can all be improved by exploiting diagonal fixing!
    
    m._DG = nx.DiGraph(G) # bidirected version of G
    m._U = U
    m._k = k
    m._population = population   
    
    if m._extended==True:
        labeling.add_labeling_cut_edges_extended_objective(m, G, k) # Z[i,j,v]=1 if {i,j} cut because i->v but j!->v
    else:
        labeling.add_labeling_cut_edges_objective(m, G, k) # Y[i,j]=1 if {i,j} is cut 
     
    labeling.build_labeling_model(m, population, L, U, k)
    
    # LFixing
    if LFixing:
        L_fixed = fixing.do_labeling_LFixing(m,G,population,L,ordered_vertices,k)
    else:
        L_fixed = 'none'
        
    m._row.append(L_fixed)    
    m.update()
    
    # UFixing
    if UFixing:
        U_fixed = fixing.do_labeling_UFixing(m,G,population,U,ordered_vertices,k)    
    else:
        U_fixed = 'none'
    
    m._row.append(U_fixed)
    m.update()
    
    '''    
    # ZFixing
    if ZFixing:
        Z_fixed = fixing.do_labeling_ZFixing(m,G,k)    
        m._row.append(Z_fixed)
        m.update()
    '''    
    if LFixing and UFixing:
        XRFixed = L_fixed + U_fixed
        m._row += [XRFixed]
        out_of = 2*G.number_of_nodes()*k
        m._row += [out_of]
        print("Number of XRFixings =",XRFixed,"out of",out_of)
    else:
        m._row += ['none']
        out_of = 2*G.number_of_nodes()*k
        m._row += [out_of]
        print("Number of XRFixings =","none","out of",out_of)
    
    # if the model is not extended, then we have no ZFixing
    if not m._extended:
        ZFixing = False
    
    # fix some Z vars
    if m._extended and ZFixing:
        Z_fixed = fixing.do_labeling_ZFixing(m,G,k)
        m._row += [Z_fixed, k*G.number_of_edges()]
    else:
        m._row += [ '-', '-' ]
        
    # give a heuristic to Hess
    if m._use_heuristic:
        minList = []
        for district in heur_districts:
            positionList = []
            for i in district:
                positionList.append(position[i])
            a = min(positionList)
            minList.append(a)
        minList.sort()

        counter = -1
        for v in minList:
            counter += 1
            for district in heur_districts:
                for i in district:
                    if position[i] == v:
                        for j in district:
                            m._X[j,counter].start=1  
        
    # use the extended formulation for the partitioning orbitope?
    if orbitope==True:
        labeling.add_orbitope_extended_formulation(m, G, k, ordered_vertices, position)
    ############################
    # Run a connectivity model
    ############################ 
    
    DG = m._DG
    mip_start = time.time() 
    #print ("model: ", m._model)
    #m.params.OutputFlag = 0
    m.params.LogToConsole = 0
    m.params.LogFile='opt_'+m._state+'_'+m._model+'_default.log'
    m.Params.timeLimit=3600
    
    if m._model=="shir":
        labeling.build_shir_model(m,DG,k,ordered_vertices,orbitope)
        m.optimize()
    elif m._model=="scf":
        labeling.build_scf_model(m,G,DG,k,ordered_vertices,orbitope)
        m.optimize()    
    elif m._model=="mcf":
        labeling.build_mcf_model(m,DG,k,population,orbitope)
        m.optimize()
    elif m._model=="cut":
        m.Params.lazyConstraints = 1
        m.optimize(separation.labeling_cut_callback)
    elif m._model=="lcut":
        m.Params.lazyConstraints = 1
        m.optimize(separation.labeling_lcut_callback)
    elif m._model=="austin_cut":
        labeling.build_Austin_model(m,DG,k,ordered_vertices,orbitope)
        m.Params.lazyConstraints = 1
        m.optimize(separation.labeling_Austin_cut_callback) 
    elif m._model=="austin_lcut":
        labeling.build_Austin_model(m,DG,k,ordered_vertices,orbitope)
        m.Params.lazyConstraints = 1
        m.optimize(separation.labeling_Austin_lcut_callback)     
    elif m._model=="labeling":
        m.optimize()
    else:
        print("ERROR: there is no model called "+m._model+".")
    mip_stop = time.time()     
    print("mip time:", mip_stop - mip_start)    
    

    if m.status == GRB.OPTIMAL or m.status == GRB.TIME_LIMIT:
        m._row.append(m.NodeCount)
        m._row.append(m.objBound)
        m._row.append(m.objVal)
        # create a new column for district labels
        m._df['district']= -1 
            
        if m._all_optima==False:
            districts = [[i for i in DG.nodes if m._X[i,j].x > 0.5] for j in range(k)]
            print("districts=",districts)
            
            for d in range(len(districts)):
                for v in districts[d]:
                    geoID = G.node[v]["GEOID10"]
                    for u in G.nodes:
                        if geoID == m._df['GEOID10'][u]:
                            m._df['district'][u] = d
            
            # display the districting map 
            oplot = m._df.plot(column='district',figsize=(10, 10)).get_figure()
            plt.axis('off')
            oplot.savefig(m._state+"_"+m._model+"_opt.png")
            
        else:
            for sol in range(m.SolCount):
                m.Params.SolutionNumber=sol
                districts = [[i for i in DG.nodes if m._X[i,j].Xn > 0.5] for j in range(k)]
                print("districts=",districts)
                
                for d in range(len(districts)):
                    for v in districts[d]:
                        geoID = G.node[v]["GEOID10"]
                        for u in G.nodes:
                            if geoID == m._df['GEOID10'][u]:
                                m._df['district'][u] = d
                
                # display the districting map 
                oplot = m._df.plot(column='district',figsize=(10, 10), linewidth=1, edgecolor='black').get_figure()
                plt.axis('off')
                oplot.savefig(m._state+"_"+m._model+"_opt"+".png")
                
    else:
        m._row.append("not optimal")    
            
    m._row.append(round(mip_stop-mip_start,2))
    
    return m._row  
