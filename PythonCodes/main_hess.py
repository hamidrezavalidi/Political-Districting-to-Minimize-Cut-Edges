from gurobipy import GRB 
import networkx as nx
import time
import fixing
import hess
import separation
import matplotlib.pyplot as plt


#############
#Apply active procedures
############# 
def build_base_hess_model(m, G, population, L, U, k, heur_districts, ordered_vertices, position):
    #print("Hello Hamid")
    m._DG = nx.DiGraph(G) # bidirected version of G
    m._U = U
    m._k = k
    m._population = population
    
    ###################
    #Hess parameters
    ###################
    DiagFixing = True  # given a vertex ordering (v1, v2, ..., vn) fix X[i,j]=0 when i comes before j in ordering
    LFixing = True     # fix X[j,j]=0 if p(R_j)<L, where R_j is the vertices reachable from j in G[V_j], and V_j is vertices that DiagFixing allows
    UFixing = True     # fix X[i,j]=0 if i and j are far apart from each other. Do not use this for base models Hess and labeling.
    ZFixing = True     # fix Z[i,j,v]=0 when X[i,v]=0 fixed
    bigM = 2  # in the SHIR model, what value(s) should be used for big-M?
            # 0: use the trivial value M = n-1
            # 1: use a smarter value M = max {|S| : p(S)<=U } - 1
            # 2: use smart values M_ij that depend on i and j.
            # FIXME the options 0,1,2 can all be improved by exploiting diagonal fixing!
    
    if m._extended==True:
        hess.add_cut_edges_extended_objective(m, G) # Z[i,j,v]=1 if {i,j} cut because i->v but j!->v
    else:
        hess.add_cut_edges_objective(m, G) # Y[i,j]=1 if {i,j} is cut 
     
    hess.build_hess_model(m, population, L, U, k)                


    
    # fix some X variables 
    XFixed = 0
    
    if DiagFixing:
        DiagFixed = fixing.do_DiagFixing(m,G,position)
        XFixed += DiagFixed
    else:
        DiagFixed = 'none'    
        
    m._row.append(DiagFixed)    
    
    if LFixing:
        LFixed = fixing.do_LFixing(m,G,population,L,ordered_vertices,position)
        XFixed += LFixed
    else:
        LFixed = 'none'
    
    m._row.append(LFixed)     
        
    m.update()
    
    if UFixing:
        UFixed = fixing.do_UFixing(m,population,U,ordered_vertices,position)
        XFixed += UFixed
    else:
        UFixed = 'none'
    
    m._row.append(UFixed)
    
    m.update()
        
    if DiagFixing or LFixing or UFixing: 
        #XFixed = DiagFixed + LFixed + UFixed
        m._row.append(XFixed)
        out_of = G.number_of_nodes()*G.number_of_nodes()
        m._row.append(out_of)
        print("Number of XFixings =",XFixed,"out of",out_of)
    else:
        m._row += ['-','-']
        
    # if the model is not extended, then we have no ZFixing
    if not m._extended:
        ZFixing = False    
    
    # fix some Z vars
    if m._extended and ZFixing:
        ZFixed = fixing.do_ZFixing(m._X,m._Z,G)
        m._row += [ZFixed, G.number_of_nodes()*G.number_of_edges()]
    else:
        m._row += [ '-', '-' ]
    
    '''       
    # display the remaining possible centers
    m._df['centers']= -1 
    for j in G.nodes:
        if m._X[j,j].UB>0.5:
            m._df['centers'][j] = 0
        else:
            m._df['centers'][j] = 1
    m._df.plot(column='centers') 
    plt.axis('off')
    '''
    # give a heuristic to Hess
    if m._use_heuristic:
        centers = []

        for district in heur_districts:    
            min_pos = min([position[v] for v in district])
            center = ordered_vertices[min_pos]
            centers.append(center)
            for i in district:
                m._X[i,center].start=1
        '''            
        labels = {}       
        label=0
        for i in ordered_vertices:
            if i in centers:
                labels[i] = label
                label += 1        
                    
        # create a new column for heuristic district labels, and fill it in
        m._df['heuristic']= -1 
        for district in heur_districts:
            min_pos = min([position[v] for v in district])
            center = ordered_vertices[min_pos]
            for v in district:
                m._df['heuristic'][v] = labels[center]         
        '''
    ############################
    # Run a connectivity model
    ############################  
    #X = m._X
    DG = m._DG
    mip_start = time.time()  
    #m.params.OutputFlag = 0
    m.params.LogToConsole = 0
    m.params.LogFile='opt_'+m._state+'_'+m._model+'.log' 
    m.Params.timeLimit=3600
    # add contiguity constraints (whichever type is requested) and solve MIP
    if m._model=="shir":
        hess.build_shir_model(m,DG,population,U,bigM)
        m.optimize()
    elif m._model=="mcf":
        hess.build_mcf_model(m,DG,population)
        m.optimize()
    elif m._model=="scf":
        hess.build_scf_model(m,G,DG,population)
        m.optimize() 
    elif m._model=="williams":
        hess.build_williams_model(m,G,population,k)
        m.optimize()    
    elif m._model=="cut":
        m.Params.lazyConstraints = 1
        m.optimize(separation.hess_cut_callback)
    elif m._model=="lcut":
        m.Params.lazyConstraints = 1
        m.optimize(separation.hess_lcut_callback)
    elif m._model=="hess":
        #m.Params.timeLimit=3600
        #m.Params.method=2
        m.optimize() 
    else:
        print("ERROR: there is no model called "+m._model+".")     
    mip_stop = time.time()     
    print("mip time:", mip_stop - mip_start)
    #m._row.append(round(mip_stop-mip_start,2))

    #############################
    # Draw Solutions
    #############################
    if m.status == GRB.OPTIMAL or m.status == GRB.TIME_LIMIT:
        m._row.append(m.NodeCount)
        m._row.append(m.objBound)
        m._row.append(m.objVal)
        # create a new column for district labels
        m._df['district']= -1 
        
        if m._all_optima==False:
            centers = []
            for j in ordered_vertices:
                if m._X[j,j].x>0.5:
                    centers.append(j)
            districts = [[i for i in DG.nodes if m._X[i,j].x > 0.5] for j in centers]
            print("districts=",districts)
        
            #print("data frame: ", m._df)
            for d in range(len(districts)):
                for v in districts[d]:
                    geoID = G.node[v]["GEOID10"]
                    #print("geoID: ", geoID)
                    for u in G.nodes:
                        if geoID == m._df['GEOID10'][u]:
                            #print("Here we got it! It is ", u)
                            m._df['district'][u] = d
        
            # display the districting map 
            oplot = m._df.plot(column='district',figsize=(10, 10)).get_figure()
            plt.axis('off')
            if m._weighted == True:
                oplot.savefig(m._state+"_"+m._model+"_opt_weighted.png")
            else:    
                oplot.savefig(m._state+"_"+m._model+"_opt.png")
        
        else:
            for sol in range(m.SolCount):
                m.Params.SolutionNumber=sol
                centers = []
                for j in ordered_vertices:
                    if m._X[j,j].Xn>0.5:
                        centers.append(j)
                districts = [[i for i in DG.nodes if m._X[i,j].Xn > 0.5] for j in centers]
                print("districts=",districts)
            
                
                print("data frame: ", m._df)
                for d in range(len(districts)):
                    for v in districts[d]:
                        geoID = G.node[v]["GEOID10"]
                        print("geoID: ", geoID)
                        for u in G.nodes:
                            if geoID == m._df['GEOID10'][u]:
                                print("Here we got it! It is ", u)
                                m._df['district'][u] = d
            
                # display the districting map 
                oplot = m._df.plot(column='district',figsize=(10, 10), linewidth=1, edgecolor='black').get_figure()
                plt.axis('off')
                #oplot = m._df.plot(column='district').get_figure()
                oplot.savefig(m._state+"_"+m._model+"_opt"+sol+".png")
                
             
    else:
        m._row.append("not optimal")    
        
    m._row.append(round(mip_stop-mip_start,2))
    m.update()
    return m._row            
