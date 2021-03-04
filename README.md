# Political Districting to Minimize Cut Edges

We consider a stylized redistricting problem.

The input is graph G=(V,E), integer k, population vector p, and bounds L and U.

The task is to find a partition of V into k subsets V_1, V_2, ..., V_k such that:
1. each G[V_i] is connected, 
2. each V_i has population within [L,U], and 
3. the number of edges between partitions is minimized.

For this problem, we consider multiple mixed integer programming (MIP) formulations. We solve them with the Gurobi solver. The code uses lots of "tricks" to speed up the computations (e.g., strong extended formulation for cut edges objective, symmetry handling techniques, safe variable fixing rules, different approaches for imposing contiguity constraints, heuristic warm start with GerryChain, etc).

The input data is duplicated from Daryl DeFord's website: 
https://people.csail.mit.edu/ddeford/dual_graphs.html

The shape files used to draw maps are duplicated from Eugene Lykhovyd's website: 
https://lykhovyd.com/files/public/districting
