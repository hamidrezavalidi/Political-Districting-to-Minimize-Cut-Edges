# Political Districting to Minimize Cut Edges

Python code for the paper "Political districting to minimize cut edges" by Hamidreza Validi and Austin Buchanan.

We consider a stylized redistricting problem.

The input is graph G=(V,E), integer k, population vector p, and bounds L and U.

The task is to find a partition of V into k subsets V_1, V_2, ..., V_k such that:
1. each G[V_i] is connected, 
2. each V_i has population within [L,U], and 
3. the number of edges between partitions is minimized.

For this problem, we consider multiple mixed integer programming (MIP) formulations. We solve them with the Gurobi solver. The code uses lots of "tricks" to speed up the computations (e.g., strong extended formulation for cut edges objective, symmetry handling techniques, safe variable fixing rules, different approaches for imposing contiguity constraints, heuristic warm start with GerryChain, etc).

To run the code, you will need installations of [Gurobi](https://www.gurobi.com/) and [GerryChain](https://gerrychain.readthedocs.io/en/latest/).

The input data is duplicated from Daryl DeFord's website: 
https://people.csail.mit.edu/ddeford/dual_graphs.html

The shape files used to draw maps are duplicated from Eugene Lykhovyd's website: 
https://lykhovyd.com/files/public/districting

You can run the code from command line, like this:

  python main.py config.json

The config file can specify a batch of runs. A particular run might look like this:
* state: OK
* level: county
* base: hess
* contiguity: scf
* symmetry: default
* extended: true
* order: B_decreasing
* heuristic: true
* lp: true

Generally, each run should pick from the following options:
* state : {AL, AK, AZ, AR, CA, ...} [see link](https://en.wikipedia.org/wiki/List_of_U.S._state_and_territory_abbreviations)
* level : {county, tract}
* base : {hess, labeling},
* contiguity : {none, lcut, scf, shir},
* symmetry : {default, aggressive, orbitope}
* extended : {True, False}
* order : {none, decreasing, B_decreasing}
* heuristic : {True, False}
* lp : {True, False} 

