# Political Districting to Minimize Cut Edges

Python code for the paper "Political districting to minimize cut edges" by Hamidreza Validi and Austin Buchanan.

We consider a stylized redistricting problem.

The input is graph G=(V,E), integer k, population vector p, and bounds L and U.

The task is to find a partition of V into k subsets V_1, V_2, ..., V_k such that:
1. each G[V_i] is connected, 
2. each V_i has population within [L,U], and 
3. the number of edges between partitions is minimized.

For this problem, we consider multiple mixed integer programming (MIP) formulations. We solve them with the Gurobi solver. The code uses lots of "tricks" to speed up the computations (e.g., strong extended formulation for cut edges objective, symmetry handling techniques, safe variable fixing rules, different approaches for imposing contiguity constraints, heuristic warm start with GerryChain, etc).

## Require
To run the code, you will need installations of [Gurobi](https://www.gurobi.com/) and [GerryChain](https://gerrychain.readthedocs.io/en/latest/).

The input data is duplicated from [Daryl DeFord's website](https://people.csail.mit.edu/ddeford/dual_graphs.html).

The shape files used to draw maps are duplicated from [Eugene Lykhovyd's website](https://lykhovyd.com/files/public/districting).

## Run
You can run the code from command line, like this:

```
python main.py config.json
```

## config.json
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

The config.json file might look like this:
```
{
  "run1": {"state": "ME", "level": "county", "base": "hess", "contiguity": "scf", "symmetry": "default", "extended": true, "order": "B_decreasing", "heuristic": false, "lp": true},
  "run2": {"state": "NM", "level": "county", "base": "hess", "contiguity": "scf", "symmetry": "default", "extended": true, "order": "B_decreasing", "heuristic": true, "lp": true},
  "run3": {"state": "ID", "level": "county", "base": "hess", "contiguity": "scf", "symmetry": "default", "extended": true, "order": "B_decreasing", "heuristic": true, "lp": true},
  "run4": {"state": "WV", "level": "county", "base": "hess", "contiguity": "scf", "symmetry": "default", "extended": true, "order": "B_decreasing", "heuristic": true, "lp": true},
  "run5": {"state": "LA", "level": "county", "base": "hess", "contiguity": "scf", "symmetry": "default", "extended": true, "order": "B_decreasing", "heuristic": false, "lp": true},
  "run6": {"state": "AL", "level": "county", "base": "hess", "contiguity": "scf", "symmetry": "default", "extended": true, "order": "B_decreasing", "heuristic": true, "lp": true},
  "run7": {"state": "AR", "level": "county", "base": "hess", "contiguity": "scf", "symmetry": "default", "extended": true, "order": "B_decreasing", "heuristic": true, "lp": true},
  "run8": {"state": "OK", "level": "county", "base": "hess", "contiguity": "scf", "symmetry": "default", "extended": true, "order": "B_decreasing", "heuristic": true, "lp": true},
  "run9": {"state": "MS", "level": "county", "base": "hess", "contiguity": "scf", "symmetry": "default", "extended": true, "order": "B_decreasing", "heuristic": true, "lp": true},
  "run10": {"state": "NE", "level": "county", "base": "hess", "contiguity": "scf", "symmetry": "default", "extended": true, "order": "B_decreasing", "heuristic": true, "lp": true},
  "run11": {"state": "IA", "level": "county", "base": "hess", "contiguity": "scf", "symmetry": "default", "extended": true, "order": "B_decreasing", "heuristic": true, "lp": true},
  "run12": {"state": "KS", "level": "county", "base": "hess", "contiguity": "scf", "symmetry": "default", "extended": true, "order": "B_decreasing", "heuristic": true, "lp": true}
}
```

## Config options
Generally, each run should pick from the following options:
* state : {AL, AK, AZ, AR, CA, ... } 
  * [See link for list of 2-letter codes](https://en.wikipedia.org/wiki/List_of_U.S._state_and_territory_abbreviations)
* level : {county, tract}
  * Either treat counties or census tracts as indivisible land units
* base : {hess, labeling} 
  * Hess model uses binary variables x_ij that equal one when vertex i is assigned to the district rooted at vertex j
  * Labeling model uses binary variables x_ij that equal one when vertex i is assigned to district number j, where j in 1..k
* contiguity : {none, lcut, scf, shir}
  * none means that contiguity is not imposed
  * LCUT imposes contiguity with length-U a,b-separator inequalities (in branch-and-cut fashion)
  * SCF imposes contiguity with a single-commodity flow model. See [Hojny et al](https://link.springer.com/article/10.1007/s12532-020-00186-3)
  * SHIR imposes contiguity with a multi-commodity flow model. See [Shirabe2005](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1538-4632.2005.00605.x) and [Shirabe2009](https://journals.sagepub.com/doi/abs/10.1068/b34104) and [Oehrlein and Haunert](http://www.josis.org/index.php/josis/article/viewArticle/379) and [Validi et al.](http://www.optimization-online.org/DB_HTML/2020/01/7582.html).
* symmetry : {default, aggressive, orbitope}
  * Default uses whatever the Gurobi MIP solver does by default
  * Aggressive is a Gurobi setting that seeks to exploit symmetry
  * Orbitope is the extended formulation for partitioning orbitopes due to [Faenza and Kaibel](https://pubsonline.informs.org/doi/abs/10.1287/moor.1090.0392)
* extended : {True, False}
  * If true, use an extended formluation to better capture the cut edges objective function. See [Ferreira et al](https://link.springer.com/article/10.1007/BF02592198)
* order : {none, decreasing, B_decreasing}
  * If none, the given vertex ordering will be used
  * If decreasing, the vertices will be sorted in terms of decreasing population
  * If B_decreasing, a vertex subset B in which all components of G[B] have population less than L will be placed at back, others placed at front by decreasing population
* heuristic : {True, False}
  * If true, will use a heuristic MIP warm start obtained from [GerryChain](https://gerrychain.readthedocs.io/en/latest/)
* lp : {True, False} 
  * If true, will create a (separate) model for the LP relaxation and solve it to evaluate LP strength

