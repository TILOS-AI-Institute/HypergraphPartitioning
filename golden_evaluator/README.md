# Golden Solution Evaluator #

We provide a Python script for evaluating partitioning solutions. The script serves the following purposes:
1. Evaluate the hyperedge cut size. 
2. Validate the partitioning solution by evaluating whether the balance constraints are satisfied. 

# Example use #

```
>>> from golden_evaluator import Evaluator
>>> Evaluator("../benchmark/IBM/ibm02.hgr", "solutions/ISPD_benchmark_solutions/hMetis/UBfactor_2/ibm02.hgr.k.2.UBfactor.2.seed.0", 2, 2)
Partition satisfies balance constraints!
cutsize   :    334
blocks_balance :   [0.5199224529360746, 0.4800775470639253]
```
