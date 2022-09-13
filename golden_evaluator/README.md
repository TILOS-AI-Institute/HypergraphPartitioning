# Golden Solution Evaluator #

We provide a Python script for evaluating partitioning solutions. The script serves the following purposes:
1. Evaluate the hyperedge cut size.
2. Validate the partitioning solution by evaluating whether the balance constraints are satisfied.

# Example use #

In order to run the golden evaluator on hypergraph solution files, follow the instructions:

```
Evaluator("benchmark file", "partition solution file", number_of_partitions, imbalance_factor)
```

We show how we run the golden evaluator on the "ibm02.hgr" benchmark with a solution file and with number of partitions 2 and imbalance factor 2.

```
>>> from golden_evaluator import Evaluator
>>> Evaluator("../benchmark/IBM/ibm02.hgr", "../solutions/ISPD_benchmark_solutions/hMetis/UBfactor_2/ibm02.hgr.k.2.UBfactor.2.seed.0", 2, 2)
Partition satisfies balance constraints!
cutsize   :    334
blocks_balance :   [0.5199224529360746, 0.4800775470639253]
```

# Leaderboard Generator #

We also provide a Python script to generate a leaderboard in markdown format. The script takes as input a folder with
solution files and the folder containing the actual instances. The output of the script can be used to update the leaderboard
in README.md.

# Example use #

```
./markdown_leaderboard_generator.py ../solutions/Titan23_benchmark_solutions ../benchmark/Titan23_benchmark
```

Output:
```
| Instance | Statistics | | Cutsize | |
|:------------:|:------------:|:------------:|:------------:|:------------:|
|  | # Vertices | # Hyperedges | ε = 2% | ε = 20% |
| sparcT1_core | 91976 | 92827 | 976<br><sub>MtKaHyPar</sub> | 903<br><sub>hMetis_SpecPart</sub> |
| neuron | 92290 | 125305 | 252<br><sub>hMetis_SpecPart</sub> | 206<br><sub>hMetis_SpecPart</sub> |
| stereo_vision | 94050 | 127085 | 170<br><sub>MtKaHyPar</sub> | 91<br><sub>hMetis_SpecPart</sub><br><sub>MtKaHyPar</sub> |
| des90 | 111221 | 139557 | 402<br><sub>hMetis_SpecPart</sub> | 358<br><sub>hMetis_SpecPart</sub> |
| SLAM_spheric | 113115 | 142408 | 1061<br><sub>hMetis</sub><br><sub>hMetis_SpecPart</sub><br><sub>MtKaHyPar</sub> | 1061<br><sub>hMetis</sub><br><sub>hMetis_SpecPart</sub><br><sub>MtKaHyPar</sub> |
| cholesky_mc | 113250 | 144948 | 285<br><sub>hMetis_SpecPart</sub> | 281<br><sub>MtKaHyPar</sub> |
| segmentation | 138295 | 179051 | 126<br><sub>hMetis_SpecPart</sub> | 78<br><sub>hMetis_SpecPart</sub> |
| dart | 202354 | 223301 | 784<br><sub>MtKaHyPar</sub> | 543<br><sub>hMetis_SpecPart</sub> |
| bitonic_mesh | 192064 | 235328 | 585<br><sub>hMetis_SpecPart</sub> | 483<br><sub>hMetis_SpecPart</sub> |
| openCV | 217453 | 284108 | 511<br><sub>hMetis_SpecPart</sub> | 518<br><sub>hMetis_SpecPart</sub> |
| stap_qrd | 240240 | 290123 | 399<br><sub>hMetis</sub> | 275<br><sub>MtKaHyPar</sub> |
| sparcT2_core | 300109 | 302663 | 1185<br><sub>MtKaHyPar</sub> | 1183<br><sub>MtKaHyPar</sub> |
| minres | 261359 | 320540 | 215<br><sub>hMetis</sub><br><sub>hMetis_SpecPart</sub> | 189<br><sub>hMetis</sub><br><sub>hMetis_SpecPart</sub> |
| cholesky_bdti | 266422 | 342688 | 1156<br><sub>hMetis_SpecPart</sub> | 1024<br><sub>hMetis</sub><br><sub>hMetis_SpecPart</sub> |
| denoise | 275638 | 356848 | 416<br><sub>hMetis_SpecPart</sub> | 224<br><sub>hMetis_SpecPart</sub> |
| gsm_switch | 493260 | 507821 | 1837<br><sub>hMetis_SpecPart</sub> | 1407<br><sub>hMetis_SpecPart</sub> |
| mes_noc | 547544 | 577664 | 634<br><sub>hMetis_SpecPart</sub> | 617<br><sub>hMetis_SpecPart</sub> |
| LU230 | 574372 | 669477 | 3273<br><sub>hMetis_SpecPart</sub> | 2677<br><sub>hMetis_SpecPart</sub> |
| LU_Network | 635456 | 726999 | 525<br><sub>hMetis_SpecPart</sub> | 524<br><sub>hMetis_SpecPart</sub><br><sub>MtKaHyPar</sub> |
| sparcT1_chip2 | 820886 | 821274 | 899<br><sub>hMetis_SpecPart</sub> | 783<br><sub>hMetis_SpecPart</sub> |
| directrf | 931275 | 1374742 | 574<br><sub>hMetis_SpecPart</sub> | 295<br><sub>hMetis_SpecPart</sub> |
| bitcoin_miner | 1089284 | 1448151 | 1566<br><sub>hMetis_SpecPart</sub> | 1225<br><sub>hMetis_SpecPart</sub> |
```

Note that when you run the script on a folder for the first time, it may take while to analyze all solutions. However, the script
caches the results and uses them when it runs again. If you change the solution files for an algorithm, you just have to delete the 'cached_results.json' file in the corresponding folder and run the script again.