# Hypergraph Partitioning for VLSI: Benchmarks, Code and Leaderboard

Hypergraph Partitioning Leaderboard: Leaderboard of minimum hyperedge cut values for different testcases with multiple imbalance factors.

SpecPart: A Supervised Spectral Framework for Hypergraph Partitioning Solution Improvement

*******************************************************************************************

## Description ##

This repository supports "Data, Benchmarking and Roadmapping" goals for the balanced hypergraph min-cut partitioning problem, which is central to chip design and the divide-and-conquer paradigm. The repository is a one-stop shop" that serves the following purposes.

1. We provide the [ISPD98 benchmarks](https://dl.acm.org/doi/10.1145/274535.274546) (with unit vertex weights and actual vertex weights) and [Titan23 benchmarks](https://www.eecg.utoronto.ca/~kmurray/titan.html) in [hMETIS](http://glaros.dtc.umn.edu/gkhome/metis/hmetis/overview) format. To understand this format please refer to the [hMETIS manual](http://glaros.dtc.umn.edu/gkhome/metis/hmetis/download). 
2. We provde open source code for the Julia implementation of a recent hypergraph partitioning code, [SpecPart](https://github.com/TILOS-AI-Institute/HypergraphPartitioning/blob/main/SpecPart/SpecPart_final_submission.pdf). We also provide the implementation of CMG (Combinatorial Multigrid) preconditioner with this package. 
3. We provide the best known partitioning solutions for the ISPD98 benchmarks (both with unit vertex weights and actual vertex weights) and the Titan Benchmarks. 
4. We provide a "Golden Evaluator" which processes a partition file to report the cutsize and the block balances. 
5. We provide a leaderboard of minimum hyperedge cutsize values found on the ISPD98 benchmarks and the Titan23 benchmarks. We encourage fellow researchers to update the leaderboard if better solutions are found. 

We hope to see pull requests with new solutions and optimization codes, and will continue to update the repository and leaderboard as we this happens.

## Table of Contents ##

1. Current file/directory tree and description
2. (SpecPart-specific) installation and run instructions
3. Leaderboard of minimum hyperedge cutsize values
4. Authors

*******************************************************************************************
## Current File/Directory Tree and Description ##

    .
    |── Leaderboard
    |   └── ISPD98_Leaderboard                # Leaderboard of minimum hyperedge cut values for ISPD98 benchmarks (unit vertex weights and actual vertex weights)
    |   └── Titan23_Leaderboard               # Leaderboard of minimum hyperedge cut values for Titan23 benchmarks
    |
    ├── SpecPart                              # SpecPart Implementation
    │   ├── SpectralCommunityDetection
    │   │   ├── cmg                           # Combinatorial Multigrid Implementation
    |
    ├── benchmark                             # hypergraph files for each benchmark
    │   ├── ISPD_benchmark                    # ISPD98 VLSI Circuit Benchmark Suite with unit vertex weights
    │   ├── ISPD_weight_benchmark             # ISPD98 VLSI Circuit Benchmark Suite with actual vertex weights
    │   └── Titan23_benchmark                 # Titan23 Benchmark Suite
    │   
    ├── golden_evaluator                      # script for evaluating partitioning solutions
    │ 
    ├── solutions                             # solutions for all the benchmarks with different imbalance factors
    │   ├── ISPD_benchmark_solutions          # solutions on ISPD98 testcases with unit vertex weights
    │   |   ├── hMetis                        # solutions using hMETIS
    │   |   |   ├── UBfactor_2                # solutions with imbalance factor 2
    │   |   |   └── UBfactor_10               # solutions with imbalance factor 10
    │   |   |   
    │   |   ├── KaHyPar                       # solutions with KaHyPar
    │   |   |   ├── UBfactor_2                # solutions with imbalance factor 2
    │   |   |   └── UBfactor_10               # solutions with imbalance factor 10
    │   |   |   
    |   |   └── SpecPart                      # solutions with SpecPart
    │   |       ├── hMetis_SpecPart           # SpecPart with initial solution generated from hMETIS
    |   │       |   ├── UBfactor_2            # solutions with imbalance factor 2
    │   |       |   └── UBfactor_10           # solutions with imbalance factor 10
    │   |       |
    │   |       └── KaHyPar_SpecPart          # SpecPart with initial solution generated from KaHyPar
    |   │           ├── UBfactor_2            # solutions with imbalance factor 2
    │   |           └── UBfactor_10           # solutions with imbalance factor 10
    │   |       
    │   ├── ISPD_weight_benchmark_solutions   # solutions on ISPD98 testcases with actual vertex weights
    │   |   ├── hMetis                        # solutions with hMETIS
    │   |   |   ├── UBfactor_2                # solutions with imbalance factor 2
    │   |   |   └── UBfactor_10               # solutions with imbalance factor 10
    │   |   |   
    │   |   ├── KaHyPar                       # solutions with KaHyPar
    │   |   |   ├── UBfactor_2                # solutions with imbalance factor 2
    │   |   |   └── UBfactor_10               # solutions with imbalance factor 10
    │   |   |   
    │   |   └── SpecPart                      # solutions with SpecPart
    │   |       ├── hMetis_SpecPart           # SpecPart with initial solution generated from hMETIS
    |   │       |   ├── UBfactor_2            # solutions with imbalance factor 2
    │   |       |   └── UBfactor_10           # solutions with imbalance factor 10
    │   |       |   
    │   |       └── KaHyPar_SpecPart          # SpecPart with initial solution generated from KaHyPar
    |   │           ├── UBfactor_2            # solutions with imbalance factor 2
    │   |           └── UBfactor_10           # solutions with imbalance factor 10
    │   |       
    └── └── Titan23_benchmark_solutions       # solutions on Titan23 testcases
            ├── hMetis                        # solutions with hMETIS
            |   ├── UBfactor_2                # solutions with imbalance factor 2
            |   └── UBfactor_20               # solutions with imbalance factor 20
            |   
            ├── hMetis_Autotune               # solutions with Autotuned hMETIS
            |   └── UBfactor_10               # solutions with imbalance factor 10
            |   
            ├── hMetis_Autotune_SpecPart      # SpecPart with initial solution generated from Autotuned hMETIS   
            |   └── UBfactor_10               # solutions with imbalance factor 10
            |   
            └── hMetis_SpecPart               # SpecPart with initial solution generated from hMETIS   
                ├── UBfactor_2                # solutions with imbalance factor 2
                └── UBfactor_20               # solutions with imbalance factor 20
    
*******************************************************************************************
  
## Installation and Run Instructions ## 

Open Julia REPL and do the following: 
``` 
include("SpectralRefinement.jl")
using Main.SpecPart
SpecPart.SpectralRefinement(hg = "Hypergraph file", pfile = "Partition file", Nparts = "Number of partitions", cycles = ζ, hyperedges_threshold = γ, ub = ε, nev = m, refine_iters = β, best_solns = δ)

ζ is the number of graph cycles (we recommend ζ=2)
γ is the threshold number of hyperedges to run ILP (we recommend γ=600)
ε is the imbalance factor (ε=1-49)
m is the number of eigenvectors (we recommend m=2)
β is the number of SpecPart iterations (we recommend β=2)
δ is the number of solutions picked from candidate solutions for overlay based clustering (we recommend δ=5)

```

Please note that current version of SpecPart supports bipartitions only. We will publish future versions of the code to tackle k-way partitions. A user guide is available [here](https://github.com/TILOS-AI-Institute/HypergraphPartitioning/tree/main/SpecPart/README.md). 

*******************************************************************************************

## Leaderboards of minimum hyperedge cut values ##

Current Leaderboard of minimum hyperedge cut values on [ISPD98 testcases](https://dl.acm.org/doi/10.1145/274535.274546) with unit vertex weights and actual vertex weights, with different imbalance factors (ε):  

|   Testcase   | Statistics |              |    Cutsize    |              |              |               |
|:------------:|:----------:|:------------:|:-------------:|:------------:|:------------:|:-------------:|
|              | # Vertices | # Hyperedges | ε  = 1 | ε = 2 | ε = 5 | ε = 10 |
|     IBM01    |    12752   |     14111    |      203      |      200     |      180     |      166      |
| IBM01_wt |    12752   |     14111    |      216      |      215     |      215     |      215      |
|     IBM02    |    19601   |     19584    |      341      |      307     |      262     |      262      |
| IBM02_wt |    19601   |     19584    |      266      |      266     |      260     |      207      |
|     IBM03    |    23136   |     27401    |      954      |      951     |      950     |      950      |
| IBM03_wt |    23136   |     27401    |      775      |      691     |      680     |      580      |
|     IBM04    |    27507   |     31970    |      580      |      573     |      514     |      388      |
| IBM04_wt |    27507   |     31970    |      496      |      475     |      438     |      393      |
|     IBM05    |    29347   |     28446    |      1719     |     1706     |     1697     |      1645     |
| IBM05_wt |    29347   |     28446    |      1724     |     1722     |     1693     |      1647     |
|     IBM06    |    32498   |     34826    |      976      |      962     |      871     |      728      |
| IBM06_wt |    32498   |     34826    |      483      |      442     |      363     |      297      |
|     IBM07    |    45926   |     48117    |      910      |      878     |      818     |      764      |
| IBM07_wt |    45926   |     48117    |      737      |      736     |      717     |      636      |
|     IBM08    |    51309   |     50513    |      1140     |     1140     |     1140     |      1120     |
| IBM08_wt |    51309   |     50513    |      1170     |     1168     |     1120     |      972      |
|     IBM09    |    53395   |     60902    |      625      |      620     |      620     |      519      |
| IBM09_wt |    53395   |     60902    |      519      |      519     |      519     |      522      |
|     IBM10    |    69429   |     75196    |      1285     |     1253     |     1244     |      1244     |
| IBM10_wt |    69429   |     75196    |      1030     |      947     |      732     |      413      |
|     IBM11    |    70558   |     81454    |      1065     |     1051     |      951     |      763      |
| IBM11_wt |    70558   |     81454    |      767      |      766     |      687     |      656      |
|     IBM12    |    71076   |     77240    |      1934     |     1919     |     1871     |      1841     |
| IBM12_wt |    71076   |     77240    |      1976     |     1973     |     1976     |      1855     |
|     IBM13    |    84199   |     99666    |      837      |      831     |      831     |      655      |
| IBM13_wt |    84199   |     99666    |      892      |      846     |      846     |      793      |
|     IBM14    |   147605   |    152772    |      1842     |     1842     |     1794     |      1509     |
| IBM14_wt |   147605   |    152772    |      1740     |     1674     |     1501     |      1282     |
|     IBM15    |   161570   |    186608    |      2743     |     2730     |     2530     |      2135     |
| IBM15_wt |   161570   |    186608    |      2014     |     1913     |     1771     |      1444     |
|     IBM16    |   183484   |    190048    |      1975     |     1827     |     1695     |      1619     |
| IBM16_wt |   183484   |    190048    |      1656     |     1634     |     1634     |      1634     |
|     IBM17    |   185495   |    189581    |      2314     |     2270     |     2186     |      1989     |
| IBM17_wt |   185495   |    189581    |      2302     |     2289     |     2187     |      2018     |
|     IBM18    |   210613   |    201920    |      1521     |     1521     |     1521     |      1520     |
| IBM18_wt |   210613   |    201920    |      1641     |     1588     |     1520     |      1520     |

*******************************************************************************************

Current Leaderboard of minimum hyperedge cut values on [Titan23 testcases](https://www.eecg.utoronto.ca/~kmurray/titan.html) with different imbalance factors (ε): 

|               | Statistics |             |    Cutsize   |              |              |               |               |
|---------------|:----------:|:-----------:|:------------:|:------------:|:------------:|:-------------:|:-------------:|
|    Testcase   |  #Vertices | #Hyperedges | ε = 1 | ε = 2 | ε = 5 | ε = 10 | ε = 20 |
|  sparcT1_core |    91976   |    92827    |     1088     |     1012     |     1001     |      1090     |      903      |
|     neuron    |    92290   |    125305   |      239     |      239     |      239     |      239      |      206      |
|  stereovision |    94050   |    127085   |      186     |      180     |      180     |       91      |       91      |
|     des90     |   111221   |    139557   |      416     |      402     |      393     |      393      |      358      |
|  SLAM_spheric |   113115   |    142408   |     1061     |     1061     |     1061     |      1061     |      1061     |
|  cholesky_mc  |   113250   |    144948   |      343     |      285     |      281     |      281      |      285      |
|  segmentation |   138295   |    179051   |      165     |      126     |      107     |      81      |       78      |
|  bitonic_mesh |   192064   |    235328   |      588     |      587     |      581     |      543      |      483      |
|      dart     |   202354   |    223301   |      807     |      807     |      807     |      719      |      540      |
|     openCV    |   217453   |    284108   |      566     |      510     |      506     |      506      |      506      |
|    stap_qrd   |   240240   |    290123   |      398     |      398     |      380     |      296      |      295      |
|     minres    |   261359   |    320540   |      241     |      215     |      215     |      199      |      189      |
| cholesky_bdti |   266422   |    342688   |     1213     |     1156     |     1156     |      1156     |      947      |
|    denoise    |   275638   |    356848   |      457     |      416     |      416     |      416      |      224      |
|  sparcT2_core |   300109   |    302663   |     1227     |     1227     |     1227     |      1227     |      1227     |
|   gsm_switch  |   493260   |    507821   |     1836     |     1827     |     1615     |      1460     |      1407     |
|    mes_noc    |   547544   |    577664   |      703     |      634     |      634     |      634      |      617      |
|     LU230     |   574372   |    669477   |     3362     |     3273     |     3339     |      3262     |      2677     |
|   LU_Network  |   635456   |    726999   |      646     |      525     |      524     |      524      |      524      |
| sparcT1_chip2 |   820886   |    821274   |     1037     |      899     |      899     |      815      |      783      |
|    directrf   |   931275   |   1374742   |      673     |      574     |      574     |      378      |      295      |
| bitcoin_miner |   1089284  |   1448151   |     1512     |     1297     |     1232     |      1232     |      1225     |

*******************************************************************************************

## Authors ##
Ismail Bustany, Andrew Kahng, Yiannis Koutis, Bodhisatta Pramanik, Zhiang Wang




## References ##
To reference SpecPart, please cite:

I. Bustany, A. B. Kahng, Y. Koutis, B. Pramanik and Z. Wang, "SpecPart: A Supervised Spectral Framework for Hypergraph Partitioning Solution Improvement", ([pdf](https://github.com/TILOS-AI-Institute/HypergraphPartitioning/blob/main/SpecPart/SpecPart_final_submission.pdf)),
*Proc. ACM/IEEE International Conference on Computer-Aided Design*, 2022, to appear.
