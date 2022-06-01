# SpecPart

SpecPart: A Supervised Spectral Framework for Hypergraph Partitioning Solution Improvement

*******************************************************************************************

## Description ##

This repository serves the following purposes:

1. We provide the ISPD98 benchmarks (with unit vertex weights and actual vertex weights) and Titan23 benchmarks in hMETIS format. To understand this format please refer to the hMETIS manual. 
2. We provde source code for the Julia implementation of SpecPart. We also provide the implementation of CMG (Combinatorial Multigrid) preconditioner with this package. 
3. We provide the best partitioning solutions for the ISPD98 benchmarks (both with unit vertex weights and actual vertex weights) and Titan23 benchmarks. 
4. We provide a "Golden Evaluator" which processes a partition file to report the cutsize and the block balances. 
5. We provide a leaderboard of best cuts found on the ISPD98 benchmarks and the Titan23 benchmarks. We encourage fellow researchers to update the leaderboard if better solutions are found. 

We acknowledge that further improvement on existing solutions is possible and we will continue to update the leaderboard and maintain the repository as we keep on doing so. 

*******************************************************************************************
## Current File/Directory Tree ##

    .
    |── Leaderboard
    |   └── ISPD98_Leaderboard                # Leaderboard of cuts for ISPD98 benchmarks (unit vertex weights and actual vertex weights)
    |   └── Titan23_Leaderboard               # Leaderboard of cuts for Titan23 benchmarks
    |
    ├── SpecPart                              # SpecPart Impelementation
    │   ├── SpectralCommunityDetection
    │   │   ├── cmg                           # Combinatorial Multigrid Implementation
    |
    ├── benchmark                             # hypergraph files for each benchmark
    │   ├── ISPD_benchmark                    # ISPD98 VLSI Circuit Benchmark Suite
    │   ├── ISPD_weight_benchmark             # ISPD98 VLSI Circuit Benchmark Suite with vertex weight
    │   └── Titan23_benchmark                 # Titan23 Benchmark Suite
    │   
    ├── golden_evaluator                      # script for evaluting partitioning solutions
    │ 
    ├── solutions                             # solutions for all the benchmarks under different scenarios
    │   ├── ISPD_benchmark_solutions
    │   |   ├── hMetis            
    │   |   |   ├── UBfactor_2
    │   |   |   └── UBfactor_10
    │   |   |   
    │   |   ├── KaHyPar         
    │   |   |   ├── UBfactor_2
    │   |   |   └── UBfactor_10
    │   |   |   
    |   |   └── SpecPart    
    │   |       ├── hMetis_SpecPart
    |   │       |   ├── UBfactor_2
    │   |       |   └── UBfactor_10
    │   |       |
    │   |       └── KaHyPar_SpecPart
    |   │           ├── UBfactor_2
    │   |           └── UBfactor_10
    │   |       
    │   ├── ISPD_weight_benchmark_solutions
    │   |   ├── hMetis            
    │   |   |   ├── UBfactor_2
    │   |   |   └── UBfactor_10
    │   |   |   
    │   |   ├── KaHyPar         
    │   |   |   ├── UBfactor_2
    │   |   |   └── UBfactor_10
    │   |   |   
    │   |   └── SpecPart    
    │   |       ├── hMetis_SpecPart
    |   │       |   ├── UBfactor_2
    │   |       |   └── UBfactor_10
    │   |       |
    │   |       └── KaHyPar_SpecPart
    |   │           ├── UBfactor_2
    │   |           └── UBfactor_10
    │   |       
    └── └── Titan23_benchmark_solutions
            ├── hMetis            
            |   ├── UBfactor_2
            |   └── UBfactor_20
            |   
            ├── hMetis_Autotune          
            |   └── UBfactor_10
            |   
            ├── hMetis_Autotune_SpecPart        
            |   └── UBfactor_10
            |   
            └── hMetis_SpecPart    
                ├── UBfactor_2
                └── UBfactor_20  
    
*******************************************************************************************
  
## Run Instructions ## 
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

Please note that current version of SpecPart supprorts bipartitions only. We will publish future versions of the code to tackle k-way partitions. 

*******************************************************************************************

## Leaderboards ##

Current Leaderboard of Cuts on ISPD98 testcases with unit vertex weights and actual weights: 

|   Testcase   | Statistics |              |    Cutsize    |              |              |               |
|:------------:|:----------:|:------------:|:-------------:|:------------:|:------------:|:-------------:|
|              | # Vertices | # Hyperedges | UBfactor  = 1 | UBfactor = 2 | UBfactor = 5 | UBfactor = 10 |
|     IBM01    |    12752   |     14111    |      203      |      200     |      180     |      166      |
| IBM01_weight |    12752   |     14111    |      216      |      215     |      215     |      215      |
|     IBM02    |    19601   |     19584    |      341      |      307     |      262     |      262      |
| IBM02_weight |    19601   |     19584    |      266      |      266     |      260     |      207      |
|     IBM03    |    23136   |     27401    |      954      |      951     |      950     |      950      |
| IBM03_weight |    23136   |     27401    |      775      |      691     |      680     |      580      |
|     IBM04    |    27507   |     31970    |      580      |      573     |      514     |      388      |
| IBM04_weight |    27507   |     31970    |      496      |      475     |      438     |      393      |
|     IBM05    |    29347   |     28446    |      1719     |     1706     |     1697     |      1645     |
| IBM05_weight |    29347   |     28446    |      1724     |     1722     |     1693     |      1647     |
|     IBM06    |    32498   |     34826    |      976      |      962     |      871     |      728      |
| IBM06_weight |    32498   |     34826    |      483      |      442     |      363     |      297      |
|     IBM07    |    45926   |     48117    |      910      |      878     |      818     |      764      |
| IBM07_weight |    45926   |     48117    |      737      |      736     |      717     |      636      |
|     IBM08    |    51309   |     50513    |      1140     |     1140     |     1140     |      1120     |
| IBM08_weight |    51309   |     50513    |      1170     |     1168     |     1120     |      972      |
|     IBM09    |    53395   |     60902    |      625      |      620     |      620     |      519      |
| IBM09_weight |    53395   |     60902    |      519      |      519     |      520     |      522      |
|     IBM10    |    69429   |     75196    |      1285     |     1253     |     1244     |      1244     |
| IBM10_weight |    69429   |     75196    |      1030     |      947     |      732     |      413      |
|     IBM11    |    70558   |     81454    |      1065     |     1051     |      951     |      763      |
| IBM11_weight |    70558   |     81454    |      767      |      766     |      687     |      656      |
|     IBM12    |    71076   |     77240    |      1934     |     1919     |     1871     |      1841     |
| IBM12_weight |    71076   |     77240    |      1976     |     1973     |     1976     |      1855     |
|     IBM13    |    84199   |     99666    |      837      |      831     |      831     |      655      |
| IBM13_weight |    84199   |     99666    |      892      |      846     |      858     |      793      |
|     IBM14    |   147605   |    152772    |      1842     |     1842     |     1794     |      1509     |
| IBM14_weight |   147605   |    152772    |      1740     |     1674     |     1501     |      1282     |
|     IBM15    |   161570   |    186608    |      2743     |     2730     |     2530     |      2135     |
| IBM15_weight |   161570   |    186608    |      2014     |     1913     |     1771     |      1444     |
|     IBM16    |   183484   |    190048    |      1975     |     1827     |     1695     |      1619     |
| IBM16_weight |   183484   |    190048    |      1656     |     1634     |     1638     |      1649     |
|     IBM17    |   185495   |    189581    |      2314     |     2270     |     2186     |      1989     |
| IBM17_weight |   185495   |    189581    |      2302     |     2289     |     2187     |      2018     |
|     IBM18    |   210613   |    201920    |      1521     |     1521     |     1521     |      1520     |
| IBM18_weight |   210613   |    201920    |      1641     |     1588     |     1520     |      1520     |

*******************************************************************************************

Current Leaderboard of Cuts on Titan23 testcases: 

|               | Statistics |             |    Cutsize   |              |              |               |               |
|---------------|:----------:|:-----------:|:------------:|:------------:|:------------:|:-------------:|:-------------:|
|    Testcase   |  #Vertices | #Hyperedges | UBfactor = 1 | UBfactor = 2 | UBfactor = 5 | UBfactor = 10 | UBfactor = 20 |
|  sparcT1_core |    91976   |    92827    |     1088     |     1012     |     1001     |      1090     |      903      |
|     neuron    |    92290   |    125305   |      239     |      252     |      252     |      245      |      206      |
|  stereovision |    94050   |    127085   |      186     |      180     |      182     |       91      |       91      |
|     des90     |   111221   |    139557   |      416     |      402     |      393     |      407      |      358      |
|  SLAM_spheric |   113115   |    142408   |     1061     |     1061     |     1061     |      1061     |      1061     |
|  cholesky_mc  |   113250   |    144948   |      343     |      285     |      286     |      353      |      345      |
|  segmentation |   138295   |    179051   |      165     |      126     |      125     |      126      |       78      |
|  bitonic_mesh |   192064   |    235328   |      588     |      587     |      581     |      543      |      483      |
|      dart     |   202354   |    223301   |      807     |      807     |      815     |      719      |      540      |
|     openCV    |   217453   |    284108   |      566     |      510     |      506     |      525      |      518      |
|    stap_qrd   |   240240   |    290123   |      398     |      399     |      380     |      296      |      295      |
|     minres    |   261359   |    320540   |      241     |      215     |      215     |      199      |      189      |
| cholesky_bdti |   266422   |    342688   |     1213     |     1156     |     1156     |      1156     |      947      |
|    denoise    |   275638   |    356848   |      457     |      416     |      416     |      416      |      224      |
|  sparcT2_core |   300109   |    302663   |     1227     |     1244     |     1239     |      1252     |      1245     |
|   gsm_switch  |   493260   |    507821   |     1836     |     1827     |     1615     |      1460     |      1407     |
|    mes_noc    |   547544   |    577664   |      703     |      634     |      717     |      634      |      617      |
|     LU230     |   574372   |    669477   |     3362     |     3273     |     3339     |      3262     |      2677     |
|   LU_Network  |   635456   |    726999   |      646     |      525     |      524     |      524      |      524      |
| sparcT1_chip2 |   820886   |    821274   |     1037     |      899     |     1137     |      815      |      783      |
|    directrf   |   931275   |   1374742   |      673     |      574     |      587     |      378      |      295      |
| bitcoin_miner |   1089284  |   1448151   |     1512     |     1297     |     1232     |      1232     |      1225     |

*******************************************************************************************

## Authors ##
Ismail Bustany, Andrew Kahng, Yiannis Koutis, Bodhisatta Pramanik, Zhiang Wang
