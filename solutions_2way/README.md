# Hypergraph partition solutions #

This directory contains the hypergraph solution files corresponding to the best published hyperedge cutsizes for different imbalance factors.

Following is a map of the directories in this folder:

```

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
    │   |   ├── MtKaHyPar                     # solutions with MtKaHyPar
    │   |   |   ├── UBfactor_2                # solutions with imbalance factor 2
    │   |   |   └── UBfactor_10               # solutions with imbalance factor 10
    │   |   |
    │   |   ├── hMetis_SpecPart               # solutions with hMetis_SpecPart
    │   |   |   ├── UBfactor_2                # solutions with imbalance factor 2
    │   |   |   └── UBfactor_10               # solutions with imbalance factor 10
    │   |   |
    │   |   └── KaHyPar_SpecPart              # solutions with KaHyPar_SpecPart
    │   |   |   ├── UBfactor_2                # solutions with imbalance factor 2
    │   |   |   └── UBfactor_10               # solutions with imbalance factor 10
    |   |   |
    │   |   └── HyperEF_2.0                   # solutions with HyperEF 2.0
    │   |       ├── UBfactor_2                # solutions with imbalance factor 2
    │   |       └── UBfactor_10               # solutions with imbalance factor 10
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
    │   |   ├── MtKaHyPar                     # solutions with MtKaHyPar
    │   |   |   ├── UBfactor_2                # solutions with imbalance factor 2
    │   |   |   └── UBfactor_10               # solutions with imbalance factor 10
    │   |   |
    │   |   ├── hMetis_SpecPart               # solutions with hMetis_SpecPart
    │   |   |   ├── UBfactor_2                # solutions with imbalance factor 2
    │   |   |   └── UBfactor_10               # solutions with imbalance factor 10
    │   |   |
    │   |   └── KaHyPar_SpecPart              # solutions with KaHyPar_SpecPart
    │   |       ├── UBfactor_2                # solutions with imbalance factor 2
    │   |       └── UBfactor_10               # solutions with imbalance factor 10
    │   |
    └── └── Titan23_benchmark_solutions       # solutions on Titan23 testcases
            ├── hMetis                        # solutions with hMETIS
            |   ├── UBfactor_2                # solutions with imbalance factor 2
            |   └── UBfactor_20               # solutions with imbalance factor 20
            |
            ├── MtKaHyPar                     # solutions with Autotuned MtKaHyPar
            |   ├── UBfactor_2                # solutions with imbalance factor 2
            |   └── UBfactor_10               # solutions with imbalance factor 10
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
```

## Leaderboards of minimum hyperedge cut values ##

Current Leaderboard of minimum hyperedge cut values on [ISPD98 testcases](https://dl.acm.org/doi/10.1145/274535.274546) with unit vertex weights, with different imbalance factors (ε):

| Instance | Statistics | | Cutsize | |
|:------------:|:------------:|:------------:|:------------:|:------------:|
|  | # Vertices | # Hyperedges | ε = 2% | ε = 10% |
| ibm01 | 12752 | 14111 | 202<br><sub>KaHyPar</sub><br><sub>KaHyPar_SpecPart</sub> | 166<br><sub>KaHyPar</sub> |
| ibm02 | 19601 | 19584 | 336<br><sub>KaHyPar</sub><br><sub>hMetis_SpecPart</sub> | 262<br><sub>hMetis</sub><br><sub>KaHyPar</sub><br><sub>hMetis_SpecPart</sub><br><sub>MtKaHyPar</sub> |
| ibm03 | 23136 | 27401 | 954<br><sub>KaHyPar</sub> | 950<br><sub>MtKaHyPar</sub> |
| ibm05 | 29347 | 28446 | 1719<br><sub>KaHyPar</sub> | 1645<br><sub>KaHyPar</sub> |
| ibm04 | 27507 | 31970 | 579<br><sub>KaHyPar</sub> | 388<br><sub>hMetis</sub><br><sub>KaHyPar</sub><br><sub>hMetis_SpecPart</sub><br><sub>KaHyPar_SpecPart</sub> |
| ibm06 | 32498 | 34826 | 966<br><sub>KaHyPar</sub> | 733<br><sub>KaHyPar_SpecPart</sub> |
| ibm07 | 45926 | 48117 | 904<br><sub>KaHyPar</sub> | 760<br><sub>KaHyPar</sub><br><sub>KaHyPar_SpecPart</sub> |
| ibm08 | 51309 | 50513 | 1140<br><sub>KaHyPar</sub><br><sub>MtKaHyPar</sub> | 1120<br><sub>KaHyPar</sub> |
| ibm09 | 53395 | 60902 | 620<br><sub>KaHyPar</sub><br><sub>KaHyPar_SpecPart</sub> | 519<br><sub>KaHyPar</sub><br><sub>KaHyPar_SpecPart</sub> |
| ibm10 | 69429 | 75196 | 1313<br><sub>hMetis</sub><br><sub>KaHyPar</sub> | 1250<br><sub>KaHyPar</sub> |
| ibm12 | 71076 | 77240 | 1920<br><sub>KaHyPar</sub><br><sub>KaHyPar_SpecPart</sub> | 1842<br><sub>KaHyPar</sub><br><sub>KaHyPar_SpecPart</sub> |
| ibm11 | 70558 | 81454 | 1062<br><sub>KaHyPar</sub><br><sub>hMetis_SpecPart</sub> | 764<br><sub>hMetis_SpecPart</sub> |
| ibm13 | 84199 | 99666 | 831<br><sub>KaHyPar</sub> | 671<br><sub>KaHyPar</sub> |
| ibm14 | 147605 | 152772 | 1850<br><sub>MtKaHyPar</sub> | 1518<br><sub>KaHyPar</sub> |
| ibm15 | 161570 | 186608 | 2728<br><sub>KaHyPar</sub> | 2135<br><sub>KaHyPar</sub> |
| ibm17 | 185495 | 189581 | 2310<br><sub>KaHyPar</sub> | 1989<br><sub>KaHyPar</sub><br><sub>KaHyPar_SpecPart</sub> |
| ibm16 | 183484 | 190048 | 1882<br><sub>KaHyPar</sub> | 1619<br><sub>KaHyPar</sub><br><sub>KaHyPar_SpecPart</sub> |
| ibm18 | 210613 | 201920 | 1539<br><sub>hMetis_SpecPart</sub> | 1550<br><sub>hMetis</sub><br><sub>hMetis_SpecPart</sub> |

Current Leaderboard of minimum hyperedge cut values on [ISPD98 testcases](https://dl.acm.org/doi/10.1145/274535.274546) with actual vertex weights, with different imbalance factors (ε):

| Instance | Statistics | | Cutsize | |
|------------|:------------:|:------------:|:------------:|:------------:|
|  | # Vertices | # Hyperedges | ε = 2% | ε = 10% |
| ibm01 | 12752 | 14111 | 216<br><sub>MtKaHyPar</sub> | 215<br><sub>MtKaHyPar</sub> |
| ibm02 | 19601 | 19584 | 267<br><sub>MtKaHyPar</sub> | 258<br><sub>MtKaHyPar</sub> |
| ibm03 | 23136 | 27401 | 724<br><sub>MtKaHyPar</sub> | 642<br><sub>MtKaHyPar</sub> |
| ibm05 | 29347 | 28446 | 1719<br><sub>KaHyPar</sub> | 1645<br><sub>hMetis</sub><br><sub>KaHyPar</sub> |
| ibm04 | 27507 | 31970 | 502<br><sub>MtKaHyPar</sub> | 442<br><sub>MtKaHyPar</sub> |
| ibm06 | 32498 | 34826 | 570<br><sub>MtKaHyPar</sub> | 377<br><sub>MtKaHyPar</sub> |
| ibm07 | 45926 | 48117 | 807<br><sub>MtKaHyPar</sub> | 695<br><sub>hMetis_SpecPart</sub> |
| ibm08 | 51309 | 50513 | 1245<br><sub>MtKaHyPar</sub> | 1120<br><sub>hMetis</sub><br><sub>KaHyPar</sub> |
| ibm09 | 53395 | 60902 | 519<br><sub>MtKaHyPar</sub> | 519<br><sub>hMetis</sub><br><sub>KaHyPar</sub><br><sub>KaHyPar_SpecPart</sub><br><sub>MtKaHyPar</sub> |
| ibm10 | 69429 | 75196 | 1191<br><sub>MtKaHyPar</sub> | 817<br><sub>MtKaHyPar</sub> |
| ibm12 | 71076 | 77240 | 1969<br><sub>MtKaHyPar</sub> | 1961<br><sub>MtKaHyPar</sub> |
| ibm11 | 70558 | 81454 | 765<br><sub>KaHyPar_SpecPart</sub> | 650<br><sub>hMetis_SpecPart</sub> |
| ibm13 | 84199 | 99666 | 949<br><sub>MtKaHyPar</sub> | 831<br><sub>MtKaHyPar</sub> |
| ibm14 | 147605 | 152772 | 1767<br><sub>MtKaHyPar</sub> | 1492<br><sub>MtKaHyPar</sub> |
| ibm15 | 161570 | 186608 | 2430<br><sub>MtKaHyPar</sub> | 1790<br><sub>MtKaHyPar</sub> |
| ibm17 | 185495 | 189581 | 2272<br><sub>MtKaHyPar</sub> | 2194<br><sub>MtKaHyPar</sub> |
| ibm16 | 183484 | 190048 | 1624<br><sub>MtKaHyPar</sub> | 1623<br><sub>MtKaHyPar</sub> |
| ibm18 | 210613 | 201920 | 1690<br><sub>MtKaHyPar</sub> | 1532<br><sub>hMetis_SpecPart</sub> |

Current Leaderboard of minimum hyperedge cut values on [Titan23 testcases](https://www.eecg.utoronto.ca/~kmurray/titan.html) with different imbalance factors (ε):

| Instance | Statistics | | Cutsize | |
|------------|:------------:|:------------:|:------------:|:------------:|
|  | # Vertices | # Hyperedges | ε = 2% | ε = 20% |
| sparcT1_core | 91976 | 92827 | 976<br><sub>MtKaHyPar</sub> | 903<br><sub>hMetis_SpecPart</sub> |
| neuron | 92290 | 125305 | 252<br><sub>hMetis_SpecPart</sub> | 206<br><sub>hMetis_SpecPart</sub> |
| stereo_vision | 94050 | 127085 | 169<br><sub>KaHyPar</sub> | 91<br><sub>KaHyPar</sub><br><sub>hMetis_SpecPart</sub><br><sub>MtKaHyPar</sub> |
| des90 | 111221 | 139557 | 402<br><sub>hMetis_SpecPart</sub> | 358<br><sub>hMetis_SpecPart</sub> |
| SLAM_spheric | 113115 | 142408 | 1061<br><sub>hMetis</sub><br><sub>KaHyPar</sub><br><sub>hMetis_SpecPart</sub><br><sub>MtKaHyPar</sub> | 1061<br><sub>hMetis</sub><br><sub>KaHyPar</sub><br><sub>hMetis_SpecPart</sub><br><sub>MtKaHyPar</sub> |
| cholesky_mc | 113250 | 144948 | 285<br><sub>hMetis_SpecPart</sub> | 281<br><sub>KaHyPar</sub><br><sub>MtKaHyPar</sub> |
| segmentation | 138295 | 179051 | 126<br><sub>hMetis_SpecPart</sub> | 78<br><sub>hMetis_SpecPart</sub> |
| dart | 202354 | 223301 | 784<br><sub>MtKaHyPar</sub> | 543<br><sub>hMetis_SpecPart</sub> |
| bitonic_mesh | 192064 | 235328 | 585<br><sub>hMetis_SpecPart</sub> | 483<br><sub>hMetis_SpecPart</sub> |
| openCV | 217453 | 284108 | 511<br><sub>hMetis_SpecPart</sub> | 518<br><sub>hMetis_SpecPart</sub> |
| stap_qrd | 240240 | 290123 | 399<br><sub>hMetis</sub> | 275<br><sub>MtKaHyPar</sub> |
| sparcT2_core | 300109 | 302663 | 1185<br><sub>MtKaHyPar</sub> | 1183<br><sub>KaHyPar</sub><br><sub>MtKaHyPar</sub> |
| minres | 261359 | 320540 | 215<br><sub>hMetis</sub><br><sub>hMetis_SpecPart</sub> | 189<br><sub>hMetis</sub><br><sub>hMetis_SpecPart</sub> |
| cholesky_bdti | 266422 | 342688 | 1156<br><sub>hMetis_SpecPart</sub> | 1024<br><sub>hMetis</sub><br><sub>hMetis_SpecPart</sub> |
| denoise | 275638 | 356848 | 416<br><sub>hMetis_SpecPart</sub> | 224<br><sub>hMetis_SpecPart</sub> |
| gsm_switch | 493260 | 507821 | 1822<br><sub>KaHyPar</sub> | 1407<br><sub>hMetis_SpecPart</sub> |
| mes_noc | 547544 | 577664 | 634<br><sub>hMetis_SpecPart</sub> | 617<br><sub>hMetis_SpecPart</sub> |
| LU230 | 574372 | 669477 | 3273<br><sub>hMetis_SpecPart</sub> | 2677<br><sub>hMetis_SpecPart</sub> |
| LU_Network | 635456 | 726999 | 525<br><sub>hMetis_SpecPart</sub> | 524<br><sub>KaHyPar</sub><br><sub>hMetis_SpecPart</sub><br><sub>MtKaHyPar</sub> |
| sparcT1_chip2 | 820886 | 821274 | 899<br><sub>hMetis_SpecPart</sub> | 783<br><sub>hMetis_SpecPart</sub> |
| directrf | 931275 | 1374742 | 574<br><sub>hMetis_SpecPart</sub> | 295<br><sub>hMetis_SpecPart</sub> |
| bitcoin_miner | 1089284 | 1448151 | 1520<br><sub>KaHyPar</sub> | 1225<br><sub>hMetis_SpecPart</sub> |
