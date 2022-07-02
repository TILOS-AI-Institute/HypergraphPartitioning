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
```              
