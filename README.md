# SpecPart: A Supervised Spectral Framework for Hypergraph Partitioning Solution Improvement

Results for SpecPart

*******************************************************************************************
In this repo, both hypergraph format and partitioning solution format
are the same as hMetis. Please refer to hMetis manual for detailed explaination
*******************************************************************************************

    .
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
    
  

