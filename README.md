# SpecPart: A Supervised Spectral Framework for Hypergraph Partitioning Solution Improvement

*******************************************************************************************

This repository serves the following purpuses:

1. We provide the ISPD98 benchmarks (with unit vertex weights and actual vertex weights) and Titan23 benchmarks in hMETIS format. To understand this format please refer to the hMETIS manual. 
2. We provde source code for the Julia implementation of SpecPart. We also provide the implementation of CMG (Combinatorial Multigrid) preconditioner with this package. 
3. We provide the best partitioning solutions for the ISPD98 benchmarks (both with unit vertex weights and actual vertex weights) and Titan23 benchmarks. 
4. We provide a "Golden Evaluator" which processes a partition file to report the cutsize and the block balances. 
5. We provide a leaderboard of best cuts found on the ISPD98 benchmarks and the Titan23 benchmarks. We encourage fellow researchers to update the leaderboard if better solutions are found. 

*******************************************************************************************

We acknowledge that further improvement on existing solutions is possible and we will continue to update the leaderboard and maintain the repository as we keep on doing so. 



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
    
  

