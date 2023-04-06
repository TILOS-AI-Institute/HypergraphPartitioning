# K-SpecPart #

K-SPecPart is an extension (and improvement) over SpecPart that can compute multi-way (K > 2) partitions. 

# Requirements #

In order to run SpecPart, the following requirements are mandatory: 
* [Julia](https://julialang.org/) version 1.6.0 or newer. 
* The following Julia libraries:
  * [Shuffle](https://docs.juliahub.com/Shuffle/X0eqg/0.1.1/)
  * [LightGraphs](https://github.com/sbromberger/LightGraphs.jl)
  * [SimpleWeightedGraphs](https://github.com/JuliaGraphs/SimpleWeightedGraphs.jl)
  * [LDLFactorizations](https://github.com/JuliaSmoothOptimizers/LDLFactorizations.jl)
  * [IterativeSolvers](https://iterativesolvers.julialinearalgebra.org/stable/)
  * [Laplacians](https://github.com/danspielman/Laplacians.jl)
  * [LinearMaps](https://github.com/JuliaLinearAlgebra/LinearMaps.jl)
* [CPLEX ILP solver](https://www.ibm.com/support/pages/downloading-ibm-ilog-cplex-optimization-studio-v1290) 
* [METIS](https://github.com/KarypisLab/METIS)
* [hMETIS](http://glaros.dtc.umn.edu/gkhome/metis/hmetis/overview)

Julia libraries can be installed by doing the following in Julia REPL:

```
julia> using Pkg
julia> Pkg.add("Package name")
```
