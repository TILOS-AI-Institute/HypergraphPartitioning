
# CombinatorialMultigrid.jl
This package implements the Combinatorial Multigrid Preconditioner. Refer to *[1]* for details on the algorithm. 


In order to run CMG we present a quick example. Lets load a a million sized example matrix ```X``` from the ```example``` directory and build the ```b``` side. ```X``` is derived from a proprietary hypergraph. 

```
## load example matrix
file = matopen("../example/X.mat"); X = read(file, "X"); close(file)
LX = lap(X);
b = rand(Float64, size(X, 1));
b = b .- sum(b)/length(b);
```
CMG needs to be built before solving a linear system. To build CMG we provide two wrapper functions: ```cmg_preconditioner_lap```, which requires the user to provide the laplacian matrix and ```cmg_preconditioner_adj``` which requires the user to provide the adjacent matrix. 
We run the following script:

```
t = @elapsed (pfunc, h) = cmg_preconditioner_lap(LX);
@info "Time Required to build CMG Solver: $(t) seconds"
t = @elapsed x = pfunc(b);
@info "Time Required to find x: $(t) seconds"
```
Both ```cmg_preconditioner_lap``` and ```cmg_preconditioner_adj``` returns two parameters: the solver function and the hierarchy.
The above script generates the following output: 
```
[ Info: Time Required to build CMG Solver: 3.258120718 seconds
[ Info: Time Required to find x: 0.194163587 seconds
```
We try to solve a linear system using CMG. For this purpose we leverage ```pcg``` from the ```Laplacians``` package. We run the following script: 
```
f = pcgSolver(LX,pfunc);
t = @elapsed x = f(b, maxits=40, tol=1e-6,verbose=true);
@info "Time Required to solve system: $(t) seconds"
```
This generates the following output: 
```
PCG stopped after: 7.828 seconds and 29 iterations with relative error 8.429320186485909e-7.
[ Info: Time Required to solve system: 8.034237971 seconds
```

For comparison we run ```approxchol_lap``` which is the fastest solver from the ```Laplacians``` package. We run the following script: 
```
solver = approxchol_lap(X; tol=1e-6, maxits=1000, maxtime=Inf, verbose=true, pcgIts=Int[], params=ApproxCholParams());
@info "Time Required to build Lap Solver: $(t) seconds"
t = @elapsed x = solver(b);
@info "Time Required to find x: $(t) seconds"
```
This generates the following output: 
```
Using greedy degree ordering. Factorization time: 25.985280990600586
Ratio of operator edges to original edges: 4.548086602013639
ratio of max to min diagonal of laplacian : 583453.0510646806
Solver build time: 26.293 seconds.
[ Info: Time Required to build Lap Solver: 26.304803977 seconds
PCG stopped after: 10.288 seconds and 26 iterations with relative error 9.697886904926194e-7.
[ Info: Time Required to find x: 12.226966502 seconds
```

```CMG``` builds the solver in ```3.26 seconds``` compared to ```30 seconds``` with ```approxchol_lap``` and solves ```x``` in ```0.19 seconds``` compared to ```12.23 seconds```. When we run ```pcg``` using the two solvers as preconditioners, we find ```cmg``` yields better performance compared to ```approxchol_lap``` by solving the linear system in 7.8 seconds and 10.29 seconds respectively. 


**Citations:**

[1] Ioannis Koutis, Gary L. Miller, David Tolliver, Combinatorial preconditioners and multilevel solvers for problems in computer vision and image processing, Computer Vision and Image Understanding, Volume 115, Issue 12, 2011, Pages 1638-1646, ISSN 1077-3142, https://doi.org/10.1016/j.cviu.2011.05.013.*
