using Test; using Laplacians; using MAT; using LinearAlgebra; using BenchmarkTools
include("../src/CombinatorialMultigrid.jl");

using Main.CombinatorialMultigrid

## load example matrix

file = matopen("../example/X.mat"); X = read(file, "X"); close(file)
LX = lap(X);
b1 = rand(Float64, size(X, 1));
b1 = b1 .- sum(b1)/length(b1);

@info "###  Running CMG ###"
t = @elapsed (pfunc, h) = cmg_preconditioner_lap(LX);
@info "Time Required to build CMG Solver: $(t) seconds"
t = @elapsed x = pfunc(b1);
@info "Time Required to find x: $(t) seconds"

@info "###  Running PCG with CMG    ### "
## solve with pcg and cmg preconditioner
f1 = pcgSolver(LX,pfunc);
t = @elapsed x = f1(b1, maxits=40, tol=1e-6,verbose=true);
@info "Time Required to solve system: $(t) seconds"

## solve with approxchol_lap from laplacians
@info "###  Running approxchol_lap    ###"
t = @elapsed solver = approxchol_lap(X; tol=1e-6, maxits=1000, maxtime=Inf, verbose=true, pcgIts=Int[], params=ApproxCholParams());
@info "Time Required to build Lap Solver: $(t) seconds"
t = @elapsed x = solver(b1);
@info "Time Required to find x: $(t) seconds"