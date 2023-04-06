using GraphLaplacians
using LinearMaps
using IterativeSolvers

include("cmg/CombinatorialMultigrid.jl")

function LinearAlgebra.ldiv!(c::AbstractVecOrMat{T}, 
                            P::CombinatorialMultigrid.lPreconditioner, 
                            b::AbstractVecOrMat{T}) where {T}
	(n, m) = size(b)
    for k in 1:m
        c[:, k] = P.f(b[:, k])
	end
end

@inline function make_b_func(vWts::Vector, 
                        pindex::__pindex__, 
                        multiplier::Vector)
	bfunc = X -> clique(X, vWts, multiplier) + 500*bi_clique(X, pindex)
	return bfunc
end

@inline function make_a_func(H::__hypergraph__, epsilon::Int)
	afunc = X -> hypl(H, X, epsilon)
	return afunc
end

function projection_step!(eigs::Array{Float64},
                        degs::Vector{Float64},
                        bfunc::Function)
    (n, m) = size(eigs)
    o = ones(Float64, n)
    for i in 1:m
        c = -dot(degs, eigs[:, i]) / dot(degs, o)
        eigs[:, i] += c*o
        eigs[:, i] /= sqrt(eigs[:, i]' * bfunc(eigs[:, i]))
    end
    for i in 1:n
        eigs[i, :] /= norm(eigs[i, :], 2)
    end
end

function solve_eigs(hgraph::__hypergraph__,
                    adj::SparseMatrixCSC,
                    pindex::__pindex__,
                    largest::Bool,
                    nev::Int,
                    solver_iters::Int,
                    bmap = nothing;
                    epsilon::Int = 1)
    d = ones(hgraph.num_vertices) ./ 1e06
    degs = GraphLaplacians.degree_matrix(adj)
    lap_matrix = spdiagm(d) + degs - adj
    multiplier = ones(size(hgraph.vwts, 2))
    afunc = make_a_func(hgraph, epsilon)
    amap = LinearMap(afunc, issymmetric=true, hgraph.num_vertices)
    if bmap == nothing
        bfunc = make_b_func(hgraph.vwts, pindex, multiplier)
        bmap = LinearMap(bfunc, hgraph.num_vertices)
    end
    evecs = Float64[]
    if size(lap_matrix, 1) < 100
        results = lobpcg(lap_matrix, largest, nev+1, maxiter=solver_iters)
        evecs = results.X[:, 2:nev+1]
    else 
        (pfunc, hierarchy) = CombinatorialMultigrid.cmg_preconditioner_lap(lap_matrix)
        results = lobpcg(amap, 
                        bmap, 
                        false, 
                        nev, 
                        tol=1e-40, 
                        maxiter=solver_iters, 
                        P = CombinatorialMultigrid.lPreconditioner(pfunc), 
                        log = true)
        evecs = results.X
        line_log = repeat("=", 60)
        @info "$line_log"
        @info "$results"
        @info "$line_log"
    end
 
    #=projection_step!(evecs, 
                    GraphLaplacians.degrees(adj), 
                    bfunc)=#
    return evecs
end