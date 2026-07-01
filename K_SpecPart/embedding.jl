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

# Relative weight of the supervision (bi-clique) term against the vertex-mass
# (clique) term in the generalized eigenproblem's B operator.
const BICLIQUE_SUPERVISION_WEIGHT = 500

# LOBPCG requires the right-hand operator B to be symmetric POSITIVE DEFINITE.
# Both `clique` and `bi_clique` are only positive *semi*-definite (they
# annihilate the constant vector), and for large k or empty/degenerate
# supervision blocks B becomes singular -> LOBPCG's CholQR throws
# "matrix is not positive definite". Adding a small diagonal term `reg*I` makes
# B strictly SPD. `reg` is increased by the retry loop in `solve_eigs` if needed.
const B_REGULARIZATION = 1e-3

@inline function make_b_func(vWts::Vector,
                        pindex::__pindex__,
                        multiplier::Vector,
                        reg::Float64 = B_REGULARIZATION)
	bfunc = X -> clique(X, vWts, multiplier) +
	             BICLIQUE_SUPERVISION_WEIGHT * bi_clique(X, pindex) +
	             reg .* X
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

    # Small problems: dense symmetric eigensolve on the (regularized) Laplacian.
    if size(lap_matrix, 1) < 100
        results = lobpcg(lap_matrix, largest, nev+1, maxiter=solver_iters)
        return results.X[:, 2:nev+1]
    end

    (pfunc, hierarchy) = CombinatorialMultigrid.cmg_preconditioner_lap(lap_matrix)

    # Try the generalized eigenproblem A x = λ B x with increasing B
    # regularization. If LOBPCG still fails (e.g. "matrix is not positive
    # definite" from a degenerate/ill-conditioned B), fall back to the
    # unsupervised Laplacian eigenvectors so the run continues instead of
    # aborting.
    for reg in (B_REGULARIZATION, 1e-1, 1.0, 10.0)
        local_bmap = bmap
        if local_bmap === nothing
            bfunc = make_b_func(hgraph.vwts, pindex, multiplier, reg)
            local_bmap = LinearMap(bfunc, hgraph.num_vertices)
        end
        try
            results = lobpcg(amap, local_bmap, false, nev,
                             tol=1e-40, maxiter=solver_iters,
                             P = CombinatorialMultigrid.lPreconditioner(pfunc),
                             log = true)
            @debug "LOBPCG finished" reg results
            return results.X
        catch e
            @warn "LOBPCG failed; retrying with larger B regularization" reg exception=e
            bmap === nothing || rethrow(e)  # caller-supplied B: don't silently re-reg
        end
    end

    @warn "LOBPCG could not converge a generalized eigenproblem; " *
          "falling back to unsupervised Laplacian eigenvectors"
    results = lobpcg(lap_matrix, false, nev+1, maxiter=solver_iters)
    return results.X[:, 2:nev+1]
end