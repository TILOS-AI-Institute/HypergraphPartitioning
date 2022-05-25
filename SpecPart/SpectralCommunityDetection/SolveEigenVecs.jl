function LinearAlgebra.ldiv!(c::AbstractVecOrMat{T}, P::CombinatorialMultigrid.lPreconditioner, b::AbstractVecOrMat{T}) where {T}
	(n, m) = size(b)

    for k in 1:m
        c[:, k] = P.f(b[:, k])
	end
end

@inline function makeBFunc(vWts::Vector, pindex::Pindex, multiplier::Vector)
	bfunc = X -> dCliques(X, vWts, multiplier) + 10*dBiClique(X, pindex)
	return bfunc
end

@inline function makeAFunc(H::Hypergraph)
	afunc = X -> hypL(H, X)
	return afunc
end

function ProjectionStep!(eigvecs::Array{Float64}, degsG::Vector{Float64}, Bfun::Function)
    (n, m) = size(eigvecs)
	o = ones(Float64, n)

	for i in 1:m
		c = -dot(degsG, eigvecs[:, i]) / dot(degsG, o)
		eigvecs[:, i] += c * o
		eigvecs[:, i] /= sqrt(eigvecs[:, i]' * Bfun(eigvecs[:, i]))
	end

	for i in 1:n
		eigvecs[i, :] /= norm(eigvecs[i, :], 2)
	end
end

function GenEigenVecs(H::Hypergraph, A::SparseMatrixCSC, pindex::Pindex, largest::Bool, nev::Int, solver_iters::Int; bmap = nothing)
	vWts = H.vwts
	d = rand(length(vWts)) ./ 1e06
    L = lap(A) + spdiagm(d)
	(~, n) = size(L)
	multiplier = ones(size(vWts, 2))
	Y = rand(n)
	Y = zeros(n)
	Y[pindex.p1] .= -1.0
	Y[pindex.p2] .= 1.0
	afunc = makeAFunc(H)
	amap = LinearMap(afunc, issymmetric=true, n)
	evec = Float64[]
	#=A_g = SimpleGraph(A)
	ii = zeros(Int, n)
	jj = similar(ii)
	vv = similar(ii)

	for i in 1:n
		ii[i] = i
		jj[i] = i
		vv[i] = degree(A_g, i)
	end

	#B = sparse(ii, jj, vv, n, n)
	B = sparse(I, n, n)=#

	if bmap == nothing
		bfunc = makeBFunc(vWts, pindex, multiplier)
		bmap = LinearMap(bfunc, n)
	end
	
	t_solver = @elapsed begin
		if size(L)[1] < 100
			results = lobpcg(L, largest, nev + 1, maxiter=solver_iters)
			evec = results.X[:,2:nev+1]
		else
			(pfunc, hierarchy) = CombinatorialMultigrid.cmg_preconditioner_lap(L)
			#results = lobpcg(amap, bmap, largest, nev + 1, tol=1e-20, maxiter=solver_iters)
			results = lobpcg(amap, bmap, largest, nev, tol=1e-40, maxiter=solver_iters, P=CombinatorialMultigrid.lPreconditioner(pfunc), log=true)
			#results = lobpcg(amap, largest, Y, nev, tol=1e-40, maxiter=solver_iters, P=CombinatorialMultigrid.lPreconditioner(pfunc), log=true)
			#results = lobpcg(amap, bmap, largest, nev, tol=1e-40, maxiter=solver_iters, log=true)
			@info "[EIGEN VECTOR DETAILS] :: $results"
			evec = results.X
		end
	end

	if nev > 1
		#ProjectionStep!(results.X, sum(A, dims=1)[1, :], bfunc)
	end

	@info "[EIGEN VECTOR DETAILS] :: SOLVER TIME :: $t_solver seconds"

	#@info "$results"
	return evec
end
