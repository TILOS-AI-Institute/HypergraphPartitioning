@inline combine(x, y) = x

function knuthshuffle!(v::AbstractVector)
    for i in length(v):-1:2
        j = rand(1:i)
        v[i], v[j] = v[j], v[i]
    end
end

function HypergraphToGraph(H::Hypergraph, nt::Int)
	n = H.n
	e = H.e
	hedges = H.hedges
	eptr = H.eptr
	s = eptr[end]
	i = zeros(Int, 2*nt*s)
	j = zeros(Int, 2*nt*s)
	w = zeros(Float64, 2*nt*s)
	w_e = zeros(Float64, 2*nt*s)
	r = 1
	hsizes = H.eptr[2:end] - H.eptr[1:end-1]
	max_pin_size = maximum(hsizes)
    ewts = H.hwts

	for k in 1:e
		ew = ewts[k]
		u_loc = eptr[k]
		v_loc = eptr[k+1]
		T = v_loc-u_loc
		if T == 2
			i[r] = hedges[u_loc]
			j[r] = hedges[u_loc+1]
			w[r] = ew
			w_e[r] = ew/max_pin_size #hsizes[k]/max_pin_size
			r += 1
		elseif T == 3
			i[r] = hedges[u_loc]
			j[r] = hedges[u_loc+1]
			w[r] = ew/2.0
			w_e[r] = ew/max_pin_size #hsizes[k]/max_pin_size

			r += 1
			i[r] = hedges[u_loc+1]
			j[r] = hedges[u_loc+2]
			w[r] = ew/2.0
			w_e[r] = ew/max_pin_size #hsizes[k]/max_pin_size
			r += 1
			i[r] = hedges[u_loc+2]
			j[r] = hedges[u_loc]
			w[r] = ew/2.0
			w_e[r] = ew/max_pin_size #hsizes[k]/max_pin_size
			r += 1
		else
			#mult = ceil(T/2)
			#cw = 1/(nt*mult)
			mult = (floor(T/2)*ceil(T/2))/(T-1)
			#mult = T-1
			cw = 1/(nt*2*mult)
			edge = hedges[u_loc:u_loc+T-1]
            p = Vector{Int}(1:T)

			for t in 1:nt
				p = randperm(T)
                #knuthshuffle!(p)
                i[r:r+T-2] = edge[p[1:end-1]]
                j[r:r+T-2] = edge[p[2:end]]
                w[r:r+T-2] .= cw*ew
                w_e[r:r+T-2] .= ew/max_pin_size #hsizes[k]/max_pin_size

                r = r + T-1
                i[r] = edge[p[end]]
                j[r] = edge[p[1]]
                w[r] = cw*ew
                w_e[r] =  ew/max_pin_size #hsizes[k]/max_pin_size
                r += 1
			end
		end
	end

	r -= 1
	i = i[1:r]
	j = j[1:r]
	w = w[1:r]	
	w_e = w_e[1:r]

	n = max(maximum(i), maximum(j))
	A = sparse(i, j, w, n, n)
	A_e = sparse(i, j, w_e, n, n)

	return (A + A'), (A_e + A_e')
end	
