function hypL(H::Hypergraph, x::AbstractArray)
    v = H.hedges
    loc = H.eptr
    N = length(x)
    m = length(loc)-1
    y = zeros(Float64, N)
    w_ = H.hwts

    for j in 1:m
        k = loc[j+1] - loc[j]
        #mult = ceil(k/2)
        mult = (floor(k/2)*ceil(k/2))/(k-1)
        sm = 0.0
        
        for t in loc[j]:loc[j+1]-1
            sm += x[v[t]]
        end

        s = sm/k

        for t in loc[j]:loc[j+1]-1
            ind = v[t]
            y[ind] += w_[j]*(x[ind] - s)/mult
        end
    end

    return y
end
