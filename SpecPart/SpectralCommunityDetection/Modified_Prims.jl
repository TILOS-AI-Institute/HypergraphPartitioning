function modified_prim_mst end
@traitfn function modified_prim_mst(g::AG::(!IsDirected), adj_mat::SparseMatrixCSC, 
    distmx::AbstractMatrix{T}=weights(g)) where {T <: Real, U, AG <: AbstractGraph{U}}
    
    nvg = nv(g)

    pq = PriorityQueue{U, T}()
    finished = zeros(Bool, nvg)
    wt = fill(typemax(T), nvg) #Faster access time
    parents = zeros(U, nv(g))

    pq[1] = typemin(T)
    wt[1] = typemin(T)

    while !isempty(pq)
        v = dequeue!(pq)
        finished[v] = true

        for u in neighbors(g, v)
            finished[u] && continue
            
            if wt[u] > distmx[u, v]
                wt[u] = distmx[u, v] + adj_mat[u, v]
                pq[u] = wt[u]
                parents[u] = v
            end
        end
    end

    return [Edge{U}(parents[v], v) for v in vertices(g) if parents[v] != 0]
end

