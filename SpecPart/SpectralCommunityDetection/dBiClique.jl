function dBiClique(X::AbstractArray, pindex::Pindex)
    n = length(X)
    y = zeros(n)
    d1 = ones(length(pindex.p1))
    d2 = ones(length(pindex.p2))
    #y[pindex.p1] = (sum(d2) .* d1 .* X[pindex.p1] - d1 * (d2' * X[pindex.p2])) 
    #y[pindex.p2] = (sum(d1) .* d2 .* X[pindex.p2] - d2 * (d1' * X[pindex.p1])) 
    
    t1 = Threads.@spawn (sum(d2) .* d1 .* X[pindex.p1] - d1 * (d2' * X[pindex.p2])) 
    t2 = Threads.@spawn (sum(d1) .* d2 .* X[pindex.p2] - d2 * (d1' * X[pindex.p1])) 

    t1 = fetch(t1)
    t2 = fetch(t2)

    y[pindex.p1] = t1
    y[pindex.p2] = t2

    return y
end
