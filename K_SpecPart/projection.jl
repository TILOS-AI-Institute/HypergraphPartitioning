function projection(evec::Matrix)
    total_vecs = size(evec, 2)
    total_iters = total_vecs/2
    side_0 = zeros(size(evec, 1))
    side_1 = zeros(size(evec, 1))
    iter = 0
    for i in 1:total_vecs
        if i%2 == 0
            side_1 += evec[:,i]
        else
            side_0 += evec[:,i]
        end
    end
    #side_0 ./= total_iters
    #side_1 ./= total_iters
    return hcat(side_0, side_1)
end

function dimensionality_reduction(evec::Matrix, target_dims::Int, seed::Int)
    Random.seed!(seed)
    n = size(evec, 2);
    rand_projection_matrix = [rand((-1, 1)) for i=1:n, j=1:target_dims]
    return evec * rand_projection_matrix
end