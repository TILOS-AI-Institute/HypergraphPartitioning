using MultivariateStats
using Statistics

function db_index(embedding::Matrix{Float64}, labels::Vector{Int})
    dims = size(embedding, 1)
    (min_id, max_id) = extrema(labels)
    num_clusters = min_id == 0 ? max_id+1 : max_id
    centroid_arr = zeros(dims, num_clusters)
    sigma_arr = zeros(num_clusters)
    # calculate centroids and sigmas
    for i in min_id:max_id
        points = findall(labels .== i)
        coords = labels[:, points]
        centroid = mean(coords, dims=2)
        centroid_arr[:, i] = centroid
        sigma = coords .- centroid
        sigma_eucl = sqrt(sigma'sigma)
        sigma_arr[i] = sigma
    end
    db_arr = Float64[]
    for i in min_id:max_id
        db_local = Float64[]
        for j in i+1:max_id
            centroid_diff = centroid_arr[:, i] - centroid_arr[:, j]
            centroid_diff_eucl = sqrt(centroid_diff'centroid_diff)
            push!(db_local, (sigma_arr[i] + sigma_arr[j])/centroid_diff_eucl)
        end
        push!(db_arr, maximum(db_local))
    end
    dbi = mean(db_arr)
    return dbi
end

function vrc_index(embedding::Matrix{Float64}, labels::Vector{Int})
    
end

function SC(points::Matrix{Float64}, labels::Vector{Int})

end

function eval_cluster(mat::Array{Float64}, labels::Vector{Int})
    mat_t = mat'
    pca = fit(PCA, mat_t; maxoutdim=2)
    Ypca = predict(pca, mat_t)
    lda = fit(MulticlassLDA, mat_t, labels; outdim=2)
    Ylda = predict(lda, mat_t)

end