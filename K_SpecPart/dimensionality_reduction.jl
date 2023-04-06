using MultivariateStats

function pca(mat::Array{Float64}, labels::Vector{Int})
    mat_t = mat'
    pca = fit(PCA, mat_t; maxoutdim=2)
    Ypca = predict(pca, mat_t)
    return Ypca
end

function lda(mat::Array{Float64}, labels::Vector{Int})
    mat_t = mat'
    lda = fit(MulticlassLDA, mat_t, labels; outdim=2)
    Ylda = predict(lda, mat_t)
end

function write_embedding_to_file(mat::Array{Float64}, embedding_file::String)
    f = open(embedding_file, "w")
    for i in 1:size(mat, 2)
        for j in 1:size(mat, 1)
            print(f, mat[j, i], " ")
        end
        print(f, "\n")
    end
    close(f)
end

function write_labels_to_file(labels::Vector{Int}, labels_file::String)
    f = open(labels_file, "w")
    for i in 1:length(labels)
        println(f, labels[i])
    end
    close(f)
end