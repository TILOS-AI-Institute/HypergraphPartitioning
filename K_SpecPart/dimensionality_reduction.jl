using MultivariateStats

function pca(mat::Array{Float64}, labels::Vector{Int})
    mat_t = mat'
    pca = fit(PCA, mat_t; maxoutdim=2)
    Ypca = predict(pca, mat_t)
    return Ypca
end

function lda(mat::Array{Float64}, labels::Vector{Int}, num_parts::Int)
    mat_t = mat'
    # MulticlassLDA produces at most (num_classes - 1) discriminative axes, and
    # cannot exceed the input dimension. Use all available axes (was hard-coded
    # to 2) so K-way tree construction gets a richer embedding for K > 3.
    # Set KSPECPART_LDA_FULL=0 to recover the old fixed 2-dimensional behavior.
    full_dims = get(ENV, "KSPECPART_LDA_FULL", "1") == "1"
    maxdim = full_dims ? num_parts - 1 : 2
    outdim = clamp(maxdim, 1, size(mat, 2))
    # MulticlassLDA throws when a class is empty/singleton (common for large k
    # with sparse hints). Fall back to the leading embedding columns so the run
    # continues rather than aborting.
    try
        model = fit(MulticlassLDA, mat_t, labels; outdim=outdim)
        return predict(model, mat_t)
    catch e
        @warn "MulticlassLDA failed (degenerate classes?); using leading embedding columns" exception=e
        return permutedims(mat[:, 1:outdim])
    end
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