function readPartitions(partition_path::String)
    j = 0
    listdir = readdir(partition_path, join=true)
    dir = listdir[1]
    listdir_x = readdir(dir)
    fname = listdir_x[findall(y-> y == match(r"([\w\d\_.]+.hgr)", y).captures[1], listdir_x)]
    fname = fname[1]
    fname = listdir[1] * "/" * fname
    (hedges, eptr, vwts, n, hE, pins) = readHmetisBenchmark(fname)
    #bn = match(r"(ibm\d+)", fname) 
    #list_dir = readdir(pname * bn.captures[1], join=true)
    bins = zeros(Int, length(listdir), n)
    #pf_r = match(r"/(ISPD98_ibm[0-9]+.weighted.hgr)", fname)
    pf_r = match(r"([\w\d_.]+.hgr)", fname)
    pf_name = pf_r.captures[1] * ".part.2"

    for i in 1:length(listdir)
        j = 0
        dir_name = listdir[i]
        pf = dir_name * "/" * pf_name
        f = open(pf)

        for ln in eachline(f)
            j += 1
            bins[i, j] = parse(Int, ln)
        end

        close(f)
    end

    return (bins, hedges, eptr, vwts, n, hE, pins)
end

function prunePartitions(nsolns::Int, bins::Matrix{Int}, H::Hypergraph, B::Incidence)
    (m, n) = size(bins)
    cutlist = Dict{Int, Vector{Int}}()
    pids = Vector{Int}()

    for i in 1:m
        bin = bins[i, :]
        sum_bin = sum(bin)
        cutsize = findCutsize(bin, H, B)
        
        if haskey(cutlist, cutsize) == true
            ptns = cutlist[cutsize]
            flag = false

            for j in 1:length(ptns)
                bin_j = bins[ptns[j], :]
                sum_bin_j = sum(bin_j)
            
                if sum_bin == sum_bin_j && sum_bin == n-sum_bin_j
                    flag = true
                end
            end
            
            if flag == false
                push!(cutlist[cutsize], i)
            end
        else
            push!(cutlist, cutsize => [i])
        end
    end

    cutlist_n = sort(collect(pairs(cutlist)), by=x->x[1])

    for i in 1:length(cutlist_n)
        append!(pids, cutlist_n[i][2])

        if length(pids) >= nsolns
            break
        end
    end

    return (cutlist_n, pids)
end

function processBestPartitions(bins::Matrix{Int}, pids::Vector{Int})
    (m, n) = size(bins)
    fixed = -ones(Int, n)
    vwts = ones(Int, n)
    bins_new = zeros(Int, length(pids), n)

    for i in 1:length(pids)
        area_split = [0, 0]
        p = pids[i]
        bin = bins[p, :]

        for j in 1:n
            area_split[bin[j]+1] += vwts[j]
        end

        if area_split[1] > area_split[2]
            bin = 1 .- bin
        end

        bins_new[i, :] = bin
    end

    vtx_history = sum(bins_new, dims=1)[1, :]
        
    for i in 1:n
        if vtx_history[i] == m || vtx_history[i] == 0
            fixed[i] = bins[1, i]
        end
    end

    return (bins_new, fixed)
end