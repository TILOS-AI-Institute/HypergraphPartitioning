function hmetis(fname::String, Nparts::Int, UBfactor::Int, Nruns::Int, CType::Int, RType::Int, Vcycle::Int, Reconst::Int, dbglvl::Int)
    hmetis_r = `./hmetis $fname $Nparts $UBfactor $Nruns $CType $RType $Vcycle $Reconst $dbglvl`
    run(hmetis_r, wait=true)
end

function RunHmetis(fname::String, NParts::Int, UBfactor::Int, nruns::Int, n::Int)
    j = 0
    Nruns = 10
    CType = 1
    RType = 1
    Vcycle = 1
    Reconst = 0
    bM = zeros(Int, nruns, n)

    for i in 1:nruns
        hmetis(fname, NParts, UBfactor, Nruns, CType, RType, Vcycle, Reconst, 0)
        j = 0
        f = open(fname*".part.2")

        for ln in eachline(f)
            j += 1
            bM[i, j] = parse(Int, ln)
        end

        close(f)
    end

    for i in 1:nruns
        b = bM[i, :]
        area = zeros(Int, 2)
        for j in 1:length(b)
            area[b[j]+1] + 1
        end

        if area[1] > area[2] 
            b = 1 .- b
        end

        bM[i, :] = b
    end

    return bM
end

function IdentifyFixedVtxs(bM::Matrix{Int})
    (m, ~) = size(bM)
    psum = sum(bM, dims=1)
    psum = psum[1, :]

    v1 = findall(x-> x==0, psum)
    v2 = findall(x-> x==m, psum)

    if isempty(v1) && isempty(v2)
        mp = mean(psum)
        stdev = std(psum)

        v2 = findall(x-> x > mp+stdev, psum)
        v1 = findall(x-> x < mp-stdev, psum)
    end

    return v1, v2
end
