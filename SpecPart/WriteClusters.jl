function WriteClusters(community::Vector{Int})
    f = open("clusters.dat", "w")
    
    for i in 1:length(community)
        println(f, i, ",", community[i])
    end

    close(f)
end