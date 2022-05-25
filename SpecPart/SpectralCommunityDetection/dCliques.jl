function dCliques(X::AbstractArray, vWts::Vector, multiplier::AbstractArray)  
    T = sum(vWts)
    n = size(X, 1)
    y = zeros(Float64, n)
    
    s = multiplier[1]/sum(vWts)
    K = vWts'X
    sd = sum(vWts)

    Threads.@threads for j in 1:n
        y[j] += T * ((vWts[j] * X[j]) - ((K * vWts[j])/sd))*s
    end
    
    return y
end 
