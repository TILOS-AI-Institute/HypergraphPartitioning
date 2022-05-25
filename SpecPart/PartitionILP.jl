import JuMP.MathOptInterface.TerminationStatusCode
import JuMP.MathOptInterface.TerminationStatus

function PartitionILP(H::Hypergraph, fixed::Vector{Int}, capacities::Vector{Int})
    k = 0
    n = H.n
    m = H.e
    vwts = H.vwts
    hwts = H.hwts
    N = 2*m
    f = zeros(N)
    tvwt = sum(vwts)
    model = Model(with_optimizer(Cbc.Optimizer, logLevel = 0))
    @variable(model, x[1:n], Bin)
    @variable(model, y[1:N], Bin)
    @constraint(model, c1, x'vwts <= capacities[2])
    @constraint(model, c2, x'vwts >= tvwt-capacities[1])

    for i in 1:n
        if fixed[i] > -1
            @constraint(model, x[i] == fixed[i])
        end
    end

    for i in 1:m
        uloc = H.eptr[i]
        vloc = H.eptr[i+1]-1
        k += 1
        for j in uloc:vloc
            v = H.hedges[j]
            @constraint(model, y[k] - x[v] <= 0)
            @constraint(model, x[v] - y[k+1] <= 0)
        end
        k += 1
    end

    f[1:2:N] = -hwts
    f[2:2:N] = hwts
    @objective(model, Min, y'f)
    optimize!(model)

    if termination_status(model) == TerminationStatusCode(2)
        return -ones(Int, n)
    end

    @show solution_summary(model)

    return Int.(round.(value.(x)))
end
