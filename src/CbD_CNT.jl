using JuMP, HiGHS

"""
    CNT2(N, pvals)

- Kujala, J. V., & Dzhafarov, E. N. (2019). Measures of contextuality and non-contextuality. Philosophical Transactions of the Royal Society A: Mathematical Physical and Engineering Sciences, 377(2157), 20190149. https://doi.org/10.1098/rsta.2019.0149
"""
function CNT2(N, pvals)
    # 1. Calculate M and P
    Ml, Mb, Mc = SM(N)
    pl, pb, pc = SP(N, pvals)
    dlength, xlength = size(Mb)

    # 2. Optimize
    cnt2 = Model(HiGHS.Optimizer)
    @variable(cnt2, x[1:xlength] >= 0)
    @variable(cnt2, d[1:dlength]   >= 0)
    @constraint(cnt2, (pb - Mb * x) .>= -d)
    @constraint(cnt2, (pb - Mb * x) .<= +d)
    @constraint(cnt2, Ml * x .== pl)
    @constraint(cnt2, Mc * x .== pc)
    @objective(cnt2, Min, sum(d))
    optimize!(cnt2)

    # 3. Results
    is_solved_and_feasible(cnt2) || error("Fail to Converge")
    if is_solved_and_feasible(cnt2)
        xvs = collect(value.(x))
        return mapreduce(abs, +, pb - Mb * xx)
    end
end
