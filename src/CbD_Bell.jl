"""
    Bell(N, pvals)

Calculate modifed CHSH value
- Dzhafarov, E. N., & Kujala, J. V. (2016). Context–content systems of random variables: The contextuality-by-default theory. Journal of Mathematical Psychology, 74, 11-33. https://doi.org/10.1016/j.jmp.2016.04.010
- Dzhafarov, E. N., Kujala, J. V., & Cervantes, V. H. (2020). Contextuality and noncontextuality measures and generalized Bell inequalities for cyclic systems. Physical Review A, 101(4). https://doi.org/10.1103/PhysRevA.101.042119 (Page 5)
- Kujala, J. V., & Dzhafarov, E. N. (2019). Measures of contextuality and non-contextuality. Philosophical Transactions of the Royal Society A: Mathematical Physical and Engineering Sciences, 377(2157), 20190149. https://doi.org/10.1098/rsta.2019.0149 

"""
function Bell(N, pvals)
    # 1. Calculate p values
    pv = pVector(N, pvals; Bernoulli = false)

    # 2. Calcualte s1eb
    eb = pv[2][:, "RXXXP"]
    s1eb = SODD(eb)

    # 3. Calculate δ
    el = pv[2][:, "RXXMX"]
    δel = mapreduce(x->abs(diff(x)[]), +, el)

    # 4. Calculate Δ
    Δ = min(N - 2 + δel, N)
    KD = s1eb - Δ
    return [s1eb KD]
end