using DataFrames

"""
    pVector(N, vals; Bernoulli = false)

Use `vals` to create a `pvector`` of rank `N`
`Bernoulli=false`

"""
function pVector(N, vals; Bernoulli = false)
    ## 0. Check
    (length(vals) == N * 4 ) || error("vals length not equal $(N*4)")

    ## 1. Create Rows and Columns
    res = Bernoulli ? [0, 1] : [-1, 1]
    inputMatrix = allcombinations(DataFrame; :ID => 1:N,
        [Pair(Symbol("R$(x)"), res) for x in 1:2]...)
    sort!(inputMatrix, [:ID, order(:R1, rev=true), order(:R2, rev=true)])
    transform!(inputMatrix,
        :ID => ByRow(x -> "R$(x)_")               => :ConteXt,
        :ID => ByRow(x -> "R_$(x)")               => :ConteNt,
        :ID => ByRow(x -> "R$(x)$(x)")            => :RXX,
        :ID => ByRow(x -> "R$(x)$(oplus(x, N))")  => :RXP,
        :ID => ByRow(x -> "R$(ominus(x, N))$(x)") => :RMX
        )
    select!(inputMatrix, [:ConteXt, :ConteNt], Not(:ID))
    insertcols!(inputMatrix, "Prob" => vals)

    ## 2. Calculate Values
    pMatrix = combine(groupby(inputMatrix, [:ConteXt, :ConteNt, :RXX, :RXP, :RMX]),
        [:R1, :R2, :Prob] => ((x, y, z) ->        x'z) => :RXXV,
        [:R1, :R2, :Prob] => ((x, y, z) ->        y'z) => :RXPV,
        [:R1, :R2, :Prob] => ((x, y, z) -> (x .* y)'z) => :RXXXP
    )

    function CouplProb(RXP, RMX, RXXV, RXPV)
        vvs = Vector{Number}[]
        for (rxxi, rxxv) in enumerate(RXXV)
            rmxi = findall(==(RMX[rxxi]), RXP)[]
            rmxv = RXPV[rmxi]
            push!(vvs, [rxxv, rmxv])
        end
        return vvs
    end

    transform!(pMatrix, [:RXP, :RMX, :RXXV, :RXPV] => (CouplProb) => :RXXMX)

    ## 3. Return
    return inputMatrix, pMatrix
end

"""
    SP(N, vals; Bernoulli = false)

return splited `p` vector, i.e., pl, pb, pc

"""
function SP(N, vals; Bernoulli = false)
    pv = pVector(N, vals; Bernoulli)
    pl = [1; pv[2][:, "RXXV"]; pv[2][:, "RXPV"]] # Low-order marginals
    pb = pv[2][:, "RXXXP"]                       # Bunch-probabilities
    pc = map(minimum, pv[2][:, "RXXMX"])         # Connection Probabilities
    return pl, pb, pc
end
