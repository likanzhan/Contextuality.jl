using DataFrames

"""
    mMatrix(N)

Create a M matrix of rank N
"""
function mMatrix(N)
    # 1. Low marginals
    RXX = ["R$(x)$(x)" for x in 1:N]
    RXP = ["R$(x)$(oplus(x, N))" for x in 1:N]
    ps = [Pair(Symbol(col), [1, 0]) for col in [RXX; RXP]]
    MMatrix = allcombinations(DataFrame; ps...)
    insertcols!(MMatrix, 1, :R00 => 1)
    sort!(MMatrix, rev = true)

    # 2. Bunch products
    for x in 1:N
        R1 = "R" * string(x) * string(x)
        RP = "R" * string(x) * string(oplus(x, N))
        transform!(MMatrix,
        [Symbol(R1), Symbol(RP)] => ByRow((x, y) -> x==+1 && y==+1 ? 1 : 0) => Symbol(R1*RP)
        )
    end

    # 3. Coupling products
    for x in 1:N
        R1 = "R" * string(x) * string(x)
        RM = "R" * string(ominus(x, N)) * string(x)
        transform!(MMatrix,
        [Symbol(R1), Symbol(RM)] => ByRow((x, y) -> x==+1 && y==+1 ? 1 : 0) => Symbol(R1*RM)
        )
    end

    # 4. return
    MMatrix
end

"""
    SM(N)

Generate Separate M Matrix, i.e., 
Ml, Mc, Mb
"""
function SM(N)
    mm = mMatrix(N)
    mmname = names(mm)

    fml(x) = (length(x) == 3)                   # Low-order marginals
    fmb(x) = (length(x) == 6) && (x[2] == x[5]) # Bunch probability matrix
    fmc(x) = (length(x) == 6) && (x[3] == x[6]) # Connection Probability matrix

    mlname = mmname[fml.(mmname)]
    mbname = mmname[fmb.(mmname)]
    mcname = mmname[fmc.(mmname)]

    mldf = mm[:, mlname]
    mbdf = mm[:, mbname]
    mcdf = mm[:, mcname]

    Ml = permutedims(Matrix(mldf))
    Mb = permutedims(Matrix(mbdf))
    Mc = permutedims(Matrix(mcdf))
    
    return Ml, Mb, Mc
end