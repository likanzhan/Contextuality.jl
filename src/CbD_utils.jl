using Combinatorics

"""
    oplus(a, N)

"""
function oplus(a, N)
    a = isa(a, Number) ? a : parse(Int, a)
    a + 1 > N ? 1 : a + 1
end

"""
    ominus(a, N)

"""
function ominus(a, N)
    a = isa(a, Number) ? a : parse(Int, a)
    a - 1 < 1 ? N : a - 1
end

"""
    SODD(xs)

Similar to `CHSH`, but more compelx.

"""
function SODD(xs)
    n = length(xs)
    mv = []
    for cs in 1:2:n
        for ngs in combinations(1:n, cs)
            lo = ones(Int, n)
            lo[ngs] .= -1
            push!(mv, lo)
        end
    end
    mat = Matrix(hcat(mv...)')
    vals = mat * xs
    maximum(vals)
end