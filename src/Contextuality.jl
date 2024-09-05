module Contextuality

include("CbD_utils.jl")
export oplus, ominus, SODD

include("CbD_M.jl")
export mMatrix, SM

include("CbD_P.jl")
export pVector, SP

include("CbD_Bell.jl")
export Bell

include("CbD_CNT.jl")
export CNT2

end
