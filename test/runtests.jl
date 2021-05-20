using VarianceComponentsHDFE
using Test



#run_full_test = get(ENV, "VCHDFE_FULL_TEST", "1") == "1" ? true : false
#include("prepare_test_data.jl")

#@testset "VarianceComponentsHDFE.jl" begin
#    include("test_matrices.jl")
#    run_full_test && include("fulltest.jl")
#end

include("test_matrices.jl")
