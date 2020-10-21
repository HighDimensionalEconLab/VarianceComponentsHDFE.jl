using VarianceComponentsHDFE
using Random, BenchmarkTools, Test
using JLD, SparseArrays, LinearAlgebra
# Make sure LDLFactorizations is version 0.5.0 and that the `multiple-rhs` branch is checked out

pkg_dir = pkgdir(VarianceComponentsHDFE)

# Obtain environment variables
run_large_benchmark = get(ENV, "VCHDFE_LARGE_BENCHMARK", "false")  == "true" ? true : false
num_rhs_large = try
        parse(Int, get(ENV, "BENCHMARK_NUM_RHS_LARGE", "2"))
    catch
        2
    end

# Entry point for the PkgBenchmarks call.  Can split into different files later.
include("prepare_benchmark_data.jl")
# NOTE: Suite below can assume that the `benchmark/data/...` has been filled

Random.seed!(1234)


# 1) Medium sized network benchmarks 

medium_data = load(pkg_dir*"/benchmark/data/medium_main.jld")

X = medium_data["X_GroundedLaplacian"]
S_xx = medium_data["S_xx"]
S_xx_sparse = sparse(S_xx)

#Loading Operator
ldli, la = computeLDLinv(S_xx_sparse)
buffs = zeros(size(la)[1],1)
P = approxcholOperator(ldli, buffs[:,1])
compute_sol = approxcholSolver(P,la)

#RHS 
m,k = size(X)
R_p = convert(Array{Float64,2}, bitrand(1,m))
rademacher!(R_p)
JLA_RHS = (R_p*X)[1,:]

#Benchmarks 
@btime compute_sol($JLA_RHS;verbose=false)

#2) Large network benchmark 
if run_large_benchmark

large_data = load(pkg_dir*"/benchmark/data/large_main.jld")

X = large_data["X_GroundedLaplacian"]
S_xx = large_data["S_xx"]
S_xx_sparse = sparse(S_xx)

#Loading Operator
ldli, la = computeLDLinv(S_xx_sparse)
buffs = zeros(size(la)[1],1)
P = approxcholOperator(ldli, buffs[:,1])
compute_sol = approxcholSolver(P,la)

#RHS 
m,k = size(X)
R_p = convert(Array{Float64,2}, bitrand(1,m))
rademacher!(R_p)
JLA_RHS = (R_p*X)[1,:]

#Benchmarks 
@btime compute_sol($JLA_RHS;verbose=false)


end