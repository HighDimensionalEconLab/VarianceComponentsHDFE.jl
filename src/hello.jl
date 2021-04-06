# In this file, we omit all observations that have their wage growth 
# for s = 1, 2, 3, 4, 5 in the top or the bottom percentiles. 
# Which makes it easier and faster to do experiments in the resulted dataset
using CSV
using VarianceComponentsHDFE
using Parameters
using SparseArrays
using LinearAlgebra
using DataFrames
using Pipe
using ShiftedArrays
using LinearMaps

# First line: the whole data, second and third lines: SSA data
# also uncomment the write section at the ends


# data = CSV.read("layoff_raw_data_mlayoff_dropped.csv", DataFrame, missingstrings = ["NA", ""])
# data = CSV.read("layoff_raw_data_mlayoff_dropped_SSA_B.csv", DataFrame, missingstrings = ["NA", ""])
# data = CSV.read("layoff_raw_data_mlayoff_dropped_SSA_B.csv", DataFrame, missingstrings = ["NA", ""])
# data = CSV.read("py_final_1985_2001_veneto_only_added_vars.csv", DataFrame, missingstrings = ["NA", ""])
# data = CSV.read("SSA_data_group_A_True.csv", DataFrame, missingstrings = ["NA", ""])
data = CSV.read("gen_data.csv", DataFrame, missingstrings = ["NA", ""])
# data = CSV.read("SSA_data_group_B_False.csv")

# ommitting the obs that are outliers in term of wage growth for at least one s
# @pipe data |> filter!(x -> x.obs_1 & x.obs_2 & x.obs_3 & x.obs_4 & x.obs_5, _)
data = data[1:100000, :]
first_id_raw = data[:,"id"]
second_id_raw = data[:, "year_by_firm"]
# y_raw = data[:, "log_dailywages"]
y_raw = data[:, "lwage"]

settings = VCHDFESettings(leverage_algorithm = JLAAlgorithm(num_simulations=0),
    print_level = 1,
    leave_out_level = "obs"
    )

@unpack obs,  y  , first_id , second_id, controls = get_leave_one_out_set(y_raw, first_id_raw, second_id_raw, settings, nothing)

print(size(data,1) == size(y, 1))

# using JLD2
# @save "tmp.jld2" obs y first_id second_id controls settings
# # @save "tmp2.jld2" obs y first_id second_id controls settings
# # @load "tmp.jld2" obs y first_id second_id controls y_raw first_id_raw second_id_raw settings
# @load "tmp3.jld2" obs y first_id second_id controls 
# # using LinearOperators
# # Dbarvar = hcat(F * opCholesky((F' *F)) * F' * D , spzeros(NT,J-1) )

@time @unpack θ_first, θ_second, θCOV, θ_firstbar, β, Dalpha, Fpsi, Pii, Bii_first, Bii_second, Bii_cov, Bii_first_bar, y, X, sigma_i = leave_out_estimation(y,first_id,second_id,controls,settings)

print(θ_firstbar)
print(θCOV)
# using JLD2
# @save "obs.jld2" obs
# @save "y.jld2" y
# @save "first_id.jld2" first_id
# @save "second_id.jld2" second_id

NT = size(y,1)
J = maximum(second_id)
N = maximum(first_id)
nparameters = N + J 
#first (worker) Dummies
D = sparse(collect(1:NT),first_id,1)
#second (firm) Dummies
F = sparse(collect(1:NT),second_id,1)
# N+J x N+J-1 restriction matrix
S= sparse(1.0I, J-1, J-1)
S=vcat(S,sparse(-zeros(1,J-1)))
X = hcat(D, -F*S)
#SET DIMENSIONS
n=size(X,1)
# ESTIMATE HIGH DIMENSIONAL MODEL
xx=X'*X
xy=X'*y
compute_sol = approxcholSolver(xx;verbose = settings.print_level > 0)
beta = compute_sol([xy...];verbose=false)

pe=D * beta[1:N]
fe=F*S * beta[N+1:N+J-1]

data_leaveout = data[obs,:]
data_leaveout
data_2 = hcat(data_leaveout, pe)
rename!(data_2, "x1" => "person_effects")
data_3 = hcat(data_2, -fe)
rename!(data_3, "x1" => "year_by_firm_effects")

CSV.write("first_step_reduced_leave_out_data_layoff_SSA_B.csv", data_3)

CSV.write("new.csv", data_3)
# CSV.write("SSA_group_A_reduced_leave_out_data.csv", data_3)
# CSV.write("SSA_group_B_reduced_leave_out_data.csv", data_3)


@unpack θ_first, θ_second, θCOV, β, Dalpha, Fpsi, Pii, Bii_first, Bii_second, Bii_cov, y, X, sigma_i = leave_out_estimation(y,first_id,second_id,controls,settings)

# @unpack θ_first, θ_second, θCOV, β, Dalpha, Fpsi, Pii, Bii_first, Bii_second, Bii_cov = leave_out_estimation(y,first_id,second_id,controls,settings)
