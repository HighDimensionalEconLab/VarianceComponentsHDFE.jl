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
# data = CSV.read("gen_data.csv", DataFrame, missingstrings = ["NA", ""])
# data = CSV.read("small_data_gen4.csv", DataFrame, missingstrings = ["NA", ""])
# data = CSV.read("reduced_data.csv", DataFrame, missingstrings = ["NA", ""])
# data = CSV.read("gen_data_homoskedastic_m0_v1.csv")
# data = CSV.read("first_step_reduced_leave_out_data_renamed.csv", DataFrame, missingstrings = ["NA", ""])
data = CSV.read("gen_data_homoskedastic_m0_v1.csv", DataFrame, missingstrings = ["NA", ""])

# data = CSV.read("first_step_reduced_leave_out_data_renamed_nge50.csv", DataFrame, missingstrings = ["NA", ""])
# CSV.write("first_step_reduced_leave_out_data_renamed_nge50.csv", data)

# CSV.write("small_data_gen3.csv", data)


### experiments
rename!(data, :firmid => :firmidg)
# data = data[obs, :]
tmp = unique(data.:firmidg)
data[!, :firmidg] = indexin(data.:firmidg, tmp)

tmp = unique(unique(data.:firmyearid))
data[!, :firmyearid] = indexin(data.:firmyearid, tmp)

tmp = unique(unique(data.:id))
data[!, :id] = indexin(data.:id, tmp)

#first we run leave on out function to find leave one out sets. aassume we have them. We then get the corresponding firmids too in addition to firmyearids. For now assume we have them.
# firmid = data[:, "firmidg"]
# firmyearid = data.:firmyearid
# year = data.:year
# # nworkers = data.:nworkers

# personid = data[:, "idg"]
# firmyearid = data[:, "firmyearid"]
# y = data[:, :lwage]
# firmid = data[:, :firmidg]
# year = data[:, :year]

# personid = [1, 1, 2, 2, 1, 2]
# firmyearid = [1, 2, 1, 3, 4, 4, 5, 6, 7, 8]
# firmid = [1, 2, 1, 1, 2, 2, 2, 2, 1, 1]
# year = [1, 2, 1, 2, 3, 3, 4, 5, 6, 7]

df = DataFrame(firmid = data.:firmidg, firmyearid = data.:firmyearid, year = data.:year)
# df = DataFrame(firmid = firmid, firmyearid = firmyearid,  year = year)
df = @pipe groupby(df, [:firmyearid]) |> combine(_, :firmid => first => :firmid, :year => first => :year, nrow => :nworkers)

#we assume that the data is sorted by firmid and year, so the sort in the next line shouldn't affect anything.
df = @pipe sort(df, [:firmid, :year]) |> groupby(_, :firmid) |> transform(_, :year => (x -> x .== (lead(x) .- 1) ) => :has_next) |> transform(_, :nworkers => lead => :nworkers_next)
df[ismissing.(df.:has_next), :has_next] = 0
df[!, :row_number] = 1:nrow(df)

# sum(skipmissing(d3.:has_next))

df2 = df[df.:has_next .== 1, :]
df2[!, :row_number2] = 1:nrow(df2)
M = size(df, 1)
N = size(df2, 1)

A1 = sparse(df2.:row_number2, df2.:row_number, -1)
A1 = hcat(A1, spzeros(N, M - size(A1, 2)))
A2 = sparse(df2.:row_number2, df2.row_number .+ 1, 1)
A2 = hcat(A2, spzeros(N, M - size(A2, 2)))

A = A1 + A2

# sum_weights = sum(df2.:nworkers_next)
# denom = nrow(data)
ncols = nrow(df2)
# numer = 2
# scall = sqrt((numer-1)/(denom-1))
# scall = 1
weights = df2.:nworkers_next

# @time weightedA = vcat([repeat(A[i,:]', inner=(weights[i],1)) for i=1:size(A,1)]...)

# sum_weights = sum(weights)
# weights = df2.:nworkers_next
# W = Diagonal(weights)
# weightedA = scall * W * A

# weights = [3 2 4]
# S = [-1 1 0 0 0; 0 -1 1 0 0; 0 0 0 -1 1]
# ncols = 3
sum_weights = sum(weights)
Wp = spzeros(sum_weights, ncols)
i = 1
j = 1
@time for weight in weights
    Wp[i:i+weight-1, j] .= 1.0
    i += weight
    j += 1
end

# Matrix(Wp)
# weightedA = Wp
weightedA = Wp * A

# rows, cols, vals = findnz(A)

y = data.:lwage
first_id = data.:id
second_id = data.:firmyearid

settings = VCHDFESettings(leverage_algorithm = JLAAlgorithm(num_simulations=0),
    print_level = 1,
    leave_out_level = "obs", 
    first_id_effects = true,
    first_id_bar_effects = true,
    cov_effects = true
    )

## NOW GO run the Leave out estimation function

# y_raw = y
# first_id_raw = first_id
# second_id_raw = second_id

# nrows = size(A, 1)
# A1 = hcat(spzeros(nrows, N), weightedA[:, 1:J-1])

# tmp = Diagonal(1 ./ sum(F, dims = 1)[1,:]) * F'
# A2 = weightedA * tmp
# W
# NT = size(personid, 1)
# J = maximum(firmyearid)
# N = maximum(personid)

#first (worker) Dummies
D = sparse(collect(1:NT),personid,1)
D = Matrix(D)
#second (firm) Dummies
F = sparse(collect(1:NT),firmyearid,1)
F = Matrix(F)
F
F * (F'*F)^(-1) * F'
### experimentsd1

# ommitting the obs that are outliers in term of wage growth for at least one s
# @pipe data |> filter!(x -> x.obs_1 & x.obs_2 & x.obs_3 & x.obs_4 & x.obs_5, _)
data = data[1:100000, :]
first_id_raw = data[:,"id"]
second_id_raw = data[:, "year_by_firm"]
# y_raw = data[:, "log_dailywages"]
y_raw = data[:, "lwage"]



S = sparse( [-1.0 1.0 0 0;0 0.0 -1.0 1.0; ])
w = [3 ; 2 ]  #3 times first row, 2 times second row
Matrix(S)
Srep = vcat([repeat(S[i,:]', inner=(w[i],1)) for i=1:size(S,1)]...)
Matrix(Srep)



@unpack obs,  y  , first_id , second_id, controls = get_leave_one_out_set(y_raw, first_id_raw, second_id_raw, settings, nothing)
data
data = data[obs, :]



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
print(θ_second)
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
