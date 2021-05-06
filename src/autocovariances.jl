using CSV
using VarianceComponentsHDFE
using Parameters
using SparseArrays
using LinearAlgebra
using DataFrames
using Pipe
using ShiftedArrays
using LinearMaps
using Random
using Statistics

# data = CSV.read("gen_data_homoskedastic_v03_v01_psibar.csv", DataFrame, missingstrings = ["NA", ""])
# data = CSV.read("small_data_gen4.csv", DataFrame, missingstrings = ["NA", ""])
data = CSV.read("first_step_reduced_leave_out_data_renamed.csv", DataFrame, missingstrings = ["NA", ""])

### some data cleaning
# rename!(data, :firmid => :firmidg)

tmp = unique(data.:firmidg)
data[!, :firmidg] = indexin(data.:firmidg, tmp)

tmp = unique(unique(data.:firmyearid))
data[!, :firmyearid] = indexin(data.:firmyearid, tmp)

tmp = unique(unique(data.:id))
data[!, :id] = indexin(data.:id, tmp)

# ################################################
# #We will find Flag
# df = DataFrame(firmid = data.:firmidg, firmyearid = data.:firmyearid, year = data.:year)
# df = @pipe groupby(df, [:firmyearid]) |> combine(_, :firmid => first => :firmid, :year => first => :year, nrow => :nworkers)

# df = @pipe sort(df, [:firmid, :year]) |> groupby(_, :firmid) |> transform(_, :year => (x -> x .== (lag(x) .+ 1) ) => :has_prev_1) |> transform(_, :nworkers => lag => :nworkers_prev)
# df[ismissing.(df.:has_prev_1), :has_prev_1] = 0
# df[!, :row_number] = 1:nrow(df)


# df2 = df[df.:has_prev_1 .== 1, :]
# df2[!, :row_number2] = 1:nrow(df2)
# M = size(df, 1)
# N = size(df2, 1)

# Flag1 = sparse(df2.:row_number2, df2.:row_number .- 1, 1)
# Flag1 = hcat(Flag1, spzeros(N, M - size(Flag1, 2)))

# Flag0 = sparse(df2.:row_number2, df2.:row_number, 1)
# Flag0 = hcat(Flag0, spzeros(N, M - size(Flag0, 2)))

# ##### Now creating W and WFlag1
# ncols = nrow(df2)
# weights = df2.:nworkers_prev
# sum_weights = sum(weights)
# W = spzeros(sum_weights, ncols)
# i = 1
# j = 1
# @time for weight in weights
#     W[i:i+weight-1, j] .= 1.0
#     i += weight
#     j += 1
# end

# W2 = spzeros(sum_weights, ncols)
# i = 1
# j = 1
# @time for weight in weights
#     W2[i:i+weight-1, j] .= 1.0
#     i += weight
#     j += 1
# end

# WFlag1 = W * Flag1
# WF = W2 * Flag0

y = data.:lwage
first_id = data.:id
second_id = data.:firmyearid


settings = VCHDFESettings(leverage_algorithm = JLAAlgorithm(num_simulations=0),
    print_level = 1,
    leave_out_level = "obs", 
    first_id_effects = false,
    first_id_bar_effects = false,
    cov_effects = false, 
    )

#We call a modified version of leave_out_estimation, so if you want to call it, you should run "] instantiate" again. Or just run that function in REPL line by line
# leave_out_estimation(y,first_id,second_id,nothing,settings, WFdelta)
firmid = data.:firmidg
year = data.:year
leave_out_AR(y, first_id, second_id, data.:firmidg, data.:year, settings)