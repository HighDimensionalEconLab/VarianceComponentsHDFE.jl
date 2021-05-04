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

data = CSV.read("gen_data_homoskedastic_v03_v01_psibar.csv", DataFrame, missingstrings = ["NA", ""])

### some data cleaning
# rename!(data, :firmid => :firmidg)

tmp = unique(data.:firmidg)
data[!, :firmidg] = indexin(data.:firmidg, tmp)

tmp = unique(unique(data.:firmyearid))
data[!, :firmyearid] = indexin(data.:firmyearid, tmp)

tmp = unique(unique(data.:id))
data[!, :id] = indexin(data.:id, tmp)

### In the following lines, we find W and Fdelta
df = DataFrame(firmid = data.:firmidg, firmyearid = data.:firmyearid, year = data.:year)
df = @pipe groupby(df, [:firmyearid]) |> combine(_, :firmid => first => :firmid, :year => first => :year, nrow => :nworkers)

#we assume that the data is sorted by firmid and year, so the sort in the next line shouldn't affect anything.
df = @pipe sort(df, [:firmid, :year]) |> groupby(_, :firmid) |> transform(_, :year => (x -> x .== (lead(x) .- 1) ) => :has_next) |> transform(_, :nworkers => lead => :nworkers_next)
df[ismissing.(df.:has_next), :has_next] = 0
df[!, :row_number] = 1:nrow(df)

df2 = df[df.:has_next .== 1, :]
df2[!, :row_number2] = 1:nrow(df2)
M = size(df, 1)
N = size(df2, 1)

A1 = sparse(df2.:row_number2, df2.:row_number, -1)
A1 = hcat(A1, spzeros(N, M - size(A1, 2)))
A2 = sparse(df2.:row_number2, df2.row_number .+ 1, 1)
A2 = hcat(A2, spzeros(N, M - size(A2, 2)))

Fdelta = A1 + A2

ncols = nrow(df2)
weights = df2.:nworkers_next
sum_weights = sum(weights)
W = spzeros(sum_weights, ncols)
i = 1
j = 1
@time for weight in weights
    W[i:i+weight-1, j] .= 1.0
    i += weight
    j += 1
end

### Finishing finding W and Fdelta
WFdelta = W * Fdelta

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

#We call a modified version of leave_out_estimation, so if you want to call it, you should run "] instantiate" again. Or just run that function in REPL line by line
leave_out_estimation(y,first_id,second_id,nothing,settings, WFdelta)
