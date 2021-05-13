using CSV
using Parameters
using SparseArrays
using LinearAlgebra
using DataFrames
using Pipe
using ShiftedArrays
using LinearMaps
using Random
using Statistics
using VarianceComponentsHDFE

# data = CSV.read("gen_data_homoskedastic_v03_v01_psibar.csv", DataFrame, missingstrings = ["NA", ""])
# data = CSV.read("data_set_firm_balanced_DGP_small_connected_set_no_error.csv", DataFrame, missingstrings = ["NA", ""])
# data = CSV.read("first_step_reduced_leave_out_data_renamed.csv", DataFrame, missingstrings = ["NA", ""])
data = CSV.read("data_set_firm_balanced.csv", DataFrame, missingstrings = ["NA", ""])
# data = CSV.read("data_set_firm_balanced_DGP_large_connected_set.csv", DataFrame, missingstrings = ["NA", ""])


### some data cleaning
# rename!(data, :firmid => :firmidg)

tmp = unique(data.:firmidg)
data[!, :firmidg] = indexin(data.:firmidg, tmp)

tmp = unique(unique(data.:firmyearid))
data[!, :firmyearid] = indexin(data.:firmyearid, tmp)

tmp = unique(unique(data.:id))
data[!, :id] = indexin(data.:id, tmp)


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

is_in_balanced = data.:is_in_balanced

# @unpack obs,  y  , first_id , second_id, controls = get_leave_one_out_set(y, first_id, second_id, settings, nothing)

# df = data[obs, :]
# CSV.write("data_set_firm_balanced_DGP_small_connected_set.csv", df)

# leave_out_AR(y, first_id, second_id, data.:firmidg, data.:year, settings)

NT = size(y,1)
J = maximum(second_id)
N = maximum(first_id)
nparameters = N + J
D = sparse(collect(1:NT),first_id,1)
F = sparse(collect(1:NT),second_id,1)
S= sparse(1.0I, J-1, J-1)
S=vcat(S,sparse(-zeros(1,J-1)))

X = hcat(D, -F*S)
Dvar = hcat(  D, spzeros(NT,J-1) )
Fvar = hcat(spzeros(NT,N), -F*S )

X = hcat(D, -F*S)

@unpack Pii , Mii  , correction_JLA , Bii_first , Bii_second , Bii_cov, Bii_first_bar, Bii_dif_cov, Bii_dif_first_bar, Bii_dif_second_bar = leverages2(settings.leverage_algorithm, X, Fvar, Fvar, Fvar, Fvar, Fvar, settings)
xx=X'*X
xy=X'*y
compute_sol = approxcholSolver(xx;verbose = settings.print_level > 0)
beta = compute_sol([xy...];verbose=false)
eta=y-X*beta

eta_h = eta ./ Mii
sigma_i = ( ( y .- mean(y) ) .* eta_h ) .* correction_JLA

θ_second = kss_quadratic_form(zeros(size(sigma_i, 1)), Fvar, Fvar, beta, Bii_second)
θ_second = kss_quadratic_form(sigma_i, Fvar, Fvar, beta, Bii_second)

# Creating matrixes
df0 = DataFrame(firmid = firmid, firmyearid = second_id, year = year, is_in_balanced = is_in_balanced)
df0 = @pipe groupby(df0, [:firmyearid]) |> combine(_, :firmid => first => :firmid, :year => first => :year, nrow => :nworkers, :is_in_balanced => maximum => :is_in_balanced)

base_year = 1986
counter = 1
for counter in 1:15
    # df = @pipe sort(df, [:firmid, :year]) |> groupby(_, :firmid) |> transform(_, :year => (x -> x .== (lag(x, counter) .+ counter) ) => :has_prev) |> transform(_, :nworkers => (x -> lag(x, counter)) => :nworkers_prev)
    # We need to be sure that the data is totally balanced, otherwise, it may not work
    df = @pipe sort(df0, [:firmid, :year]) |> groupby(_, :firmid) |> transform(_, :year => (x -> x .== (base_year + counter) ) => :has_prev) |> transform(_, :nworkers => (x -> lag(x, counter)) => :nworkers_prev)
    df[ismissing.(df.:has_prev), :has_prev] = 0
    df[!, :row_number] = 1:nrow(df)


    df2 = df[df.:has_prev .== 1 .& df.:is_in_balanced .== 1, :]
    df2[!, :row_number2] = 1:nrow(df2)
    df2[!, :nWorkers_mean] = floor.(Int, (df2.:nworkers .+ df2.:nworkers_prev) ./ 2)
    M = size(df, 1)
    L = size(df2, 1)

    Flag = sparse(df2.:row_number2, df2.:row_number .- counter, 1)
    Flag = hcat(Flag, spzeros(L, M - size(Flag, 2)))

    Flag0 = sparse(df2.:row_number2, df2.:row_number, 1)
    Flag0 = hcat(Flag0, spzeros(L, M - size(Flag0, 2)))

    ##### Now creating W and WFlag1
    ncols = nrow(df2)
    weights = df2.:nWorkers_mean
    sum_weights = sum(weights)
    W = spzeros(sum_weights, ncols)
    i = 1
    j = 1
    @time for weight in weights
        W[i:i+weight-1, j] .= 1.0
        i += weight
        j += 1
    end

    WFlag = W * Flag
    WF = W * Flag0


    Fvar = hcat(spzeros(size(WF, 1), N), WF[:, 1:J-1])
    FlagVar = hcat(spzeros(size(WFlag, 1), N), WFlag[:, 1:J-1])

    @time @unpack Bii_lag_var, Bii_current_var, Bii_lag_cov = leverages3(X, Fvar, FlagVar, settings)

    # θ_lag_cov = kss_quadratic_form(sigma_i[is_in_balanced, :], Fvar[is_in_balanced, :], FlagVar[is_in_balanced, :], beta[is_in_balanced], Bii_lag_cov[is_in_balanced])
    # θ_lag_var = kss_quadratic_form(sigma_i[is_in_balanced, :], FlagVar[is_in_balanced, :], FlagVar[is_in_balanced, :], beta[is_in_balanced], Bii_lag_var[is_in_balanced])
    θ_lag_var = kss_quadratic_form(sigma_i[is_in_balanced .== 1, :], FlagVar, FlagVar, beta, Bii_lag_var[is_in_balanced .== 1])
    θ_lag_cov = kss_quadratic_form(zeros(size(sigma_i, 1))[is_in_balanced .== 1, :], Fvar, FlagVar, beta, zeros(size(Bii_lag_var, 1))[is_in_balanced .== 1, :])
    θ_current_var = kss_quadratic_form(sigma_i[is_in_balanced .== 1, :], Fvar, Fvar, beta, Bii_current_var[is_in_balanced .== 1])
    println(counter, ":")
    
    print("cov: ")
    println(θ_lag_cov)
    print("var lag: ")
    println(θ_lag_var)
    print("var current: ")
    println(θ_current_var)
    print("corr:")
    println(θ_lag_cov/(sqrt(θ_lag_var) * sqrt(θ_second)))
    # println(size(WF, 1))
    println(" ")
end