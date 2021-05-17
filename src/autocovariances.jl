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


function find_firm_balanced_set(id, firmid, year)
    df = DataFrame(id = id, firmid = firmid, year = year)

    num_periods = maximum(year) - minimum(year) + 1

    df_firmyear_collapsed = @pipe unique(df, [:firmid, :year]) |> select(_, [:firmid, :year])

    df_firmyear_collapsed = @pipe groupby(df_firmyear_collapsed, :firmid) |> transform(_, nrow => :num_years_present) |> transform(_, :num_years_present => (x -> (x .== num_periods)) => :is_in_balanced)

    new_df = leftjoin(df_firmyear_collapsed, df, on = [:firmid, :year])

    sum(new_df.:is_in_balanced)
    return new_df.:is_in_balanced
end

function leave_out_AR(y, first_id, second_id, year, settings; autocorr_table = false, lags = nothing)

    @unpack obs,  lwage  , id , firmid, controls = get_leave_one_out_set(y, first_id, second_id, settings, controls)
    year = year[obs, :]

    data = DataFrame(lwage = y, id = id, firmidg = firmid, year = year)
    data = @pipe data |> transform(_, [:firmidg, :year] => ((x, y) -> (string.(x, "_", y))) => :firmyearid)
    
    tmp = unique(data.:firmidg)
    data[!, :firmidg] = indexin(data.:firmidg, tmp)

    tmp = unique(unique(data.:firmyearid))
    data[!, :firmyearid] = indexin(data.:firmyearid, tmp)

    tmp = unique(unique(data.:id))
    data[!, :id] = indexin(data.:id, tmp)

    data.:is_in_balanced = find_firm_balanced_set(data.:id, data.:firmidg, data.:year)

    y = data.lwage
    first_id = data.:id
    second_id = data.:firmyearid
    firmid = data.:firmidg
    year = data.:year
    is_in_balanced = data.:is_in_balanced

    
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

    first_year = minimum(data.:year)
    last_year =maximum(data.:year)

    lags = (lags == nothing) ? (1:(last_year - first_year)) : lags

    println("autocorrelation plot")
    #TODO check sanity
    for counter in lags
        if (counter + first_year) <= last_year
            df = @pipe sort(df0, [:firmid, :year]) |> groupby(_, :firmid) |> transform(_, :year => (x -> x .== (lag(x, counter) .+ counter) ) => :has_prev) |> transform(_, :nworkers => (x -> lag(x, counter)) => :nworkers_prev)
            # We need to be sure that the data is totally balanced, otherwise, it may not work
            # df = @pipe sort(df0, [:firmid, :year]) |> groupby(_, :firmid) |> transform(_, :year => (x -> x .== (base_year + counter) ) => :has_prev) |> transform(_, :nworkers => (x -> lag(x, counter)) => :nworkers_prev)
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
            println(θ_lag_cov/(sqrt(θ_lag_var) * sqrt(θ_current_var)))
            print("corr2:")
            println(θ_lag_cov/(sqrt(θ_second) * sqrt(θ_current_var)))
            # println(size(WF, 1))
            println(" ")
        end
    end

    if autocorr_table == true
        for base_year in first_year:last_year
            print("base year:")
            println(base_year)
            for counter in lags
                if (counter + first_year) <= last_year
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
                    println(θ_lag_cov/(sqrt(θ_lag_var) * sqrt(θ_current_var)))
                    print("corr2:")
                    println(θ_lag_cov/(sqrt(θ_second) * sqrt(θ_current_var)))
                    # println(size(WF, 1))
                    println(" ")
                end
            end
        end
    end
end

# data = CSV.read("gen_data_homoskedastic_v03_v01.csv", DataFrame, missingstrings = ["NA", ""])
# data = CSV.read("data_set_firm_balanced_DGP_small_connected_set_no_error.csv", DataFrame, missingstrings = ["NA", ""])
data = CSV.read("first_step_reduced_leave_out_data_renamed.csv", DataFrame, missingstrings = ["NA", ""])
# data = CSV.read("data_set_firm_balanced.csv", DataFrame, missingstrings = ["NA", ""])
# data = CSV.read("data_set_firm_balanced_DGP_large_connected_set.csv", DataFrame, missingstrings = ["NA", ""])

tmp = unique(data.:firmidg)
data[!, :firmidg] = indexin(data.:firmidg, tmp)

tmp = unique(unique(data.:firmyearid))
data[!, :firmyearid] = indexin(data.:firmyearid, tmp)

tmp = unique(unique(data.:id))
data[!, :id] = indexin(data.:id, tmp)

id = data.:id
firmid = data.:firmidg
year = data.:year

data.:is_in_balanced = find_firm_balanced_set(id, firmid, year)
# sum(data.:is_in_balanced)

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

base_year = 1987
counter = 1
for counter in 1:14
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