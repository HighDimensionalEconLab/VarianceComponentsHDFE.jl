# using DocStringExtensions: groupby
using CSV
using Parameters
using SparseArrays
using LinearAlgebra
using DataFrames
using Pipe
using ShiftedArrays
using Random
using Statistics
using VarianceComponentsHDFE
using DocStringExtensions

controls = nothing
lags = []
autocorr_plot = false
settings = VCHDFESettings(leverage_algorithm = JLAAlgorithm(num_simulations=200),
print_level = 2,
leave_out_level = "obs"
)

# data_raw = CSV.read("data_set_firm_balanced_churn_rate_split.csv", DataFrame)

# data_raw = CSV.read("DGP_data_churn_rates_balanced_error3.csv", DataFrame)
data_raw = CSV.read("data/data_set_firm_balanced_churn_rate_split_reduced.csv", DataFrame)
# data_raw = CSV.read("very_small_hand_G.csv", DataFrame)
# time_id = data_raw.tvar

y = data_raw.lwage
first_id = data_raw.id
second_id = data_raw.firmidg
time_id = data_raw.year

flag_vector = data_raw.cr_ge_09

#Create time-second_id identifier (naming it second_id) and rename the second_id to firmid as it is the second id in the time varying AKM model
data = @pipe DataFrame(y = y, first_id = first_id, firmid = second_id, time_id = time_id, flag1 = flag_vector) |> sort(_, [:firmid, :time_id])

data.:second_id = compute_matchid(data.:time_id, data.:firmid)
# data = @pipe data |> transform(_, [:firmid, :time_id] => ((x, y) -> (string.(x, "_", y))) => :second_id)

tmp = unique(data.:firmid)
data[!, :firmid] = indexin(data.:firmid, tmp)

tmp = unique(data.:second_id)
data[!, :second_id] = indexin(data.:second_id, tmp)

tmp = unique(data.:first_id)
data[!, :first_id] = indexin(data.:first_id, tmp)

y = data.:y
first_id = data.:first_id
second_id = data.second_id

#Find the leave-out connected set
tmp = get_leave_one_out_set(y, first_id, second_id, settings, controls)
obs = tmp.:obs

#Limiting data to the leave out connected set
data = data[obs, :]

tmp = unique(data.:firmid)
data[!, :firmid] = indexin(data.:firmid, tmp)

tmp = unique(data.:second_id)
data[!, :second_id] = indexin(data.:second_id, tmp)

tmp = unique(data.:first_id)
data[!, :first_id] = indexin(data.:first_id, tmp)

#Find the balanced set
data.:is_in_balanced = find_balanced_set(data.:first_id, data.:firmid, data.:time_id; settings)

data = sort(data, [:first_id, :time_id])

y = (data.:y)
first_id = (data.:first_id)
second_id = data.:second_id
firmid = data.:firmid
time_id = data.:time_id
is_in_balanced = data.:is_in_balanced
flag_vector = data.:flag1


NT = size(y,1)
J = maximum(second_id)
N = maximum(first_id)
K = controls == nothing ? 0 : size(controls,2)
nparameters = N + J + K
D = sparse(collect(1:NT),first_id,1)
F = sparse(collect(1:NT),second_id,1)
S= sparse(1.0I, J-1, J-1)
S=vcat(S,sparse(-zeros(1,J-1)))

X = hcat(D, -F*S)
Dvar = hcat(  D, spzeros(NT,J-1) )
Fvar = hcat(spzeros(NT,N), -F*S )




@time @unpack θ_first, θ_second, θCOV, β, Dalpha, Fpsi, Pii, Bii_first, Bii_second, Bii_cov, y, X, sigma_i = leave_out_estimation(y,first_id,second_id,controls,settings)
#save("first_round.jld2",  Dict("θ_first" => θ_first, "θ_second" => θ_second, "θCOV" => θCOV, "β" => β, "Dalpha" => Dalpha, "Fpsi" => Fpsi, "Pii" => Pii, "Bii_first" => Bii_first, "Bii_second" => Bii_second, "Bii_cov" => Bii_cov, "y" => y, "X" => X, "sigma_i" => sigma_i))
# tmp = load("first_round.jld2")
# θ_first, θ_second, θCOV, β, Dalpha, Fpsi, Pii, Bii_first, Bii_second, Bii_cov, y, X, sigma_i = tmp["θ_first"], tmp["θ_second"], tmp["θCOV"], tmp["β"], tmp["Dalpha"], tmp["Fpsi"], tmp["Pii"], tmp["Bii_first"], tmp["Bii_second"], tmp["Bii_cov"], tmp["y"], tmp["X"], tmp["sigma_i"]
beta = β
df0 = DataFrame(firmid = firmid, firmyearid = second_id, year = time_id, is_in_balanced = is_in_balanced, flag1 = flag_vector)
df0 = @pipe DataFrames.groupby(df0, :firmyearid) |> combine(_, :firmid => first => :firmid, :year => first => :year, nrow => :nworkers, :is_in_balanced => maximum => :is_in_balanced, :flag1 => maximum => :flag1) |> sort(_, [:firmid, :year])

first_year = minimum(df0.:year)
last_year =maximum(df0.:year)

lags = (isempty(lags)) ? (1:(last_year - first_year)) : lags

acp = nothing


acf = Array{Union{Missing, Float64}}(missing, last_year - first_year + 1, last_year - first_year + 1)
covss = Array{Union{Missing, Float64}}(missing, last_year - first_year + 1, last_year - first_year + 1)
varss = Array{Union{Missing, Float64}}(missing, last_year - first_year + 1, last_year - first_year + 1)
varlagss = Array{Union{Missing, Float64}}(missing, last_year - first_year + 1, last_year - first_year + 1)
for t in 1:(last_year - first_year + 1)
    acf[t, t] = 1.0::Float64
end

## Finding weights to use later:
#TODO
df_weights = @pipe df0[df0.is_in_balanced .== 1, :] |> DataFrames.groupby(_, :firmid) |> combine(_, :nworkers => mean => :nworkers)
df_weights.:nworkers = floor.(Int, df_weights.:nworkers)
weights = df_weights.:nworkers
sum_weights = sum(weights)
W = spzeros(sum_weights, size(weights, 1))
i = 1
j = 1
for weight in weights
    W[i:i+weight-1, j] .= 1.0
    i += weight
    j += 1
end


# for base_year in (first_year):last_year
for base_year=first_year:last_year
    # (settings.print_level > 0) && @info "base year:" base_year
    for counter in lags
        # (settings.print_level > 1) && @info typeof(counter)
        if (counter + base_year) <= last_year
            # df = @pipe sort(df0, [:firmid, :year]) |> groupby(_, :firmid) |> transform(_, :year => (x -> x .== (lag(x, counter) .+ counter) ) => :has_prev; ungroup = false) |> transform(_, :nworkers => (x -> lag(x, counter)) => :nworkers_prev)
            df = @pipe df0 |> DataFrames.groupby(_, :firmid) |> transform(_, :year => (x -> x .== (base_year + counter) ) => :has_prev; ungroup = false) |> transform(_, :nworkers => (x -> lag(x, counter)) => :nworkers_prev)
            df[ismissing.(df.:has_prev), :has_prev] .= 0
            df[!, :row_number] = 1:nrow(df)


            # df2 = df[(df.:has_prev .== 1) .& (df.:is_in_balanced .== 1) .& (df.:flag1 .== 1), :]
            df2 = df[(df.has_prev .== 1) .& (df.:is_in_balanced .== 1), :]

            df2[!, :row_number2] = 1:nrow(df2)
            # df2[!, :nWorkers_mean] = floor.(Int, (df2.:nworkers .+ df2.:nworkers_prev) ./ 2)
            M = size(df, 1)
            L = size(df2, 1)

            Flag = sparse(df2.:row_number2, df2.:row_number .- counter, 1)
            Flag = hcat(Flag, spzeros(L, M - size(Flag, 2)))

            Flag0 = sparse(df2.:row_number2, df2.:row_number, 1)
            Flag0 = hcat(Flag0, spzeros(L, M - size(Flag0, 2)))

            ##### Now creating W and WFlag1
            # ncols = nrow(df2)
            # weights = df2.:nWorkers_mean
            # sum_weights = sum(weights)
            # W = spzeros(sum_weights, ncols)
            # i = 1
            # j = 1
            # for weight in weights
            #     W[i:i+weight-1, j] .= 1.0
            #     i += weight
            #     j += 1
            # end

            WFlag = W * Flag
            WF = W * Flag0

            Fvar = hcat(spzeros(size(WF, 1), N), WF[:, 1:J-1])
            FlagVar = hcat(spzeros(size(WFlag, 1), N), WFlag[:, 1:J-1])

            if isempty(Fvar) || isempty(FlagVar)
                @warn "1"
                continue
            end

            @time @unpack Bii_lag_var, Bii_current_var, Bii_lag_cov = leverages3(X, Fvar, FlagVar, settings)

            # θ_lag_var = kss_quadratic_form(sigma_i[is_in_balanced .== 1, :], FlagVar, FlagVar, beta, Bii_lag_var[is_in_balanced .== 1])
            # println(Bii_current_var)
            if size(Bii_current_var, 1) == 0 || isempty(Fvar) || isempty(FlagVar) || isempty(Bii_current_var) 
                @warn "2"
                continue
            end

            # data2 = @pipe sort(data, [:firmid, :time_id]) |> groupby(_, :firmid) |> transform(_, :time_id => (x -> x .== (base_year + counter) ) => :has_prev)

            # selector =  (data.:is_in_balanced .== true ) .| (data.is_in_balanced .== false) 
            # df2 = df[(df.:has_prev .== 1) .& (df.:is_in_balanced .== 1) .& (df.:flag1 .== 1), :] .& (data.:flag1 .== 1)
            # selector = (data.:time_id .== base_year + counter) 

            θ_lag_cov = kss_quadratic_form(sigma_i, Fvar, FlagVar, beta, zeros(size(Bii_current_var, 1)))
            θ_current_var = kss_quadratic_form(sigma_i, Fvar, Fvar, beta, Bii_current_var)
            θ_lag_var = kss_quadratic_form(sigma_i, FlagVar, FlagVar, beta, Bii_lag_var)
            sizes = size(Fvar, 1) 

            # θ_lag_cov = kss_quadratic_form(zeros(size(sigma_i, 1))[selector, :], Fvar, FlagVar, beta, zeros(size(Bii_current_var, 1))[selector, :])
            # θ_current_var = kss_quadratic_form(sigma_i[selector], Fvar, Fvar, beta, Bii_current_var[selector])
            # θ_lag_var = kss_quadratic_form(sigma_i[selector, :], FlagVar, FlagVar, beta, Bii_lag_var[selector])
            # dof = size(sigma_i[selector, :], 1)  
            # sizes = size(Fvar, 1) 

            cor_kss_corrected = θ_lag_cov/(sqrt(θ_lag_var) * sqrt(θ_current_var))
            settings.print_level > 0 && @info "autocorrelation between $(base_year) and $(base_year+counter):" cov=θ_lag_cov  var=θ_current_var var_lag=θ_lag_var correlation=cor_kss_corrected sizes=sizes


            acf[base_year - first_year + 1, base_year - first_year + 1 + counter ] = cor_kss_corrected
            varss[base_year - first_year + 1, base_year - first_year + 1 + counter ] = θ_current_var
            varlagss[base_year - first_year + 1, base_year - first_year + 1 + counter ] = θ_lag_var
            covss[base_year - first_year + 1, base_year - first_year + 1 + counter ] = θ_lag_cov
        end
    end
end

