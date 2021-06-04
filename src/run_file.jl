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

data_raw = CSV.read("data_set_firm_balanced_churn_rate_split.csv", DataFrame)

y = data_raw.lwage
first_id = data_raw.id
second_id = data_raw.firmidg
time_id = data_raw.year

flag_vector = data_raw.cr_ge_09

#Create time-second_id identifier (naming it second_id) and rename the second_id to firmid as it is the second id in the time varying AKM model
data = DataFrame(y = y, first_id = first_id, firmid = second_id, time_id = time_id, flag1 = flag_vector)

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




@unpack θ_first, θ_second, θCOV, β, Dalpha, Fpsi, Pii, Bii_first, Bii_second, Bii_cov, y, X, sigma_i = leave_out_estimation(y,first_id,second_id,controls,settings)

beta = β
df0 = DataFrame(firmid = firmid, firmyearid = second_id, year = time_id, is_in_balanced = is_in_balanced, flag1 = flag_vector)
df0 = @pipe groupby(df0, [:firmyearid]) |> combine(_, :firmid => first => :firmid, :year => first => :year, nrow => :nworkers, :is_in_balanced => maximum => :is_in_balanced, :flag1 => maximum => :flag1)

first_year = minimum(df0.:year)
last_year =maximum(df0.:year)

lags = (isempty(lags)) ? (1:(last_year - first_year)) : lags

acp = nothing


acf = Array{Union{Missing, Float64}}(missing, last_year - first_year + 1, last_year - first_year + 1)
for t in 1:(last_year - first_year + 1)
    acf[t, t] = 1.0::Float64
end

for base_year in (first_year):last_year
    (settings.print_level > 0) && @info "base year:" base_year
    for counter in lags
        (settings.print_level > 1) && @info typeof(counter)
        if (counter + base_year) <= last_year
            df = @pipe sort(df0, [:firmid, :year]) |> groupby(_, :firmid) |> transform(_, :year => (x -> x .== (base_year + counter) ) => :has_prev) |> transform(_, :nworkers => (x -> lag(x, counter)) => :nworkers_prev)
            df[ismissing.(df.:has_prev), :has_prev] .= 0
            df[!, :row_number] = 1:nrow(df)


            df2 = df[(df.:has_prev .== 1) .& (df.:is_in_balanced .== 1) .& (df.:flag1 .== 1), :]
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
            for weight in weights
                W[i:i+weight-1, j] .= 1.0
                i += weight
                j += 1
            end

            WFlag = W * Flag
            WF = W * Flag0

            Fvar = hcat(spzeros(size(WF, 1), N), WF[:, 1:J-1])
            FlagVar = hcat(spzeros(size(WFlag, 1), N), WFlag[:, 1:J-1])

            if isempty(Fvar) || isempty(FlagVar)
                @warn "1"
                continue
            end

            @unpack Bii_lag_var, Bii_current_var, Bii_lag_cov = leverages3(X, Fvar, FlagVar, settings)

            # θ_lag_var = kss_quadratic_form(sigma_i[is_in_balanced .== 1, :], FlagVar, FlagVar, beta, Bii_lag_var[is_in_balanced .== 1])
            # println(Bii_current_var)
            if size(Bii_current_var, 1) == 0 || isempty(Fvar) || isempty(FlagVar) || isempty(Bii_current_var) 
                @warn "2"
                continue
            end
            θ_lag_cov = kss_quadratic_form(zeros(size(sigma_i, 1))[(is_in_balanced .== 1) .& (flag_vector .== 1), :], Fvar, FlagVar, beta, zeros(size(Bii_current_var, 1))[(is_in_balanced .== 1) .& (flag_vector .== 1), :])
            θ_current_var = kss_quadratic_form(sigma_i[(is_in_balanced .== 1) .& (flag_vector .== 1), :], Fvar, Fvar, beta, Bii_current_var[(is_in_balanced .== 1) .& (flag_vector .== 1)])
            θ_lag_var = kss_quadratic_form(sigma_i[(is_in_balanced .== 1) .& (flag_vector .== 1), :], FlagVar, FlagVar, beta, Bii_lag_var[(is_in_balanced .== 1) .& (flag_vector .== 1)])
            sizes = size(sigma_i[(is_in_balanced .== 1) .& (flag_vector .== 1), :], 1)

            cor_kss_corrected = θ_lag_cov/(sqrt(θ_lag_var) * sqrt(θ_current_var))
            settings.print_level > 0 && @info "autocorrelation between $(base_year) and $(base_year+counter):" cov=θ_lag_cov  var=θ_current_var var_lag=θ_lag_var correlation=cor_kss_corrected sizes=sizes


            acf[base_year - first_year + 1, base_year - first_year + 1 + counter ] = cor_kss_corrected
        end
    end
end

