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

"""
$(SIGNATURES)

Returns a vector indicating if an observation belongs to the second_id-balanced set. 

### Arguments
* `first_id`: first identifier (e.g. worker id)
* `second_id`: second identifier (e.g. firm id)
* `time_id`: time identifier 
* `settings`: settings based on `VCHDFESettings`
"""
function find_balanced_set(first_id, second_id, time_id; settings = VCHDFESettings())
    df = DataFrame(first_id = first_id, second_id = second_id, time_id = time_id)

    num_periods = maximum(time_id) - minimum(time_id) + 1
    (settings.print_level > 0) && @info "Number of time periods in the balanced dataset:" num_periods

    df_firmyear_collapsed = @pipe unique(df, [:second_id, :time_id]) |> select(_, [:second_id, :time_id])

    df_firmyear_collapsed = @pipe groupby(df_firmyear_collapsed, :second_id) |> transform(_, nrow => :num_years_present) |> transform(_, :num_years_present => (x -> (x .== num_periods)) => :is_in_balanced)

    new_df = leftjoin(df_firmyear_collapsed, df, on = [:second_id, :time_id])
    (settings.print_level >= 1) && @info "Number of observations belonging to the balanced dataset:" sum(new_df.:is_in_balanced)

    return new_df.:is_in_balanced
end

"""
$(SIGNATURES)

Returns the values for autocorrelation plot and autocorrelation table of the balanced set of the leave out connected set of the input dataset. 

### Arguments
* `y`: outcome vector
* `first_id`: first identifier (e.g. worker id)
* `second_id`: second identifier (e.g. firm id)
* `time_id`: time identifier 
* `settings`: settings based on `VCHDFESettings`
* `autocorr_table`: should compute the autocorrelation table? The default is false
* `lags`: lag vectors, if nothing (default), computes for all lags. 
"""
function leave_out_AR(y, first_id, second_id, time_id, controls = nothing, settings = VCHDFESettings(); autocorr_plot = false, lags = nothing)
    #Create time-second_id identifier (naming it second_id) and rename the second_id to firmid as it is the second id in the time varying AKM model
    data = DataFrame(y = y, first_id = first_id, firmid = second_id, time_id = time_id)

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

    #TODO need to clean out the leverages method
    # @unpack Pii , Mii  , correction_JLA , Bii_first , Bii_second , Bii_cov, Bii_first_bar, Bii_dif_cov, Bii_dif_first_bar, Bii_dif_second_bar = leverages2(settings.leverage_algorithm, X, Fvar, Fvar, Fvar, Fvar, Fvar, settings)


    @unpack θ_first, θ_second, θCOV, β, Dalpha, Fpsi, Pii, Bii_first, Bii_second, Bii_cov, y, X, sigma_i = leave_out_estimation(y,first_id,second_id,controls,settings)

    beta = β
    # xx=X'*X
    # xy=X'*y
    # compute_sol = approxcholSolver(xx;verbose = settings.print_level > 0)
    # beta = compute_sol([xy...];verbose=false)
    # eta=y-X*beta

    # eta_h = eta ./ Mii
    # sigma_i = ( ( y .- mean(y) ) .* eta_h ) .* correction_JLA

    # θ_second = kss_quadratic_form(zeros(size(sigma_i, 1)), Fvar, Fvar, beta, Bii_second)
    # θ_second = kss_quadratic_form(sigma_i, Fvar, Fvar, beta, Bii_second)

    # Creating matrixes
    # we use the term firmyearid to name second_id and year for time_id
    df0 = DataFrame(firmid = firmid, firmyearid = second_id, year = time_id, is_in_balanced = is_in_balanced)
    df0 = @pipe groupby(df0, [:firmyearid]) |> combine(_, :firmid => first => :firmid, :year => first => :year, nrow => :nworkers, :is_in_balanced => maximum => :is_in_balanced)

    first_year = minimum(df0.:year)
    last_year =maximum(df0.:year)

    lags = (lags == nothing) ? (1:(last_year - first_year)) : lags

    acp = nothing
    
    if autocorr_plot
        acp = Array{Union{Missing, Float64}}(missing, last_year - first_year, 1)
        (settings.print_level > 0) && @info "autocorrelation plot values:"
        #TODO check sanity
        for counter in lags
            if (counter + first_year) < last_year
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
                    continue
                end

                @unpack Bii_lag_var, Bii_current_var, Bii_lag_cov = leverages3(X, Fvar, FlagVar, settings)
                
                if size(Bii_current_var, 1) == 0 || isempty(Fvar) || isempty(FlagVar) || isempty(Bii_current_var) 
                    continue
                end
                # θ_lag_cov = kss_quadratic_form(sigma_i[is_in_balanced, :], Fvar[is_in_balanced, :], FlagVar[is_in_balanced, :], beta[is_in_balanced], Bii_lag_cov[is_in_balanced])
                # θ_lag_var = kss_quadratic_form(sigma_i[is_in_balanced, :], FlagVar[is_in_balanced, :], FlagVar[is_in_balanced, :], beta[is_in_balanced], Bii_lag_var[is_in_balanced])
                # θ_lag_var = kss_quadratic_form(sigma_i[is_in_balanced .== 1, :], FlagVar, FlagVar, beta, Bii_lag_var[is_in_balanced .== 1])
                θ_lag_cov = kss_quadratic_form(zeros(size(sigma_i, 1))[is_in_balanced .== 1, :], Fvar, FlagVar, beta, zeros(size(Bii_lag_var, 1))[is_in_balanced .== 1, :])
                θ_current_var = kss_quadratic_form(sigma_i[is_in_balanced .== 1, :], Fvar, Fvar, beta, Bii_current_var[is_in_balanced .== 1])

                if settings.print_level > 0
                    corr_kss_corrected = θ_lag_cov/(sqrt(θ_second) * sqrt(θ_current_var))
                    @info "autocorrelation " counter θ_lag_cov θ_current_var corr_kss_corrected
                    # println(counter, ":")
                    # print("cov: ")
                    # println(θ_lag_cov)
                    # # print("var lag: ")
                    # # println(θ_lag_var)
                    # print("var current: ")
                    # println(θ_current_var)
                    # # print("corr:")
                    # # println(θ_lag_cov/(sqrt(θ_lag_var) * sqrt(θ_current_var)))
                    # print("corr:")
                    
                    # println(corr_kss_corrected)
                    # # println(size(WF, 1))
                    # println(" ")
                end

                acp[counter] = corr_kss_corrected

            end
        end
    end

    acf = Array{Union{Missing, Float64}}(missing, last_year - first_year + 1, last_year - first_year + 1)
    for t in 1:(last_year - first_year + 1)
        acf[t, t] = 1.0::Float64
    end

    for base_year in first_year:last_year
        (settings.print_level > 1) && @info typeof(base_year)
        (settings.print_level > 0) && @info "base year:" base_year
        for counter in lags
            (settings.print_level > 1) && @info typeof(counter)
            if (counter + first_year) < last_year
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
                    continue
                end

                @unpack Bii_lag_var, Bii_current_var, Bii_lag_cov = leverages3(X, Fvar, FlagVar, settings)

                # θ_lag_var = kss_quadratic_form(sigma_i[is_in_balanced .== 1, :], FlagVar, FlagVar, beta, Bii_lag_var[is_in_balanced .== 1])
                # println(Bii_current_var)
                if size(Bii_current_var, 1) == 0 || isempty(Fvar) || isempty(FlagVar) || isempty(Bii_current_var) 
                    continue
                end
                θ_lag_cov = kss_quadratic_form(zeros(size(sigma_i, 1))[is_in_balanced .== 1, :], Fvar, FlagVar, beta, zeros(size(Bii_current_var, 1))[is_in_balanced .== 1, :])
                θ_current_var = kss_quadratic_form(sigma_i[is_in_balanced .== 1, :], Fvar, Fvar, beta, Bii_current_var[is_in_balanced .== 1])
                    
                if settings.print_level > 0
                    cor_kss_corrected = θ_lag_cov/(sqrt(θ_second) * sqrt(θ_current_var))
                    @info "autocorrelations: " counter θ_lag_cov θ_current_var cor_kss_corrected
                    # println(counter, ":")
                    # println("cov: $(θ_lag_cov)")
                    # # print("var lag: ")
                    # # println(θ_lag_var)
                    # println("var current: $(θ_current_var)")
                    # cor_kss_corrected = θ_lag_cov/(sqrt(θ_second) * sqrt(θ_current_var))
                    # println("corr: $(cor_kss_corrected)")
                    # # print("corr2:")
                    # # println(θ_lag_cov/(sqrt(θ_second) * sqrt(θ_current_var)))
                    # # println(size(WF, 1))
                    # println(" ")
                end

                acf[base_year - first_year + 1, base_year - first_year + 1 + counter ] = cor_kss_corrected
            end
        end
    end

    return (obs = obs  , first_id = first_id , second_id = second_id, θ_first = θ_first, θ_second = θ_second, θCOV = θCOV, β = β, Dalpha = Dalpha, Fpsi = Fpsi, Pii = Pii, Bii_first = Bii_first,
    Bii_second = Bii_second, Bii_cov = Bii_cov, y = y , X = X, sigma_i = sigma_i, acf = acf, acp = acp)
end