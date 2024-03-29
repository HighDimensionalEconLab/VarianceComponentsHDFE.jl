module VarianceComponentsHDFE

# using DataDeps
using CSV
using DataFrames, Parameters
using LinearAlgebra, SparseArrays, Random, Statistics
using SparseArrays, LightGraphs
using Distributions
using LinearOperators, FastClosures, Krylov
using ArgParse
using FixedEffects

#include("init.jl")
include("leave_out_correction.jl")
include("parameters_settings.jl")
include("laplacians/Laplacians.jl")

using .Laplacians
include("solvers.jl")

export find_connected_set,prunning_connected_set,drop_single_obs
export compute_movers, accumarray
export lincom_KSS, compute_matchid, leave_out_estimation, get_leave_one_out_set
export VCHDFESettings, JLAAlgorithm, ExactAlgorithm, AbstractLeverageAlgorithm
export leave_out_KSS, leverages, kss_quadratic_form, sigma_for_stayers

# Exporting these functions for ease of benchmarking/testing
export computeLDLinv, approxcholOperator, approxcholSolver

# Executable parse 
function parse_commandline()
    argparsesettings_obj = ArgParseSettings()

    @add_arg_table! argparsesettings_obj begin
        "path"
            help = "path to CSV file containing data"
            required = true
        "--first_id"
            help = "column index in CSV file for the first ID (e.g. Person).  Use the most granular type."
            arg_type = Int
            default = 1
        "--second_id"
            help = "column index in CSV file for the second ID (e.g. Firm).  Use the less granular type."
            arg_type = Int
            default = 2
        "--outcome_id"
            help = "column index in CSV file for outcome (e.g. Wage)."
            arg_type = Int
            default = 4
        "--no_first_id_effects"
            help = "No computing and showing of first_id effects"
            action = :store_true
        "--no_cov_effects"
            help = "No computing and showing of covariace effects"
            action = :store_true
        "--leave_out_level"
            help = "leave out level: obs or match. It defaults to match."
            arg_type = String
            default = "match"
        "--algorithm"
            help = "type of algorithm: Exact or JLA. It defaults to be Exact if the number of observations is less than 5000, and JLA otherwise."
            arg_type = String
            default = "Default"
        "--simulations"
            help = "number of simulations in the JLA algorithm. If 0, defaults to 200"
            arg_type = Int
            default = 0
        "--header"
            help = "CSV file contains header"
            action = :store_true
        "--first_id_display"
            help = "The display text associated with first_id (e.g. Person)."
            arg_type = String
            default = "Person"
        "--second_id_display"
            help = "The display text associated with second_id (e.g. Firm)"
            arg_type = String
            default = "Firm"
        "--outcome_id_display"
            help = "The display text associated with outcome_id (e.g. Wage)"
            arg_type = String
            default = "Wage"
        "--detailed_csv_path"
            help = "path to the CSV for the detailed output for each observable"
            arg_type = String
            default = joinpath(pwd(), "variance_components.csv")
        "--results_path"
            help = "path to the results of the output"
            arg_type = String
            default = joinpath(pwd(), "results.txt")
        "--write_detailed_csv"
            help = "write the detailed output to a CSV"
            action = :store_true
        "--write_results"
            help = "write the results to a file"
            action = :store_true
        #adding the print level command
        "--print_level"
            help = "Level of verbosity of output."
            arg_type = Int
            default = 1
        "--do_lincom"
            help = "perform linear combination based on observables Z"
            action = :store_true
        #Select the number of columns or put labels
        "--covariates"
            help = "Covariates. You can set column numbers or String name."
            nargs = '*'  
            default = []
        #Select the number of columns or put labels
        "--lincom_covariates" 
            help = "Lincom covariates. You can set column numbers or String name."
            nargs = '*'  
            default = []
    end

    return parse_args(argparsesettings_obj)
end

function julia_main()::Cint
    try
        real_main()
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end

function real_main()
    parsed_args = parse_commandline()

    println("Present working directory is $(pwd())")

    path = parsed_args["path"]
    #TODO: not sure about the paths
    if isabspath(path) == false
        path = joinpath(pwd(), path)
    end
    header = parsed_args["header"]
    first_idx = parsed_args["first_id"]
    second_idx = parsed_args["second_id"]
    outcome_idx = parsed_args["outcome_id"]
    leave_out_level = parsed_args["leave_out_level"]
    algorithm = parsed_args["algorithm"]
    first_id_effects = parsed_args["no_first_id_effects"] == 1 ? 0 : 1
    cov_effects = parsed_args["no_cov_effects"] == 1 ? 0 : 1
    simulations = parsed_args["simulations"]
    first_id_display = parsed_args["first_id_display"]
    second_id_display = parsed_args["second_id_display"]
    outcome_id_display = parsed_args["outcome_id_display"]
    print_level = parsed_args["print_level"]
    covariates = parsed_args["covariates"]
    lincom_covariates = parsed_args["lincom_covariates"]

    algorithm = lowercase(algorithm)
    leave_out_level = lowercase(leave_out_level)
    first_id_display_small = lowercase(first_id_display)
    second_id_display_small = lowercase(second_id_display)

    print_level > 0 && println("Number of threads: $(Threads.nthreads())")

    data  = DataFrame(CSV.File(path; header=header))
    first_id = data[:,first_idx]
    second_id = data[:,second_idx]
    y = data[:,outcome_idx]

    if covariates == []
        controls = nothing
    else
        #It creates a dictionary of the covariates 
        if typeof(covariates[1]) == String 
            #Build controls matrix 
            controls = data[:,covariates[1]]
            if length(covariates)>=2
                for k=2:length(covariates)
                    hcat(controls, data[:,covariates[k]])
                end
            end    
        else
            println("WARNING: Elements of covariates are not defined correctly. No partialling out will be performed.")
            controls = nothing 
        end
    end

    if algorithm == "default"
        if length(y) > 5000
            algorithm = "jla"
        else
            algorithm = "exact"
        end
    end
    

    if algorithm == "exact"
        #todo work with settings arguments
        settings = VCHDFESettings(leverage_algorithm =  ExactAlgorithm(),
            first_id_effects =first_id_effects,
            cov_effects = cov_effects,
            first_id_display = first_id_display,
            first_id_display_small = lowercase(first_id_display),
            second_id_display = second_id_display,
            second_id_display_small = lowercase(second_id_display),
            outcome_id_display= outcome_id_display,
            outcome_id_display_small = lowercase(outcome_id_display),
            leave_out_level = leave_out_level,
            print_level = print_level
        )
    else
        settings = VCHDFESettings(leverage_algorithm = JLAAlgorithm(num_simulations=simulations),
            first_id_effects =first_id_effects,
            cov_effects = cov_effects,
            first_id_display = first_id_display,
            first_id_display_small = lowercase(first_id_display),
            second_id_display = second_id_display,
            second_id_display_small = lowercase(second_id_display),
            outcome_id_display= outcome_id_display,
            outcome_id_display_small = lowercase(outcome_id_display),
            leave_out_level = leave_out_level,
            print_level = print_level
        )
    end

    println("Running KSS correction with the following options")
    println("Leave Out Strategy : Leave $(settings.leave_out_level) out")
    if algorithm == "jla"
        algo = "JLA"
        simul = settings.leverage_algorithm.num_simulations == 0 ? 200 : settings.leverage_algorithm.num_simulations
        println("Algorithm for Computation of Statistical Leverages : $(algo) with $(simul) simulations")
    else
        algo = "Exact"
        println("Algorithm for Computation of Statistical Leverages : $(algo) \n")
    end

    first_id_old = first_id 
    second_id_old = second_id
    y_untransformed = y

    @unpack obs,  y  , first_id , second_id, controls = get_leave_one_out_set(y, first_id, second_id, settings, controls)

    if settings.print_level>1
        num_movers = length(unique(compute_movers(first_id,second_id).movers .* first_id)) - 1 

        summary = """
        Information on Leave Out Connected Sample \n
        Number of observations: $(length(obs)) 
        Number of $(settings.first_id_display_small)s: $(maximum(first_id)) 
        Number of $(settings.second_id_display_small)s: $(maximum(second_id)) 
        Number of Movers : $(num_movers)
        Mean of $(settings.outcome_id_display): $(mean(y)) 
        Variance of $(settings.outcome_id_display): $(var(y))
        """ 
    println(summary)
    end

    Z_lincom = nothing 
    if  parsed_args["do_lincom"]
      #Construct Z_lincom 
      if lincom_covariates == [] 
        println("\n User asked for lincom but no covariates were specified. This step will not be performed.")
      else
        if typeof(lincom_covariates[1]) == String 
            #Build lincom controls matrix 
            lincom_labels = []
            push!(lincom_labels,String.(lincom_covariates[1]))
            
            Z_lincom = data[obs,lincom_covariates[1]]
            if length(lincom_covariates)>=2
                for k=2:length(lincom_covariates)
                    hcat(Z_lincom, data[obs,lincom_covariates[k]])
                    push!(lincom_labels,String.(lincom_covariates[k]))
                end
            end     
        else
            println("WARNING: Elements of lincom covariates are not defined correctly. No inference will be performed.")
        end
    end 
    end

    #Residualize outcome variable 
    if controls != nothing  
        println("\nPartialling out controls from $(settings.outcome_id_display)...")

        #We'll use FixedEffects.jl to partial-out FEs and then obtain coefs for controls 
        yres = deepcopy(y)
        Xres = Float64.(Matrix(controls))

        # partial out fixed effects
        fe1 = FixedEffect(first_id)
        fe2 = FixedEffect(second_id)
        cols = eachcol(Xres)
        sumsquares_pre = [sum(abs2, x) for x in cols]

        solve_residuals!(yres, [fe1, fe2])
        for col in cols
            solve_residuals!(col, [fe1, fe2])
        end 

        # re-compute 2-norm (sum of squares) for each variable
        sumsquares_post = [sum(abs2, x) for x in cols]

        # mark variables that are likely to be collinear with the fixed effects (e.g. an intercept)
        collinear_tol = min(1e-6, 1e-8 / 10) #Don't want to add too many tuning parameters, might change later
        collinear_fe = (sumsquares_post ./ sumsquares_pre) .< collinear_tol 
        # set variables that are likely to be collinear with the fixed effects to zero
        for i in 1:size(Xres,2)
            if collinear_fe[i] 
                    @info "Control variable on position $(i) is collinear with the fixed effects. Code will exclude it."
                    # set to zero so that removed when taking basis
                    Xres[:,i] .= 0.0
            end
        end

        XXres = Xres'Xres
        invXXres  =invsym!(Symmetric(deepcopy(XXres)), setzeros=true)
        basis = diag(invXXres).<0

        if !all(basis)
            Xres = Xres[:, basis]
            XXres = XXres[basis, basis]
        end
        Xhat = Xres
        XhatXhat = Symmetric(XXres)

        Xyres = Symmetric(hvcat(2, XhatXhat, Xhat'reshape(yres, length(yres), 1), 
                            zeros(size(Xhat, 2))', [0.0]))
        invsym!(Xyres; diagonal = 1:size(Xhat, 2))
        coef = Xyres[1:(end-1),end] #Row vector 

        y = y - Matrix(controls)[:,basis]*coef

        # NOW WE CAN GO BACK TO NO CONTROLS CASE
        controls = nothing
        println("Partialling out completed.")
    end

    # @unpack θ_first, θ_second, θCOV, obs, β, Dalpha, Fpsi, Pii, Bii_first, Bii_second, Bii_cov = compute_whole(y,first_id,second_id,controls,settings)
    @unpack θ_first, θ_second, θCOV, β, Dalpha, Fpsi, Pii, Bii_first, Bii_second, Bii_cov, y, X, sigma_i = leave_out_estimation(y,first_id,second_id,controls,settings)

    if  parsed_args["do_lincom"]   
        #Construct Transform 
        F = sparse(collect(1:length(second_id)),second_id,1)
        J = size(F,2)
        S = sparse(1.0I, J-1, J-1)
        S = vcat(S,sparse(-zeros(1,J-1)))
        Transform = hcat(spzeros(length(second_id),maximum(first_id)), -F*S )

        @unpack test_statistic, linear_combination , SE_linear_combination_KSS, SE_naive = lincom_KSS(y,X, Z_lincom, Transform, sigma_i; lincom_labels)
    end 

    if parsed_args["write_detailed_csv"]

        y_output = y
        first_id_output = first_id
        second_id_output = second_id

        max_length = length(obs)

        if leave_out_level == "match"
            #We need to take the output to the person-year space 
            match_id = compute_matchid(second_id, first_id)
            weights = accumarray(match_id, 1)

            Pii = [Pii...]
            Pii = vcat(fill.(Pii, weights)...)
            Bii_first = Bii_first == nothing ? nothing : vcat(fill.(Bii_first, weights)...)
            Bii_second = Bii_second == nothing ? nothing : vcat(fill.(Bii_second, weights)...)
            Bii_cov = Bii_cov == nothing ? nothing : vcat(fill.(Bii_cov, weights)...)
            Dalpha = vcat(fill.(Dalpha, weights)...)
            Fpsi = vcat(fill.(Fpsi, weights)...)
        end


        #Make sure that we can put them in a DataFrame 
        obs = [obs...]
        first_id_old = [first_id_old[obs]...]
        first_id_output = [first_id_output...]
        second_id_old = [second_id_old[obs]...]
        second_id_output = [second_id_output...]
        y_output = [y_untransformed[obs]...]
        Dalpha = [Dalpha...] 
        #Fpsi has the opposite sign, so we fix it here.
        Fpsi = [-1.0*Fpsi...]
        Pii = [Pii...]

        #todo rename the DataFrame arguments
        output = DataFrame(observation = obs , first_id_old = first_id_old, first_id = first_id_output ,
                           second_id_old = second_id_old, second_id = second_id_output, y = y_output,
                           Dalpha = Dalpha , Fpsi = Fpsi, Pii = Pii, 
                           Bii_first = Bii_first === nothing ? missings(max_length) : Bii_first,
                           Bii_second = Bii_second === nothing ? missings(max_length) :  Bii_second,
                           Bii_cov = Bii_cov === nothing ? missings(max_length) : Bii_cov
                        )
        output_path = parsed_args["detailed_csv_path"]
        #TODO not sure about the paths
        if isabspath(output_path) == false
            output_path = joinpath(pwd(), output_path)
        end
        CSV.write(output_path,output)
    end

    if parsed_args["write_results"]
        y_py = y_untransformed[obs]
        var_den = var(y_py)
        movers , T = compute_movers(first_id, second_id)
        num_movers = length(unique(movers .* first_id)) - 1 

        output_template = """
            Number of observations (Leave Out Sample): $(length(obs)) \n 
            Number of $(first_id_display_small)s: $(maximum(first_id)) \n 
            Number of $(second_id_display_small)s: $(maximum(second_id)) \n
            Number of movers: $(num_movers) \n
            Mean Outcome: $(mean(y_py)) \n
            Variance Outcome: $(var(y_py)) \n
            Bias Corrected Variance of $(second_id_display) Effects: $θ_second \n
        """ 

        if parsed_args["write_results"]
            try 
                output_path = parsed_args["results_path"]
                if isabspath(output_path) == false
                    output_path = joinpath(pwd(), output_path)
                end
                open(output_path, "w") do io
                   write(io, output_template)
                   if first_id_effects == 1 
                        write(io, "    Bias Corrected Variance of  $(first_id_display) Effects:  $θ_first \n\n")
                   end

                   if cov_effects == 1 
                        write(io, "    Bias Corrected Covariance of  $(first_id_display)-$(second_id_display) Effects: $θCOV \n\n")
                   end

                   if (first_id_effects == 1 ) && (cov_effects == 1 )
                        θ_first = abs(θ_first) #This is only required for the very small dataset that does the precompilation
                        θ_second = abs(θ_second) #This is only required for the very small dataset that does the precompilation
                        corr = θCOV/(sqrt(θ_first)*sqrt(θ_second)) 
                        r2 = (θ_second+2*θCOV+θ_first)/var_den
                        write(io, "    Bias Corrected Correlation of  $(first_id_display)-$(second_id_display) Effects: $corr \n\n")                   
                        write(io, "    Bias Corrected Fraction of Variance explained by  $(first_id_display)-$(second_id_display) Effects: $r2 \n\n")                   
                   end
                   

                    if Z_lincom != nothing 
                        #Establish basis 
                        zz = Z_lincom'*Z_lincom
                        invzz  =invsym!(Symmetric(deepcopy(zz)), setzeros=true)
                        basis = diag(invzz).<0

                        # Modify labels if they are present 
                        lincom_labels = lincom_labels == nothing ? nothing : lincom_labels[basis[2:end]]

                        #Check if the basis includes the intercept
                        has_cons = basis[1] == true ? 1 : 0   

                        #Write the output of inference 
                        r = size(Z_lincom_col,2)
                        write(io,"    Results of High Dimensional Lincom \n\n")
                        if lincom_labels == nothing 
                            start = has_cons == 1 ? 2 : 1
                            for q=start:r
                                if q <= r
                                    ncol = has_cons ==  1 ? q-1 : q
                                    output_inference = """
                                        Coefficient of Column $(ncol): $(linear_combination[q]) \n
                                        Traditional HC Standard Error of Column $(ncol): $(SE_naive[q]) \n 
                                        KSS Standard Error of Column $(ncol): $(SE_linear_combination_KSS[q]) \n 
                                        T-Statistic of Column $(ncol): $(test_statistic[q]) \n
                                    """
                                    write(io,output_inference)
                                end
                            end
                        else
                            start = has_cons == 1 ? 2 : 1
                            for q=start:r
                                pos = has_cons ==  1 ? q-1 : q
                                tell_me = lincom_labels[pos]
                                output_inference = """
                                    Coefficient on $(tell_me): $(linear_combination[q-1]) \n
                                    Traditional HC Standard Error on $(tell_me): $(SE_naive[q-1]) \n 
                                    KSS Standard Error on $(tell_me): $(SE_linear_combination_KSS[q-1]) \n 
                                    T-Statistic on $(tell_me): $(test_statistic[q-1]) \n
                                """
                                write(io,output_inference)  
                            end
                        end
                    end
                 end
            catch e
                println("Error writing output to $(output_path)")
                display(e)
           end
        end
    end
    
end

end

