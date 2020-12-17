module VarianceComponentsHDFE

# using DataDeps
using CSV
using DataFrames, DataFramesMeta, Parameters
using LinearAlgebra, SparseArrays, Random, Statistics
using SparseArrays, LightGraphs, DataFramesMeta
using Distributions
using LinearOperators, FastClosures, Krylov
using ArgParse

#include("init.jl")
include("leave_out_correction.jl")
include("parameters_settings.jl")
include("laplacians/Laplacians.jl")

using .Laplacians
include("solvers.jl")

export find_connected_set,prunning_connected_set,drop_single_obs, index_constr
export compute_movers, check_clustering, eff_res
export do_Pii, lincom_KSS, compute_matchid, leave_out_estimation, get_leave_one_out_set
export VCHDFESettings, JLAAlgorithm, ExactAlgorithm, AbstractLeverageAlgorithm

# Exporting these functions for ease of benchmarking/testing
export computeLDLinv, approxcholOperator, approxcholSolver

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
        "--observation_id"
            help = "column index in CSV file for observation (e.g. Wage)."
            arg_type = Int
            default = 4
        "--no_first_id_effects"
            help = "No computing and showing of first_id effects"
            action = :store_false
        "--no_cov_effects"
            help = "No computing and showing of covariace effects"
            action = :store_false
        "--algorithm"
            help = "type of algorithm: Exact or JLA. It defaults to be Exact if the number of observations is less than 5000, and JLA otherwise."
            arg_type = String
            default = "Default"
        "--simulations"
            help = "number of simulations in the JLA algorithm. If 0, defaults to 100 * log(#total fixed effect)"
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
        "--observation_id_display"
            help = "The display text associated with observable_id (e.g. Wage)"
            arg_type = String
            default = "Wage"
        # "--write_CSV"
        #     help = "write output to a CSV"
        #     action = :store_true
        # "--output_path"
        #     help = "path to output CSV"
        #     arg_type = String
        #     default = "VarianceComponents.csv"
        
        # in support of fixing up the output:
        "--detailed_output_path"
            help = "path to the CSV for the detailed output for each observable"
            arg_type = String
            default = joinpath(pwd(), "variance_components.csv")
        "--results_path"
            help = "path to the results of the output"
            arg_type = String
            default = joinpath(pwd(), "results.txt")
        "--write_detailed_CSV"
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

    println("present working directory is $(pwd())")

    path = parsed_args["path"]
    #TODO: not sure about the paths
    if isabspath(path) == false
        path = joinpath(pwd(), path)
    end
    header = parsed_args["header"]
    first_idx = parsed_args["first_id"]
    second_idx = parsed_args["second_id"]
    observation_idx = parsed_args["observation_id"]
    algorithm = parsed_args["algorithm"]
    first_id_effects = parsed_args["no_first_id_effects"] == 1 ? 0 : 1
    cov_effects = parsed_args["no_cov_effects"] == 1 ? 0 : 1
    simulations = parsed_args["simulations"]
    first_id_display = parsed_args["first_id_display"]
    second_id_display = parsed_args["second_id_display"]
    observation_id_display = parsed_args["observation_id_display"]
    print_level = parsed_args["print_level"]

    first_id_display_small = lowercase(first_id_display)
    second_id_display_small = lowercase(second_id_display)

    print_level > 0 && println("Number of threads: $(Threads.nthreads())")

    data  = DataFrame!(CSV.File(path; header=header))
    first_id = data[:,first_idx]
    second_id = data[:,second_idx]
    #todo rename y to observations
    y = data[:,observation_idx]

    controls = nothing

    if algorithm == "Default"
        if length(y) > 5000
            algorithm = "JLA"
        else
            algorithm = "Exact"
        end
    end
    

    if algorithm == "Exact"
        #todo work with settings arguments
        settings = VCHDFESettings(leverage_algorithm =  ExactAlgorithm(),
            first_id_effects =first_id_effects,
            cov_effects = cov_effects,
            first_id_display = first_id_display,
            first_id_display_small = lowercase(first_id_display),
            second_id_display = second_id_display,
            second_id_display_small = lowercase(second_id_display),
            observation_id_display= observation_id_display,
            observation_id_display_small = lowercase(observation_id_display),
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
            observation_id_display= observation_id_display,
            observation_id_display_small = lowercase(observation_id_display),
            print_level = print_level
        )
    end

    # @assert settings.first_id_effects == true && settings.cov_effects == true

    @unpack obs,  y  , first_id , second_id, controls = get_leave_one_out_set(y, first_id, second_id, settings, nothing)

    # @unpack θ_first, θ_second, θCOV, obs, β, Dalpha, Fpsi, Pii, Bii_first, Bii_second, Bii_cov = compute_whole(y,first_id,second_id,controls,settings)
    @unpack θ_first, θ_second, θCOV, β, Dalpha, Fpsi, Pii, Bii_first, Bii_second, Bii_cov = leave_out_estimation(y,first_id,second_id,controls,settings)

    if parsed_args["write_detailed_CSV"]

        y_output = y
        first_id_output = first_id
        second_id_output = second_id

        max_length = length(obs)

        #todo rename the DataFrame arguments
        output = DataFrame(observation = replace(obs, nothing => missing),
                           first_id = replace(first_id_output, nothing => missing),
                           second_id = replace(second_id_output, nothing => missing),
                           y = replace(y_output, nothing => missing),
                        #    beta = vcat(β,missings(max(max_length-length(β),0))),
                           D_alpha = replace(Dalpha, nothing => missing),
                           F_psi = replace(Fpsi, nothing => missing),
                           Pii = replace(Pii, nothing => missing),
                           Bii_first = replace(Bii_first, nothing => missing),
                           Bii_second = replace(Bii_second, nothing => missing),
                           Bii_cov = replace(Bii_cov, nothing => missing),
                        #    variance_comp_second_effects = [θ_second; missings(max_length-1)],
                        #    variance_comp_first_effects = [θ_first; missings(max_length-1)],
                        #    covariance_comp_effects = [θCOV; missings(max_length-1)]
                        )
        output_path = parsed_args["detailed_output_path"]
        #TODO not sure about the paths
        if isabspath(output_path) == false
            output_path = joinpath(pwd(), output_path)
        end
        CSV.write(output_path,output)
    end

    if parsed_args["write_results"]
        output_template = """
            Number of observations (Leave-out Sample): $(length(obs)) \n 
            Number of $(first_id_display_small)s: $(length(first_id)) \n 
            Number of $(second_id_display_small)s: $(length(second_id)) \n
            Mean Outcome: $(mean(y)) \n
            Variance Outcome: $(var(y)) \n
            Max Pii: $(maximum(Pii)) \n
            Variance of $(second_id_display) Effects: $θ_second \n
            Covariance of  $(first_id_display)-$(second_id_display) Effects: $θCOV \n
            Variance of  $(first_id_display) Effects:  $θ_first 
        """ 

        if parsed_args["write_results"]
            try 
                output_path = parsed_args["results_path"]
                if isabspath(output_path) == false
                    output_path = joinpath(pwd(), output_path)
                end
                open(output_path, "w") do io
                   write(io, output_template)
                 end
            catch e
                println("Error writing output to $(output_path)")
                display(e)
           end
        end
    end
    
end

end
