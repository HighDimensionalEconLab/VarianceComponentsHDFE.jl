using LinearAlgebra, SparseArrays,Statistics,DataFrames
#To install Laplacians you need to put in Julia's Prompt "add Laplacians#master"
#For some reason Pkg.add() doesn't work with this.
#Pkg.add("MATLAB")
#using Laplacians
#using Pkg
#include(string(Pkg.dir("Laplacians") , "/src/matlabSolvers.jl"))
using Test
using Parameters
using Random

run_full_test = get(ENV, "JLA_TESTCSV", "false")  == "true" ? true : false
pkg_dir = pkgdir(VarianceComponentsHDFE)


@testset "PruneNetwork" begin
    #Full Network
    y_full = [0.146297; 0.29686 ;  0.54344; 0.432677 ; 0.464866 ; 0.730622 ; 0.619239; 0.753429; 0.0863208; 0.372084 ;  0.958089]
    id_full = [1; 1; 2;2 ; 3; 4;4 ;5;5 ;6;6]
    firmid_full = [1;2;1;1;1;1;2;2;2;2;3]
    obs_id_full = collect(1:11)

    #Pruned Network
    y_p = [0.146297; 0.29686 ;  0.54344; 0.432677 ; 0.464866 ; 0.730622 ; 0.619239; 0.753429; 0.0863208]
    id_p = [1; 1; 2;2 ; 3; 4;4 ;5;5]
    firmid_p = [1;2;1;1;1;1;2;2;2]
    obs_id_p = collect(1:9)

    #Drop-Single Network
    y_lo = [0.146297; 0.29686 ;  0.54344; 0.432677 ; 0.730622 ; 0.619239; 0.753429; 0.0863208]
    id_lo = [1; 1; 2;2 ; 3;3 ;4;4]
    firmid_lo = [1;2;1;1;1;2;2;2]
    obs_id_lo1 = [1;2;3;4;6;7;8;9]
    obs_id_lo = collect(1:8)
    # @test find_connected_set(y_full,id_full,firmid_full; verbose=true) == (obs_id = obs_id_full , y = y_full , id = id_full, firmid = firmid_full )
    @test find_connected_set(y_full,id_full,firmid_full,  VCHDFESettings(print_level = 1)) == (obs_id = obs_id_full , y = y_full , first_id = id_full, second_id = firmid_full )

    @test prunning_connected_set(y_full,id_full,firmid_full, obs_id_full,  VCHDFESettings(print_level = 1)) == (obs_id = obs_id_p , y = y_p , first_id = id_p, second_id = firmid_p )

    @test drop_single_obs(y_p,id_p,firmid_p,obs_id_p) == (obs_id = obs_id_lo1 , y = y_lo , first_id = id_lo, second_id = firmid_lo )

end

@testset "Movers & Matches" begin
    #Drop-Single Network
    y_lo = [0.146297; 0.29686 ;  0.54344; 0.432677 ; 0.730622 ; 0.619239; 0.753429; 0.0863208]
    id_lo = [1; 1; 2;2 ; 3;3 ;4;4]
    firmid_lo = [1;2;1;1;1;2;2;2]
    obs_id_lo1 = [1;2;3;4;6;7;8;9]
    obs_id_lo = collect(1:8)

    @test  compute_matchid(firmid_lo,id_lo) == [1; 2 ; 3;3 ; 4; 5; 6; 6]
    @test  compute_movers(id_lo, firmid_lo) == (movers = Bool[1,1,0,0,1,1,0,0], T= [2, 2, 2, 2, 2, 2, 2, 2])
end

# @testset "Clustering" begin
#     id_lo = [1; 1; 2;2 ; 3;3 ;4;4]
#     obs_id_lo = collect(1:8)

#     #Leave-Out Network Matchid
#     match_id = [1; 2 ; 3;3 ; 4; 5; 6; 6]
# # Xtest = sparse([1;2;3;4;5;6;7;8;1;3;4;5],[1;1;2;2;3;3;4;4;5;5;5;5],[1.0;1.0;1.0;1.0;1.0;1.0;1.0;1.0;-1.0;-1.0;-1.0;-1.0])
# # xxtest = Xtest'*Xtest
# # compute_sol = approxcholSolver(xxtest;verbose=true)
# # settings_exact = VCHDFESettings(leverage_algorithm = ExactAlgorithm(), first_id_effects=true, cov_effects=true)
# # settings_JLA = VCHDFESettings(leverage_algorithm = JLAAlgorithm(num_simulations=3), first_id_effects=true, cov_effects=true)

#     @test check_clustering(id_lo).nnz_2 == 12

#     @test check_clustering(match_id).nnz_2 ==10

#     @test index_constr(obs_id_lo,id_lo, obs_id_lo) == [obs_id_lo  obs_id_lo  obs_id_lo   obs_id_lo]

#     @test index_constr(obs_id_lo,id_lo, match_id) == [obs_id_lo  obs_id_lo  [1;2;3;3;4;5;6;6]   obs_id_lo]

# end

@testset "Solvers" begin
    Random.seed!(1234)
    Xtest = sparse([1;2;3;4;5;6;7;8;1;3;4;5],[1;1;2;2;3;3;4;4;5;5;5;5],[1.0;1.0;1.0;1.0;1.0;1.0;1.0;1.0;-1.0;-1.0;-1.0;-1.0])
    xxtest = Xtest'*Xtest
    compute_sol = approxcholSolver(xxtest;verbose=true)

    rademacher1 = rand(1,8)  .> 0.5
    rademacher1 = rademacher1 - (rademacher1 .== 0)

    @test compute_sol( [Xtest[1,:]...] ;verbose=true) ≈ [ 0.25, -0.5, -0.25, 0.0, -0.5]
    @test compute_sol([rademacher1*Xtest...];verbose=true) ≈  [-0.5, 0.0, -1.5, 0.0, -1.0]
end

@testset "A throughout Test" begin
    y_full = [0.146297; 0.29686 ;  0.54344; 0.432677 ; 0.464866 ; 0.730622 ; 0.619239; 0.753429; 0.0863208; 0.372084 ;  0.958089]
    id_full = [1; 1; 2;2 ; 3; 4;4 ;5;5 ;6;6]
    firmid_full = [1;2;1;1;1;1;2;2;2;2;3]
    obs_id_full = collect(1:11)
    @testset "Exact Algorithm, person_effects = true, cov_effects = true" begin
        settings_exact = VCHDFESettings(leverage_algorithm = ExactAlgorithm(), first_id_effects=true, cov_effects=true, leave_out_level = "obs")

        @unpack obs,  y  , first_id , second_id, controls = get_leave_one_out_set(y_full, id_full, firmid_full, settings_exact, nothing)

    # @unpack θ_first, θ_second, θCOV, obs, β, Dalpha, Fpsi, Pii, Bii_first, Bii_second, Bii_cov = compute_whole(y,first_id,second_id,controls,settings)
        @unpack θ_first, θ_second, θCOV, β, Dalpha, Fpsi, Pii, Bii_first, Bii_second, Bii_cov = leave_out_estimation(y,first_id,second_id,controls,settings_exact)

        # @unpack θ_first, θ_second, θCOV, obs, β, Dalpha, Fpsi, Pii, Bii_first, Bii_second, Bii_cov = compute_whole(y_full,id_full,firmid_full,nothing,settings_exact)

        # println("θ_second: $θ_second \n θ_first:$θ_first \n, θCOV: $θCOV \n obs: $obs \n β: $β \n Dalpha: $Dalpha \n Fpsi: $Fpsi \n Pii: $Pii \n Bii_first: $Bii_first \n Bii_second: $Bii_second \n Bii_cov: $Bii_cov \n")
        @test β ≈ [0.2313735, 0.5076485000000001, 0.6847255, 0.4198749, 0.019589999999999996]
        @test [θ_second, θ_first, θCOV] ≈ [0.004791473893999999, 0.03292163227823429, -0.03087035341429857]
        @test Dalpha ≈  [0.2313735, 0.2313735, 0.5076485000000001, 0.5076485000000001, 0.6847255, 0.6847255, 0.4198749, 0.4198749] 
        @test Fpsi ≈  [0.019589999999999996, 0.0, 0.019589999999999996, 0.019589999999999996, 0.019589999999999996, 0.0, 0.0, 0.0]
        @test Pii ≈ diag(diagm([0.75, 0.75, 0.5, 0.5, 0.75, 0.75, 0.5, 0.5]))
        @test Bii_first ≈ diag(diagm([ -0.25, -0.25, 0.0, 0.0, -0.25, -0.25, 0.0, 0.0]))
        @test Bii_second ≈ diag(diagm([0.5, 0.5, 0, 0, 0.5, 0.5, 0, 0]))
        @test Bii_cov ≈ diag(diagm( [0.625, 0.625, 0.375, 0.375, 0.625, 0.625, 0.375, 0.375]))
    end
    @testset "Exact Algorithm, person_effects = false, cov_effects = false" begin
        settings_exact = VCHDFESettings(leverage_algorithm = ExactAlgorithm(), first_id_effects=false, cov_effects=false, leave_out_level = "obs")

        @unpack obs,  y  , first_id , second_id, controls = get_leave_one_out_set(y_full, id_full, firmid_full, settings_exact, nothing)

    # @unpack θ_first, θ_second, θCOV, obs, β, Dalpha, Fpsi, Pii, Bii_first, Bii_second, Bii_cov = compute_whole(y,first_id,second_id,controls,settings)
        @unpack θ_first, θ_second, θCOV, β, Dalpha, Fpsi, Pii, Bii_first, Bii_second, Bii_cov = leave_out_estimation(y,first_id,second_id,controls,settings_exact)

        # @unpack θ_first, θ_second, θCOV, obs, β, Dalpha, Fpsi, Pii, Bii_first, Bii_second, Bii_cov = compute_whole(y_full,id_full,firmid_full,nothing,settings_exact)

        # println("θ_second: $θ_second \n θ_first:$θ_first \n, θCOV: $θCOV \n obs: $obs \n β: $β \n Dalpha: $Dalpha \n Fpsi: $Fpsi \n Pii: $Pii \n Bii_first: $Bii_first \n Bii_second: $Bii_second \n Bii_cov: $Bii_cov \n")
        @test β ≈ [0.2313735, 0.5076485000000001, 0.6847255, 0.4198749, 0.019589999999999996]
        @test θ_second ≈ 0.004791473893999999
        @test [θ_first, θCOV] == [nothing, nothing]
        @test Dalpha ≈  [0.2313735, 0.2313735, 0.5076485000000001, 0.5076485000000001, 0.6847255, 0.6847255, 0.4198749, 0.4198749] 
        @test Fpsi ≈  [0.019589999999999996, 0.0, 0.019589999999999996, 0.019589999999999996, 0.019589999999999996, 0.0, 0.0, 0.0]
        @test Pii ≈ diag(diagm([0.75, 0.75, 0.5, 0.5, 0.75, 0.75, 0.5, 0.5]))
        @test Bii_first == nothing
        @test Bii_second ≈ diag(diagm([0.5, 0.5, 0, 0, 0.5, 0.5, 0, 0]))
        @test isequal(Bii_cov, nothing)
    end
    @testset "JLA Algorithm with three simulations, person and cov effects = true" begin
        Random.seed!(1234)
        settings_JLA = VCHDFESettings(leverage_algorithm = JLAAlgorithm(num_simulations=3), first_id_effects=true, cov_effects=true, leave_out_level = "obs")

        # @unpack θ_first, θ_second, θCOV, obs, β, Dalpha, Fpsi, Pii, Bii_first, Bii_second, Bii_cov = compute_whole(y_full,id_full,firmid_full,nothing,settings_JLA)

        @unpack obs,  y  , first_id , second_id, controls = get_leave_one_out_set(y_full, id_full, firmid_full, settings_JLA, nothing)

        # @unpack θ_first, θ_second, θCOV, obs, β, Dalpha, Fpsi, Pii, Bii_first, Bii_second, Bii_cov = compute_whole(y,first_id,second_id,controls,settings)
        @unpack θ_first, θ_second, θCOV, β, Dalpha, Fpsi, Pii, Bii_first, Bii_second, Bii_cov = leave_out_estimation(y,first_id,second_id,controls,settings_JLA)
    

        println("θ_second: $θ_second \n θ_first:$θ_first \n, θCOV: $θCOV \n obs: $obs \n β: $β \n Dalpha: $Dalpha \n Fpsi: $Fpsi \n Pii: $Pii \n Bii_first: $Bii_first \n Bii_second: $Bii_second \n Bii_cov: $Bii_cov \n")
        # @test β ≈ [0.2313735, 0.5076485000000001, 0.6847255, 0.4198749, 0.019589999999999996]
        @test isapprox([θ_second, θ_first, θCOV] ,  [0.023848627361476332, 0.0014221007208401193, 0.00017059846351962014], atol=1e-1) #Rewrite this

        @test Pii ≈ diag(diagm([ 0.8749999999999999
            0.7499999999999999
            0.33333333333333337
            0.33333333333333337
            0.875
            0.7499999999999999
            0.33333333333333337
            0.33333333333333337]))
        @test isapprox(Bii_first, diag(diagm([0.12499999999999999, 1.1249999999999998, 0.041666666666666664, 0.041666666666666664, 0.12499999999999999, 0.45833333333333326, 0.7083333333333336, 0.7083333333333336])), atol = 1e-4)	
        @test isapprox(Bii_second, diag(diagm([0.8333333333333335, 0.8333333333333335, 0.0, 0.0, 0.8333333333333335, 0.8333333333333335, 0.0, 0.0])), atol = 1e-4)	
        @test isapprox(Bii_cov, diag(diagm([-0.08333333333333331, -0.24999999999999997, 0.0, 0.0, -0.08333333333333331, -0.24999999999999997, 0.0, 0.0])), atol = 1e-4)
    end
    @testset "JLA Algorithm with three simulations, person and cov effects = false" begin
        Random.seed!(1234)
        settings_JLA = VCHDFESettings(leverage_algorithm = JLAAlgorithm(num_simulations=3), first_id_effects=false, cov_effects=false, leave_out_level = "obs")

        # @unpack θ_first, θ_second, θCOV, obs, β, Dalpha, Fpsi, Pii, Bii_first, Bii_second, Bii_cov = compute_whole(y_full,id_full,firmid_full,nothing,settings_JLA)

        @unpack obs,  y  , first_id , second_id, controls = get_leave_one_out_set(y_full, id_full, firmid_full, settings_JLA, nothing)

        # @unpack θ_first, θ_second, θCOV, obs, β, Dalpha, Fpsi, Pii, Bii_first, Bii_second, Bii_cov = compute_whole(y,first_id,second_id,controls,settings)
        @unpack θ_first, θ_second, θCOV, β, Dalpha, Fpsi, Pii, Bii_first, Bii_second, Bii_cov = leave_out_estimation(y,first_id,second_id,controls,settings_JLA)
    

        println("θ_second: $θ_second \n θ_first:$θ_first \n, θCOV: $θCOV \n obs: $obs \n β: $β \n Dalpha: $Dalpha \n Fpsi: $Fpsi \n Pii: $Pii \n Bii_first: $Bii_first \n Bii_second: $Bii_second \n Bii_cov: $Bii_cov \n")
        # @test β ≈ [0.2313735, 0.5076485000000001, 0.6847255, 0.4198749, 0.019589999999999996]
        @test isapprox(θ_second ,  0.023848627361476332, atol=1e-2) #Rewrite this

        @test Pii ≈ diag(diagm([ 0.8749999999999999
            0.7499999999999999
            0.33333333333333337
            0.33333333333333337
            0.875
            0.7499999999999999
            0.33333333333333337
            0.33333333333333337]))	
        @test isequal(Bii_first, nothing)	
        @test isapprox(Bii_second, diag(diagm([0.8333333333333335, 0.8333333333333335, 0.0, 0.0, 0.8333333333333335, 0.8333333333333335, 0.0, 0.0])), atol = 1e-4)	
        @test isequal(Bii_cov,nothing)

    end
end

if run_full_test
    using CSV
    @testset "JLA Algorithm on Test.csv" begin
        Random.seed!(1234)    
        settings_JLA = VCHDFESettings(leverage_algorithm = JLAAlgorithm(num_simulations=3155), first_id_effects=true, cov_effects=true, leave_out_level = "obs")
        data = CSV.read(joinpath(pkg_dir,"test.csv"); header=false)
        id = data[:,1]
        firmid = data[:,2]
        y = data[:,4]

        @unpack  obs,  y  , first_id , second_id, controls  = get_leave_one_out_set(y,id, firmid,settings_JLA, nothing)

        @unpack θ_first, θ_second, θCOV, β, Dalpha, Fpsi, Pii, Bii_first, Bii_second, Bii_cov = leave_out_estimation(y,first_id,second_id,nothing,settings_JLA)

        @test isapprox([θ_second, θ_first, θCOV],  [0.010505964725588434 , 0.08534310058596604, 0.004550851841563358],atol=1e-2)

    
    end
end

@testset "Autocorrelation" begin
    using CSV
    settings_JLA = VCHDFESettings(leverage_algorithm = JLAAlgorithm(num_simulations=200), first_id_effects=false, cov_effects=false, leave_out_level = "obs", print_level = 0)
    data = CSV.read(joinpath(pkg_dir, "autocorrelation_test2.csv"), DataFrame; header = false)
    id = data[:, 3]
    firmid = data[:, 2]
    year = data[:, 1]
    y = data[:, 4]

    lags = [4]
    @unpack obs, first_id , second_id, θ_first, θ_second, θCOV, β, Dalpha, Fpsi, Pii, Bii_first, Bii_second, Bii_cov, y, X, sigma_i, acf, acp = leave_out_AR(y, id, firmid, year, nothing, settings_JLA, lags = lags)
    @test isapprox(acf[1, 5], 0.42889505736, atol = 1e-4)
end

# @testset "Autocorrelation_small_test" begin
#     using CSV
#     settings_JLA = VCHDFESettings(leverage_algorithm = JLAAlgorithm(num_simulations=200), first_id_effects=false, cov_effects=false, leave_out_level = "obs", print_level = 0)
#     data = CSV.read(joinpath(pkg_dir, "autocorrelation_test2.csv"), DataFrame; header = false)
#     id = data[:, 3]
#     firmid = data[:, 2]
#     year = data[:, 1]
#     y = data[:, 4]

#     lags = [4]
#     @unpack obs, first_id , second_id, θ_first, θ_second, θCOV, β, Dalpha, Fpsi, Pii, Bii_first, Bii_second, Bii_cov, y, X, sigma_i, acf, acp = leave_out_AR(y, id, firmid, year, nothing, settings_JLA, lags = lags)
#     @test isapprox(acf[1, 5], 0.394973, atol = 1e-4)
# end