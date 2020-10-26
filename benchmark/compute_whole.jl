using VarianceComponentsHDFE
using Random, BenchmarkTools, Test
using CSV
using ArtifactUtils, Pkg.Artifacts
# Make sure LDLFactorizations is version 0.5.0 and that the `multiple-rhs` branch is checked out

pkg_dir = pkgdir(VarianceComponentsHDFE)
artifacts_path = joinpath(pkgdir(VarianceComponentsHDFE), "benchmark\\Artifacts.toml")


#Set Random seed
Random.seed!(1234)


# 1) Medium sized network benchmarks 
data = CSV.read(joinpath(pkg_dir, "test.csv");header=false)

y = [data[:,4]...]
id = [data[:,1]...]
firmid = [data[:,2]...]
obs_id = collect(1:length(y))

settings_JLA = VCHDFESettings(leverage_algorithm = JLAAlgorithm(num_simulations=3155), first_id_effects=true, cov_effects=true)
@btime compute_whole(y,id,firmid,nothing,settings_JLA) #55.669s Julia vs 152.97s Matlab (32GB Ram and 4 Threads)


# 2) Large sized network benchmarks 
data_hash = artifact_hash("large_nocontrols_pruned",artifacts_path)
data = CSV.read(joinpath(artifact_path(data_hash), "large_nocontrols_pruned.csv"))


y = [data[:,4]...]
id = [data[:,1]...]
firmid = [data[:,2]...]
obs_id = collect(1:length(y))

settings_JLA = VCHDFESettings(leverage_algorithm = JLAAlgorithm(num_simulations=100), first_id_effects=true, cov_effects=true)

@btime compute_whole($y,$id,$firmid,nothing,settings_JLA) #978.056s Julia vs 606.07s Matlab (32GB Ram and 4 Threads, 100 simulations)

settings_JLA = VCHDFESettings(leverage_algorithm = JLAAlgorithm(num_simulations=200), first_id_effects=true, cov_effects=true)

@btime compute_whole($y,$id,$firmid,nothing,settings_JLA) #1702.217 ss Julia vs 1115.04s Matlab (32GB Ram and 4 Threads, 200 simulations)


settings_JLA = VCHDFESettings(leverage_algorithm = JLAAlgorithm(num_simulations=500), first_id_effects=true, cov_effects=true)

@btime compute_whole(y,id,firmid,nothing,settings_JLA) # 4280.219173 s in Julia vs 3328.10s in Matlab (32GB Ram and 4 Threads, 500 simulations)


#Simulations on cluster for 250 rademachers: 949.138618 seconds Julia vs 379.719083  (256GB and 32cores )
#Simulations on cluster for 500 rademachers: 1463.095847 seconds Julia vs 628.135439 seconds Matlab (256GB and 32cores)
#Simulations on cluster for heuristic - 1441 rademachers : 3392.993448 seconds Julia vs 1651.819884 seconds Matlab  (256GB and 32cores)




# 3) Huge sized network benchmarks 
data_hash = artifact_hash("huge_nocontrols_pruned",artifacts_path)
data = CSV.read(joinpath(artifact_path(data_hash), "huge_pruned_main.csv"))

y = [data[:,4]...]
id = [data[:,1]...]
firmid = [data[:,2]...]
obs_id = collect(1:length(y))

settings_JLA = VCHDFESettings(leverage_algorithm = JLAAlgorithm(num_simulations=100), first_id_effects=true, cov_effects=true)




