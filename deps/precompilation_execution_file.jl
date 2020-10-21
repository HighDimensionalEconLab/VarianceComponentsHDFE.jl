using VarianceComponentsHDFE, CSV, DataFrames

# Precompile the compute_whole method

y = [0.146297; 0.29686 ;  0.54344; 0.432677 ; 0.464866 ; 0.730622 ; 0.619239; 0.753429; 0.0863208; 0.372084 ;  0.958089]
id = [1; 1; 2;2 ; 3; 4;4 ;5;5 ;6;6]
firmid = [1;2;1;1;1;1;2;2;2;2;3]
obs_id = collect(1:11)

controls = nothing

settings_exact = VCHDFESettings(leverage_algorithm = ExactAlgorithm(), first_id_effects=true, cov_effects=true)
settings_JLA = VCHDFESettings(leverage_algorithm = JLAAlgorithm(num_simulations=5), first_id_effects=true, cov_effects=true)

compute_whole(y,id,firmid,controls,settings_exact)
compute_whole(y,id,firmid,controls,settings_JLA)

# Precompile reading and writing CSV files

output = DataFrame(y = y, id = id)
CSV.write("csv_file.csv",output)
data = DataFrame!(CSV.File("csv_file.csv"))
