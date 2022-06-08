#Load the required packages
using VarianceComponentsHDFE, DataFrames, CSV
using SparseArrays
using Parameters

#Load dataset
data = DataFrame(CSV.File("test.csv"; header=false))

#Extract vectors of outcome, workerid, firmid
id = data[:,1]
firmid = data[:,2]
y = data[:,6]
year = data[:, 5]

#You can define the settings using our structures
JL = JLAAlgorithm(num_simulations = 300)
mysettings = VCHDFESettings(leverage_algorithm = JL, first_id_effects=true, cov_effects=true)

#Run KSS with no controls 
θ_first, θ_second, θCOV = leave_out_KSS(y,id,firmid)

#Create some controls and run the routine where we partial out them
controls = indexin(year,unique(sort(year)))
controls = sparse(collect(1:size(y,1)), controls, 1, size(y,1), maximum(controls))
controls = controls[:,1:end-1]

θ_first, θ_second, θCOV = leave_out_KSS(y,id,firmid; controls)

@unpack obs,  y  , first_id , second_id, controls = get_leave_one_out_set(y, id, firmid, mysettings, controls)
# X = leave_out_estimation(y, id, firmid, controls, mysettings)
@unpack θ_first, θ_second, θCOV, β, Dalpha, Fpsi, Pii, Bii_first, Bii_second, Bii_cov, y, X, sigma_i = leave_out_estimation(y,first_id, second_id,controls,mysettings)
#Perform Lincom Inference using a Region Dummy
data = DataFrame!(CSV.File("lincom.csv"; header=false))
id = data[:,1]
firmid = data[:,2]
y = data[:,5]
region = data[:,4] 
region[findall(region.==-1)].=0

θ_first, θ_second, θCOV = leave_out_KSS(y,id,firmid; do_lincom = true , Z_lincom = region, lincom_labels = ["Region Dummmy"] )
