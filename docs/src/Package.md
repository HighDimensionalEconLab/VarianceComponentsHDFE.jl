# Contents

```@contents
Pages = ["Package.md"]
Depth = 3
```

# Setting up Development Enviroment
1. Install [Julia](https://julialang.org/downloads/) and [GitHub Desktop](https://desktop.github.com/) - not strictly required but never hurts to have it!
2. Install `vscode` and follow basic instructions in https://github.com/ubcecon/tutorials/blob/master/vscode.md
   - In particular, https://github.com/ubcecon/tutorials/blob/master/vscode.md#julia, making sure to do the code formatter step.
   - and the git settings in https://github.com/ubcecon/tutorials/blob/master/vscode.md#general-packages-and-setup
3. Clone the repo by either:
   - Clicking on the `Code` then `Open in GitHub Desktop`.
   - Alternatively, you can go `] dev https://github.com/HighDimensionalEconLab/VarianceComponentsHDFE.jl` in a Julia REPL and it will clone it to the `.julia/dev` folder.
   - If you wanted to, you could then drag that folder back into github desktop.
4. Open it in vscode by right-clicking on the folder it installed to, and then opening a vscode project.
5. Open the Julia repl in vscode  (`Ctrl-Shift-P` and then go `Julia REPL` or something to find it.
6. type `] instantiate` to install all of the packages.  Get coffee.
6. In the REPL run `] test` and it should do the full unit test.

# Functions in this package

## Main Function 


```@docs
    leave_out_KSS(y,first_id,second_id;controls, do_lincom , Z_lincom , lincom_labels , settings)       
```

## Auxiliary Functions

```@docs
    find_connected_set(y, first_idvar, second_idvar, settings)
    get_leave_out_set(yvec, first_idvar, second_idvar, obs_id, settings)
    leave_out_estimation(y,first_id,second_id,controls,settings)
    compute_movers(first_id,second_id)
    compute_matchid(second_id,first_id)  
    lincom_KSS(y,X, Z, Transform, sigma_i; lincom_labels ) 
```

# Datatypes in this package

```@docs
    ExactAlgorithm
    JLAAlgorithm
    VCHDFESettings
    
```

# Typical Julia Workflow

```
#Load the required packages
using VarianceComponentsHDFE, DataFrames, CSV

#Load dataset
data = DataFrame!(CSV.File("test.csv"; header=false))

#Extract vectors of outcome, workerid, firmid
id = data[:,1]
firmid = data[:,2]
y = data[:,3]

#You can define the settings using our structures
JL = JLAAlgorithm(num.simulations = 300)
mysettings = VCHDFESettings(leverage_algorithm = JL, first_id_effects=true, cov_effects=true)

#Run KSS with no controls 
θ_first, θ_second, θCOV = leave_out_KSS(y,id,firmid)

#Create some controls and run the routine where we partial out them
controls = indexin(year,unique(sort(year)))
controls = sparse(collect(1:size(y,1)), controls, 1, size(y,1), maximum(controls))
controls = controls[:,1:end-1]

θ_first, θ_second, θCOV = leave_out_KSS(y,id,firmid; controls)

#Perform Lincom Inference using a Region Dummy
data = DataFrame!(CSV.File("lincom.csv"; header=false))
id = data[:,1]
firmid = data[:,2]
y = data[:,5]
region = data[:,4] 
region[findall(region.==-1)].=0

θ_first, θ_second, θCOV = leave_out_KSS(y,id,firmid; do_lincom = true , Z_lincom = region, lincom_labels = ["Region Dummmy"] )


```


