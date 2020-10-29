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

## Main Functions 


```@docs
    leave_out_estimation(y,first_id,second_id,controls,settings)
    get_leave_one_out_set(y, first_id, second_id, settings, controls)    
```

## Auxiliary Functions

```@docs
    find_connected_set(y, first_idvar, second_idvar, settings)
    prunning_connected_set(yvec, first_idvar, second_idvar, obs_id, settings)
    drop_single_obs(yvec, first_idvar, second_idvar,obs_id)
    compute_movers(first_id,second_id)
    eff_res(::ExactAlgorithm, X,first_id,second_id,match_id, K, settings)
    eff_res(lev::JLAAlgorithm, X,first_id,second_id,match_id, K, settings)
    compute_matchid(second_id,first_id)   
```

# Datatypes in this package

```@docs
    ExactAlgorithm
    JLAAlgorithm
    VCHDFESettings
    
```


