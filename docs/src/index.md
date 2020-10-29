

# VarianceComponentsHDFE

This package estimates a two-way fixed effects model and computes its associated variance components using the methodology developed by Kline, Saggio and Sølvsten (KSS). This is achieved by running an executable (app) in the terminal/powershell. The user needs to input to the app the path to the original data (in .csv format) and indicate the corresponding column in the .csv that contains the first identier (e.g. the worker id), the second identier (e.g the firm id) and the outcome (e.g. log wage). The link to the repository is [this](https://github.com/HighDimensionalEconLab/VarianceComponentsHDFE.jl).


## About the executable/app

The algorithm prints the plug-in and the bias-corrected variance components estimators for the first identifier effects (e.g. variance of worker effects), the second identifier effects (e.g. variance of firm effects), and the covariance of both these effects (e.g. covariance of worker-firm effects). The user may choose to compute only a subset of these three components. Additionally, the executable will create a CSV file that stores the outcome variable, the first and second set of identifiers, and the estimated fixed effects for every observation that belongs to the leave-out connected set as defined in KSS. The algorithm also saves the statistical leverages, Pii, as well the weighting terms defined as Bii in KSS for a given variance component.

## Installation (Windows)

1. Open up a powershell terminal (Windows + R, then type "powershell"), and press Enter.

Note: It is recommended to run powershell as administrator for installation. To do this, open Windows menu, type "powershell". Right-click on the powershell, click "run as administrator"

2. Change the current directory by typing `cd "your_desired_path"`. 

 Hint : To copy-paste into the terminal use the standard Ctrl+C and paste into the powershell by using right click.

3. In the powershell, install the latest version:

```
wget https://github.com/HighDimensionalEconLab/VarianceComponentsHDFE.jl/releases/download/v0.1.5.1/vchdfe-v0.1.5.1-ubuntu-latest.tar.gz -O vchdfe-windows-latest.tar.gz
 
tar -xvf vchdfe-windows-latest.tar.gz
```

 Note that to be able to use wget on Windows Internet Explorer should have been launched at least once. 

4. (OPTIONAL): Add the current directory to PATH:

```
setx PATH "$env:path;$pwd\vchdfe\bin" -m
```
Note: To permanently change the path, you need to run powershell as administrator. 

5. Test using the sample test file:

```
vchdfe vchdfe\bin\test.csv
```

## Executable guide

Before this, make sure you've completed the installation steps shown above. The basic syntax for this command is 

```
vchdfe path_to_data [--option option_value]
```
where `[]` denote optional arguments. 

Note: You can use this simple syntax if you added the binary path to your system path (Step 4 in above). Otherwise, you need to point to the binary directory to call vchdfe. 

To use the executable you only need to open powershell (Windows + R, then type "powershell" and press Enter). You can see a list of examples provided in a section below or type in the powershell

```
vchdfe --help
```

to see a complete list of available arguments:
  
    positional arguments:
      path                  path to CSV file containing data

    optional arguments:
      --first_id FIRST_ID   column index in CSV file for the first ID
                            The user should specify the most "granular"
                            id here. For instance, in the AKM context,
                            this id would correspond to worker ids.
                            (type: Int64, default: 1)
      --second_id SECOND_ID
                            column index in CSV file for the second ID
                            Use the less granular type (e.g. Firm in the 
                            AKM example)
                            (type: Int64, default: 2)
      --observation_id OBSERVATION_ID
                            column index in CSV file for observation (e.g.
                            Wage). (type: Int64, default: 4)
      --no_first_id_effects No computing and showing of first_id effects
      --no_cov_effects      No computing and showing of covariance effects
      --algorithm ALGORITHM
                            type of algorithm: Exact or JLA. It defaults
                            to be Exact if the number of observations is
                            less than 5000, and JLA otherwise. (default:
                            "Default")
      --simulations SIMULATIONS
                            number of simulations in the JLA algorithm. It
                            defaults to 100 * log(#total fixed effect)
                            (type: Int64, default: 0)
      --header              CSV file contains header
      --first_id_display FIRST_ID_DISPLAY
                            The display text associated with first_id
                            (e.g. Person). (default: "Person")
      --second_id_display SECOND_ID_DISPLAY
                            The display text associated with second_id
                            (e.g. Firm) (default: "Firm")
      --observation_id_display OBSERVATION_ID_DISPLAY
                            The display text associated with observable_id
                            (e.g. Wage) (default: "Wage")
      --detailed_output_path DETAILED_OUTPUT_PATH
                            path to the CSV for the detailed output for
                            each observable (default:
                            "C:\\Users\\owner\\Desktop\\vchdfe\\variance_components.csv")
      --results_path RESULTS_PATH
                            path to the results of the output (default:
                            "C:\\Users\\owner\\Desktop\\vchdfe\\results.txt")
      --write_detailed_CSV  write the detailed output to a CSV
      --write_results       write the results to a file
      --print_level PRINT_LEVEL
                            Level of verbosity of output. (type: Int64,
                            default: 1)
      -h, --help            show this help message and exit
      
      

The executable can provide two different type of outputs. The first (that we refer to as Results) is a log file, in txt format, that writes a summary of the main results : bias-corrected components and some statistics from the leave-out sample. The second type of output (that we refer to as Detailed Output) is a DataFrame, stored in csv format, that saves important objects computed throughout the algorithm such a subset of the original data corresponding to the leave-out sample, the corresponding Pii and Bii, and the vector of coefficients used to compute the plug-in variance components. 

The program provides two outputs: a file contains the main results and a detailed output file. 

### Results output

By using the argument `--write_results` followed by `--results_path MYFILE`, you can specify that you want to write the results in `MYFILE`. 


### Detailed output

You can save all the output details in a CSV file. You need to use `--write_detailed_CSV` followed by `--detailed_output_path MYFILE`, where `MYFILE` is the path to the detailed output file. 

The detailed output file includes:                                                                                                                                                                                                                              
                           
- `observation` : observation identifier in the original dataset that belong to the Leave-out connected set. For instance, if the number 3 appears in this vector, it means that the third observation in the original dataset belongs to the leave-out connected set. 
- `first_id`: the first identifier corresponding to each observation in `obs` (e.g. worked ids in the leave-out connected set).
- `second_id`: the second identifier corresponding to each observation in `obs` (e.g. firm ids in the leave-out connected set).
- `beta`: the vector of coefficients used to compute the plug-in variance components.
- `D_alpha`: the fixed effect for the first identifier corresponding to each observation. 
- `F_psi`: the fixed effect for the second identifier corresponding to each observation. 
- `Pii`: statistical leverage corresponding to each observation in `obs`.
- `Bii_first`: The Bii for the variance of `first_id` effects corresponding to each observation in `obs`.
- `Bii_second`: The Bii for the variance of `second_id` effects corresponding to each observation in `obs`.
- `Bii_cov`: The Bii for the co-variance of `first_id` and `second_id` effects corresponding to each observation in `obs`.

## Examples 

Suppose we have a dataset `my_data.csv` that is stored in `"C:\\Users\\owner\\Desktop\\vchdfe\\"`. The structure of the data is such that the fifth column corresponds to the worker identifier and the third column corresponds to the firm identifier. Moreover, we want to store a summary of the results in this location `"C:\\Users\\owner\\Desktop\\vchdfe\\summary.txt"`, and the detailed output (Pii and Bii matrices, stored fixed effects, etc.) here  `"C:\\Users\\owner\\Desktop\\vchdfe\\output.csv"`. 


1. To obtain all three bias-corrected components (variance of worker effects, variance of firm effects, and covariance of worker-firm effects), we only need to type in the powershell 

```
vchdfe "C:\\Users\\owner\\Desktop\\vchdfe\\my_data.csv" --first_id 5 --second_id 3 --write_results --write_detailed_CSV --detailed_output_path "C:\\Users\\owner\\Desktop\\vchdfe\\output.csv" --results_path "C:\\Users\\owner\\Desktop\\vchdfe\\summary.txt"
```

2. To run the same thing while specifying 1000 simulations for the JLA algorithm that estimates (Pii,Bii) described in the computational appendix of KSS, we type in the powershell 

```
vchdfe "C:\\Users\\owner\\Desktop\\vchdfe\\my_data.csv" --first_id 5 --second_id 3 --write_results  --write_detailed_CSV --detailed_output_path "C:\\Users\\owner\\Desktop\\vchdfe\\output.csv" --results_path "C:\\Users\\owner\\Desktop\\vchdfe\\summary.txt" --algorithm JLA --simulations 1000
```

3. To only obtain the bias-correction for the variance of firm effects, we type in the powershell 

```
vchdfe "C:\\Users\\owner\\Desktop\\vchdfe\\my_data.csv" --first_id 5 --second_id 3 --no_first_effects --no_cov_effects --write_results  --write_detailed_CSV --detailed_output_path "C:\\Users\\owner\\Desktop\\vchdfe\\output.csv" --results_path "C:\\Users\\owner\\Desktop\\vchdfe\\summary.txt" --algorithm JLA --simulations 1000
```

## Functions in this package

### Main Functions 


```@docs
    leave_out_estimation(y,first_id,second_id,controls,settings)
    get_leave_one_out_set(y, first_id, second_id, settings, controls)    
```

### Auxiliary Functions

```@docs
    find_connected_set(y, first_idvar, second_idvar, settings)
    prunning_connected_set(yvec, first_idvar, second_idvar, obs_id, settings)
    drop_single_obs(yvec, first_idvar, second_idvar,obs_id)
    compute_movers(first_id,second_id)
    eff_res(::ExactAlgorithm, X,first_id,second_id,match_id, K, settings)
    eff_res(lev::JLAAlgorithm, X,first_id,second_id,match_id, K, settings)
    compute_matchid(second_id,first_id)   
```

## Datatypes in this package

```@docs
    ExactAlgorithm
    JLAAlgorithm
    VCHDFESettings
    
```


## About the current version


- The bias-correction currently only runs on a two-way fixed effects model without controls.

- The user can manually preadjust the outcome. For instance, in an AKM context, the user can run first
    $$y = pe + fe + Xb + e$$ , where $$pe$$ are person effects and $$fe$$ are firm effects, and feed into the routine $$y-Xb$$ as the outcome.
