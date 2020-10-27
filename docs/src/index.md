

# VarianceComponentsHDFE

This package runs Kline, Saggio and Sølvsten (KSS) bias correction of variance components in two-way fixed effects models. The link to the repository is [this]( https://github.com/HighDimensionalEconLab/VarianceComponentsHDFE.jl).


## About the executable/app

The algorithm's purpose is to find two diagonal matrices that store all $$P_{ii}$$ and $$B_{ii}$$, and then proceed to compute the bias-correction of the (biased) plug-in estimators. To do so, you will only be required to run an executable file (app) where the user only needs to input the original data (in .csv format) that contains a first identifier (e.g. worker id), a second identifier (e.g firm id) and outcome (e.g. log wage).

With the use of random projections techniques, it is possible to run leave out estimation of variance components on large datasets. To give an idea, for a dataset containing 5 million person year observations, 1.3 million person effects, 90 thousand firm effects the code takes approximately 20 minutes to compute the relevant leave out matrices on a HPC system with 32 cores assigned. For more information on this please refer to Appendix B in the main paper. 

The algorithm prints the plug-in and the bias-corrected variance components estimators for the first identifier effects (e.g. variance of worker effects), the second identifier effects (e.g. variance of firm effects), and the covariance of both these effects (e.g. covariance of worker-firm effects). The user may choose to compute only a subset of these three components. Additionally, the executable will create a CSV file that stores vector of coefficients, the fixed effects for every observation in the leave-out connected set, as well as the diagonal matrices for $$P_{ii}$$s and $$B_{ii}$$s.

Use 

```
vchdfe --help
```

to see a complete list of available arguments:
  
    positional arguments:
      path                  path to CSV file containing data

    optional arguments:
      --first_id FIRST_ID   column index in CSV file for the first ID
                            (e.g. Person).  Use the most granular type.
                            (type: Int64, default: 1)
      --second_id SECOND_ID
                            column index in CSV file for the second ID
                            (e.g. Firm).  Use the less granular type.
                            (type: Int64, default: 2)
      --observation_id OBSERVATION_ID
                            column index in CSV file for observation (e.g.
                            Wage). (type: Int64, default: 4)
      --first_id_effects    Computing and showing first_id effects
      --cov_effects         Computing and showing covariace effects
      --algorithm ALGORITHM
                            type of algorithm: Exact or JLA. It defaults
                            to be Exact if the number of observations is
                            less than 5000, and JLA otherwise. (default:
                            "Default")
      --simulations SIMULATIONS
                            number of simulations in the JLA algorithm. If
                            0, defaults to 100 * log(#total fixed effect)
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
      
      
A detailed list of the output is provided below.
                      
                           
- `obs` : observation identifier in the original dataset that belong to the Leave-out connected set. For instance, if the number 3 appears in this vector, it means that the third observation in the original dataset belongs to the leave-out connected set. 
- `first_id`: the first identifier corresponding to each observation in `obs` (e.g. worked ids in the leave-out connected set).
- `second_id`: the second identifier corresponding to each observation in `obs` (e.g. firm ids in the leave-out connected set).
- `beta`: the vector of coefficients
- `D_alpha`: the fixed effect for the first identifier corresponding to each observation. We can obtain the plug-in estimator of the variance component of this effect by computing the variance of this vector.
- `F_psi`: the fixed effect for the second identifier corresponding to each observation. We can obtain the plug-in estimator of the variance component of this effect by computing the variance of this vector.
- `Pii`: statistical leverage corresponding to each observation in `obs`.
- `Bii_first`: The Bii to correct the variance of `first_id` effects corresponding to each observation in `obs`.
- `Bii_second`: The Bii to correct the variance of `second_id` effects corresponding to each observation in `obs`.
- `Bii_cov`: The Bii to correct the co-variance of `first_id` and `second_id` effects corresponding to each observation in `obs`.

## Functions in this package

```@docs
    leave_out_estimation(y,first_id,second_id,controls,settings)
    compute_movers(first_id,second_id)
    get_leave_one_out_set(y, first_id, second_id, settings, controls)
    
```

## Datatypes in this package

```@docs
    ExactAlgorithm
    JLAAlgorithm
    VCHDFESettings
    
```


## About the current version


- The bias-correction currently only runs on a model without controls.
- If the user wants, they can manually preadjust on her own the outcome. For instance, in an AKM context, the user can run first
    $$y = person effect + firm effect + Xb + e$$ and feed into the routine $$y-Xb$$ as the outcome.
