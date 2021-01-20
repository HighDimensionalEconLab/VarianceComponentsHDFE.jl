# Contents

```@contents
Pages = ["Executable.md"]
Depth = 3
```

# About the executable/app

The algorithm prints the plug-in and the bias-corrected variance components estimators for the first identifier effects (e.g. variance of worker effects), the second identifier effects (e.g. variance of firm effects), and the covariance of both these effects (e.g. covariance of worker-firm effects). The user may choose to compute only a subset of these three components. Additionally, the executable will create a CSV file that stores the outcome variable, the first and second set of identifiers, and the estimated fixed effects for every observation that belongs to the leave-out connected set as defined in KSS. The algorithm also saves the statistical leverages, Pii, as well the weighting terms defined as Bii in KSS for a given variance component. More details of this are provided below at the Results output and Detailed csv output sections. 

# Installation 

## Windows

1. Download our latest version of the package from the following [link](https://github.com/HighDimensionalEconLab/VarianceComponentsHDFE.jl/releases/download/v0.1.7.5/vchdfe--windows-latest.tar.gz). Move this file to `desired_installation_path`.
2. Open up a powershell terminal. We recommend to run powershell as administrator for installation. To do this, open Windows menu, type "powershell". Right-click on the powershell, click "run as administrator". 

3. Change the current directory to to the desired installation path by typing  in the powershell

   ```
   cd "desired_installation_path"
   ```

   Hint : To copy-paste into the terminal use the standard Ctrl+C and paste into the powershell by using right click.

4. In the powershell, install the latest version by running:

   ```
   tar -xvf vchdfe--windows-latest.tar.gz
   ```

5. (RECOMMENDED) A folder named `vchdfe/bin` now exists inside your installation path. Add the path of that folder to [PATH](#PATH-in-Windows). This will allow us to run the program everytime without specifying where the program is installed. 
6. (OPTIONAL) You can test the program using the sample test file provided with the executable. If you ran the previous step you may run:

   ```
   vchdfe vchdfe\bin\test.csv
   ```

   Otherwise, you will have to run

   ```
   .\vchdfe\bin\vchdfe .\vchdfe\bin\test.csv
   ```
   
   or you could test with some other file you store at the bin folder.

7. To set the number of threads used for parallel computing in the code, you need to use the set command before running `vchdfe` command. This might be very important if you intend to run bias-correction for very large datasets. Before running the program you may set the number of threads to, say 4, by entering the following line in the Powershell:

   ```
   $env:JULIA_NUM_THREADS=4
   ```
    
   Typically, you want to set that number to the number of cores in your computer. You can now proceed to close the terminal.

## PATH in Windows

The PATH is an important concept when working on the command line. It's a list of directories that tell your operating system where to look for programs, so that you can just write `program` instead of `some_folders\program`. But different operating systems have different ways to add a new directory to it. 

You can take a look at your current `PATH` by running the following line in the Powershell:

```
$env:path
```

We will describe two ways to add our program to `PATH`.

1. The first way consists on running  this line from the installation folder
   

   ```
   setx PATH "$env:path;$pwd\vchdfe\bin" -m
   ```

   Note: This change will be permanent only if you ran powershell as administrator. You may need to restart the powershell for this action to take place. You can confirm the change by running `$env:path` again, and checking whether the last part of this string contains `vchdfe\bin`.

2. If the previous step didn't work for you, you can still modify the `PATH` variable manually. To do so you can press (`Windows Key` + `R`) to open the Run Dialog, then type

   ```
   rundll32.exe sysdm.cpl,EditEnvironmentVariables
   ``` 
   and press the `Enter Key`. Find the `PATH` variable, select it, and click `Edit`. Then you can click `New` to add a new program path to the `PATH` list. You will write the path `desired_installation_path\vchdfe\bin`, where `desired_installation_path` is the folder where you installed the program. 

## MacOS

1. Download our latest version of the package from the following [link](https://github.com/HighDimensionalEconLab/VarianceComponentsHDFE.jl/releases/download/v0.1.7.5/vchdfe-v0.1.7.5-macos-latest.tar.gz). Move this file to the desired installation path.

2. You can unpack the .tar.gz file automatically by double-clicking the icon. 

3. A folder named `vchdfe/bin` now exists inside your installation path. You can double click the icon named vchdfe inside this folder. If the system warns you that it's not a recognized app you must click OK, to give permission.

4. (OPTIONAL) You can test the program using the sample test file provided with the executable. You may run the following line:

   ```
   vchdfe vchdfe/bin/test.csv
   ``` 

   or you can test with some other file you store at the bin folder. 

5. To set the number of threads used for parallel computing in the code, you need to use the set command before running `vchdfe` command. This might be very important if you intend to run bias-correction for very large datasets. Before running the program you may set the number of threads to, say 4, by entering the following line in the Powershell:

   ```
   set JULIA_NUM_THREADS=4
   ```
    
   Typically, you want to set that number to the number of cores in your computer. You can now proceed to close the terminal.


# Syntax of the program

Before this, make sure you've completed the installation steps shown above. To use the executable you only need to open the terminal. Windows users can press the keys Windows + R, then type "powershell" and press Enter. MacOS users can click to open the app manually. You can also see some examples in the Typical Workflow sections below.

The basic syntax of this command is 

```
vchdfe path_to_data [--option option_value]
```
where `[]` denote optional arguments. 

You can see a list of examples provided in a section below or type in the powershell

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
      --outcome_id OUTCOME_ID
                            column index in CSV file for outcome (e.g.
                            Wage). (type: Int64, default: 4)
      --no_first_id_effects No computing and showing of first_id effects
      --no_cov_effects      No computing and showing of covariance effects
      --algorithm ALGORITHM
                            type of algorithm: Exact or JLA. It defaults
                            to be Exact if the number of observations is
                            less than 5000, and JLA otherwise. (default:
                            "Default")
      --covariates NAME1 NAME2 ... 
                            Column names of the covariates that will be 
                            partialled out from the outcome. If dataset
                            contains no headers use Column1, Column2, etc.
                            to refer to such columns. (default: [])
      --do_lincom           Perform linear combination inference.
      --lincom_covariates NAME1 NAME2 ...
                            Column names of the covariates used in the linear
                            combination inference step. If dataset
                            contains no headers use Column1, Column2, etc.
                            to refer to such columns. (default: [])
      --simulations SIMULATIONS
                            number of simulations in the JLA algorithm. If
                            no number or zero is provided, it will default
                            to 200 simulations within the algorithm.
                            (type: Int64, default: 0)
      --header              CSV file contains header
      --first_id_display FIRST_ID_DISPLAY
                            The display text associated with first_id
                            (e.g. Person). (default: "Person")
      --second_id_display SECOND_ID_DISPLAY
                            The display text associated with second_id
                            (e.g. Firm) (default: "Firm")
      --outcome_id_display OUTCOME_ID_DISPLAY
                            The display text associated with outcome_id
                            (e.g. Wage) (default: "Wage")
      --detailed_csv_path DETAILED_CSV_PATH
                            path to the CSV for the detailed output for
                            each observable (default:
                            "current_directory/variance_components.csv")
      --results_path RESULTS_PATH
                            path to the results of the output (default:
                            "current_directory/results.txt")
      --write_detailed_csv  write the detailed output to a CSV
      --write_results       write the results to a file
      --print_level PRINT_LEVEL
                            Level of verbosity of output. (type: Int64,
                            default: 1)
      -h, --help            show this help message and exit
      
      

The executable can provide two different type of outputs. The first (that we refer to as Results) is a log file, in txt format, that writes a summary of the main results : bias-corrected components and some statistics from the leave-out sample. The second type of output (that we refer to as Detailed Output) is a DataFrame, stored in csv format, that saves important objects computed throughout the algorithm such a subset of the original data corresponding to the leave-out sample, the corresponding Pii and Bii, and the fixed effects that can be used to compute the plug-in variance components. 

The program provides two outputs: a file contains the main results and a detailed output file. 

## Results output

By using the argument `--write_results` followed by `--results_path MYFILE`, you can specify that you want to write the table with the main results in `MYFILE`. This will store some summary statistics of the leave out connected sample, the bias-corrected variance components. If `do_lincom` is activated, it will also store the results of the inference part, where it regresses second effects (e.g. firm effects) against some covariates.


## Detailed csv output

You can save all the output details in a CSV file. You need to use `--write_detailed_csv` followed by `--detailed_csv_path MYFILE`, where `MYFILE` is the path to the detailed output file. 

The detailed output file includes:                                                                                                                                                                                                                              
                           
- `observation` : observation identifier in the original dataset that belong to the Leave-out connected set. For instance, if the number 3 appears in this vector, it means that the third observation in the original dataset belongs to the leave-out connected set. 
- `first_id_old`: the first identifier of the original dataset corresponding to each observation in `obs` (e.g. worked ids in the leave-out connected set).
- `second_id_old`: the second identifier of the original dataset corresponding to each observation in `observation` (e.g. firm ids in the leave-out connected set).
- - `first_id`: the first identifier computed in the routine corresponding to each observation in `obs` (e.g. worked ids in the leave-out connected set). The maximum of this vector corresponds to the number of first effects (e.g. worker fixed effects).
- `second_id`: the second identifier computed in the routine corresponding to each observation in `observation` (e.g. firm ids in the leave-out connected set). The maximum of this vector corresponds to the number of second effects (e.g. firm fixed effects).
- `Dalpha`: the fixed effect for the first identifier corresponding to each observation. 
- `Fpsi`: the fixed effect for the second identifier corresponding to each observation. 
- `Pii`: statistical leverage corresponding to each observation in `observation`. If the leave out strategy was done at the match level this number is repeated for every observation that belongs to the same match.
- `Bii_first`: The Bii for the variance of `first_id` effects corresponding to each observation in `observation`. If the leave out strategy was done at the match level this number is repeated for every observation that belongs to the same match.
- `Bii_second`: The Bii for the variance of `second_id` effects corresponding to each observation in `observation`. If the leave out strategy was done at the match level this number is repeated for every observation that belongs to the same match.
- `Bii_cov`: The Bii for the co-variance of `first_id` and `second_id` effects corresponding to each observation in `observation`. If the leave out strategy was done at the match level this number is repeated for every observation that belongs to the same match.

# Typical Executable Workflow (Windows)

You begin by opening the Powershell (as administrator), and typing 

```
cd "path_to_dataset"
```    

At this point we recommend to exploit parallelization of the underlying code by changing an environment variable to be equal to the number of cores in your computer. Let `num_cores` be the number of cores. You can run the following line

```
$env:JULIA_NUM_THREADS=num_cores
```

Then, if you completed all instalation steps you may run

```
vchdfe my_data.csv --OPTIONS
```

where `OPTIONS` depend on the structure of your data. 

If you didn't complete step 5 of the installation, you will have to run instead 

```
installation_folder\vchdfe\bin\vchdfe my_data.csv --OPTIONS
```

## Detailed examples

Suppose we have a dataset `my_data.csv` that is stored in `"project_path"`. The structure of the data is such that the fifth column corresponds to the worker identifier, third column corresponds to the firm identifier, and the fourth column corresponds to the outcome. Moreover, we want to store a summary of the results in this location `"project_path/summary.txt"`, and the detailed output (Pii and Bii matrices, stored fixed effects, etc.) here  `"project_path/output.csv"`. 


1. To obtain all three bias-corrected components (variance of worker effects, variance of firm effects, and covariance of worker-firm effects), we only need to type in the terminal 

   ```
   cd project_path
   
   vchdfe my_data.csv --first_id 5 --second_id 3 --write_results --write_detailed_CSV --detailed_output_path output.csv --results_path summary.txt
   ```

2. To run the same thing while but we want to partial out the covariates in columns 7 and 8, we type in the terminal 

   ```
   cd project_path

   vchdfe my_data.csv --first_id 5 --second_id 3 --covariates Column7 Column8--write_results --write_detailed_CSV --detailed_output_path output.csv --results_path summary.txt 
   ```

3. To only obtain the bias-correction and regress the second effects (e.g. firm effects) onto the observables in columns 7 and 8 instead, we type in the terminal 

   ```
   cd project_path

   vchdfe my_data.csv --first_id 5 --second_id 3 --do_lincom --lincom_covariates Column7 Column8 --write_results --write_detailed_CSV --detailed_output_path output.csv --results_path summary.txt  
   ```


# Typical Executable Workflow (MacOS)

You begin by clicking the app that can be found inside the `vchdfe/bin` folder. You may need to click OK to give permission to use the app.  

At this point we recommend to exploit parallelization of the underlying code by changing an environment variable to be equal to the number of cores in your computer. Let `num_cores` be the number of cores. You can run the following line

```
set JULIA_NUM_THREADS=num_cores
```

Then, if you completed all instalation steps you may run

```
cd path_to_data
vchdfe my_data.csv --OPTIONS
```

where `OPTIONS` depend on the structure of your data. 


## Detailed examples

Suppose we have a dataset `my_data.csv` that is stored in `"project_path"`. The structure of the data is such that the fifth column corresponds to the worker identifier, the third column corresponds to the firm identifier and the fourth column corresponds to the outcome. Moreover, we want to store a summary of the results in this location `"project_path/summary.txt"`, and the detailed output (Pii and Bii matrices, stored fixed effects, etc.) here  `"project_path/output.csv"`.  You begin by opening the app as shown in the previous section.


1. To obtain all three bias-corrected components (variance of worker effects, variance of firm effects, and covariance of worker-firm effects), we only need to type in the terminal 

   ```
   cd project_path
   
   vchdfe my_data.csv --first_id 5 --second_id 3 --write_results --write_detailed_CSV --detailed_output_path output.csv --results_path summary.txt
   ```

2. To run the same thing while but we want to partial out the covariates in columns 7 and 8, we type in the terminal 

   ```
   cd project_path

   vchdfe my_data.csv --first_id 5 --second_id 3 --covariates Column7 Column8--write_results --write_detailed_CSV --detailed_output_path output.csv --results_path summary.txt 
   ```

3. To only obtain the bias-correction and regress the second effects (e.g. firm effects) onto the observables in columns 7 and 8 instead, we type in the terminal 

   ```
   cd project_path

   vchdfe my_data.csv --first_id 5 --second_id 3 --do_lincom --lincom_covariates Column7 Column8 --write_results --write_detailed_CSV --detailed_output_path output.csv --results_path summary.txt  
   ```


# About the current version

- The bias-correction currently only runs on a two-way fixed effects model without controls.

- The user can manually preadjust the outcome. For instance, in an AKM context, the user can run first
    $$y = pe + fe + Xb + e$$ , where $$pe$$ are person effects and $$fe$$ are firm effects, and feed into the routine $$y-Xb$$ as the outcome.
