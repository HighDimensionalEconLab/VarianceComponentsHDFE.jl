# Contents

```@contents
Pages = ["Executable.md"]
Depth = 3
```

# About the executable/app

The algorithm prints the plug-in and the bias-corrected variance components estimators for the first identifier effects (e.g. variance of worker effects), the second identifier effects (e.g. variance of firm effects), and the covariance of both these effects (e.g. covariance of worker-firm effects). The user may choose to compute only a subset of these three components. Additionally, the executable will create a CSV file that stores the outcome variable, the first and second set of identifiers, and the estimated fixed effects for every observation that belongs to the leave-out connected set as defined in KSS. The algorithm also saves the statistical leverages, Pii, as well the weighting terms defined as Bii in KSS for a given variance component.

# Installation 

## Windows

1. Download our latest version of the package from the following [link](https://github.com/HighDimensionalEconLab/VarianceComponentsHDFE.jl/releases/download/v0.1.5.6/vchdfe--windows-latest.tar.gz). Move this file to the desired installation path.
2. Open up a powershell terminal. We recommend to run powershell as administrator for installation. To do this, open Windows menu, type "powershell". Right-click on the powershell, click "run as administrator". 

3. Change the current directory to where you want to install the executable by typing  in the powershell

    ```
    cd "desired_installation_path"
    ```

    Hint : To copy-paste into the terminal use the standard Ctrl+C and paste into the powershell by using right click.


4. In the powershell, install the latest version by running:

    ```
    tar -xvf vchdfe--windows-latest.tar.gz
    ```

5. Add the installation directory to PATH. This will allow us to run the program everytime without specifying where the program is installed. To do so copy and paste the following line: 

    ```
    setx PATH "$env:path;$pwd\vchdfe\bin" -m
    ```

    Note: This change will be permanent only if you ran powershell as administrator. Otherwise, everytime you need to run the program you need to specify the installation folder : we would have to type  `"installation_path"\\vchdfe\\bin\\vchdfe` instead of `vchdfe` everytime we want to run the program. 


6. (OPTIONAL) You can test the program using the sample test file provided with the executable:

    ```
    vchdfe vchdfe\bin\test.csv
    ```

7. To set the number of threads used for parallel computing in the code, you need to use the set command before running `vchdfe` command. For example, to set the number of threads to 4, you may run the following code in the powershell:

    ```
    $env:JULIA_NUM_THREADS=4
    ```
    
You can now proceed to close the terminal.


## MacOS

1. Download our latest version of the package from the following [link](https://github.com/HighDimensionalEconLab/VarianceComponentsHDFE.jl/releases/download/v0.1.5.6/vchdfe-v0.1.5.6-macos-latest.tar.gz). Move this file to the desired installation path.
2. Open Terminal: Press COMMAND + SPACE to open spotlight search, and type terminal and hit RETURN.

3. You may unpack the .tar.gz file automatically when you double-click the icon. Otherwise, you may run the following code:

    ```
    cd desired_installation_path

    gunzip -c vchdfe-v0.1.5.3-macos-latest.tar.gz | tar xopft -
    ```

    And then close the terminal.

4. Add the installation directory to PATH. To do so you must open again a terminal window as in Step 1. Run 

    ```
    touch ~/.bash_profile; open ~/.bash_profile
    ```

    This will open a file known as `.bash_profile`, where you can copy and paste the following line  

    ```
    export PATH="desired_installation_path/vchdfe/bin:$PATH"
    ```
  
  
    where `desired_installation_path` is the folder where you installed the executable in the previous steps. Next we save the file we just modified and close it. Finally, to make sure that it will load those changes without rebooting the computer, run the following line in the terminal


    ```
    source ~/.bash_profile
    ```

5. (OPTIONAL) You can test the program using the sample test file provided with the executable:

    ```
    vchdfe vchdfe\bin\test.csv
    ```
6. To set the number of threads used for parallel computing in the code, you need to use the set command before running `vchdfe` command. For example, to set the number of threads to 4, you may run the following code in the command line:

      ```
      set JULIA_NUM_THREADS=4
      ```

You can now proceed to close the terminal.

# Executable guide

Before this, make sure you've completed the installation steps shown above. To use the executable you only need to open the terminal. Windows users can press the keys Windows + R, then type "powershell" and press Enter. MacOS users can press COMMAND + SPACE to open spotlight search, and type terminal and hit RETURN. 

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
      --outcome_id_display OUTCOME_ID_DISPLAY
                            The display text associated with outcome_id
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
      
      

The executable can provide two different type of outputs. The first (that we refer to as Results) is a log file, in txt format, that writes a summary of the main results : bias-corrected components and some statistics from the leave-out sample. The second type of output (that we refer to as Detailed Output) is a DataFrame, stored in csv format, that saves important objects computed throughout the algorithm such a subset of the original data corresponding to the leave-out sample, the corresponding Pii and Bii, and the fixed effects that can be used to compute the plug-in variance components. 

The program provides two outputs: a file contains the main results and a detailed output file. 

## Results output

By using the argument `--write_results` followed by `--results_path MYFILE`, you can specify that you want to write the table with the main results in `MYFILE`. 


## Detailed output

You can save all the output details in a CSV file. You need to use `--write_detailed_CSV` followed by `--detailed_output_path MYFILE`, where `MYFILE` is the path to the detailed output file. 

The detailed output file includes:                                                                                                                                                                                                                              
                           
- `observation` : observation identifier in the original dataset that belong to the Leave-out connected set. For instance, if the number 3 appears in this vector, it means that the third observation in the original dataset belongs to the leave-out connected set. 
- `first_id`: the first identifier corresponding to each observation in `obs` (e.g. worked ids in the leave-out connected set).
- `second_id`: the second identifier corresponding to each observation in `obs` (e.g. firm ids in the leave-out connected set).
- `D_alpha`: the fixed effect for the first identifier corresponding to each observation. 
- `F_psi`: the fixed effect for the second identifier corresponding to each observation. 
- `Pii`: statistical leverage corresponding to each observation in `obs`.
- `Bii_first`: The Bii for the variance of `first_id` effects corresponding to each observation in `obs`.
- `Bii_second`: The Bii for the variance of `second_id` effects corresponding to each observation in `obs`.
- `Bii_cov`: The Bii for the co-variance of `first_id` and `second_id` effects corresponding to each observation in `obs`.

# Examples 

Suppose we have a dataset `my_data.csv` that is stored in `"project_path"`. The structure of the data is such that the fifth column corresponds to the worker identifier and the third column corresponds to the firm identifier. Moreover, we want to store a summary of the results in this location `"project_path/summary.txt"`, and the detailed output (Pii and Bii matrices, stored fixed effects, etc.) here  `"project_path/output.csv"`. 


1. To obtain all three bias-corrected components (variance of worker effects, variance of firm effects, and covariance of worker-firm effects), we only need to type in the terminal 

```
cd project_path
vchdfe my_data.csv --first_id 5 --second_id 3 --write_results --write_detailed_CSV --detailed_output_path output.csv --results_path summary.txt
```

2. To run the same thing while specifying 1000 simulations for the JLA algorithm that estimates (Pii,Bii) described in the computational appendix of KSS, we type in the terminal 

```
cd project_path
vchdfe my_data.csv --first_id 5 --second_id 3 --write_results --write_detailed_CSV --detailed_output_path output.csv --results_path summary.txt --algorithm JLA --simulations 1000
```

3. To only obtain the bias-correction for the variance of firm effects, we type in the powershell 

```
cd project_path
vchdfe my_data.csv --first_id 5 --second_id 3 --no_first_effects --no_cov_effects --write_results --write_detailed_CSV --detailed_output_path output.csv --results_path summary.txt  --algorithm JLA --simulations 1000
```

# About the current version


- The bias-correction currently only runs on a two-way fixed effects model without controls.

- The user can manually preadjust the outcome. For instance, in an AKM context, the user can run first
    $$y = pe + fe + Xb + e$$ , where $$pe$$ are person effects and $$fe$$ are firm effects, and feed into the routine $$y-Xb$$ as the outcome.
