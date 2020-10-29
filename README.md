# VarianceComponentsHDFE

[![Build Status](https://github.com/HighDimensionalEconLab/VarianceComponentsHDFE.jl/workflows/CI/badge.svg)](https://github.com/HighDimensionalEconLab/VarianceComponentsHDFE.jl/actions)
[![Coverage](https://codecov.io/gh/HighDimensionalEconLab/VarianceComponentsHDFE.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/HighDimensionalEconLab/VarianceComponentsHDFE.jl)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://HighDimensionalEconLab.github.io/VarianceComponentsHDFE.jl/dev)

- [Development Setup](develop.md)


The instructions below show how to manually install, and run the executable using a test sample. Follow the appropriate installation guide based on the operating system in use.

## Windows

### Installation Instructions for the Executable

1. Choose where you want to install the package, open up a powershell terminal (windows + R, then type "powershell"), and go to that location.

2. In the powershell, install the latest version:

```
wget https://github.com/HighDimensionalEconLab/VarianceComponentsHDFE.jl/releases/download/v0.1.5.1/vchdfe--windows-latest.tar.gz -O vchdfe-windows-latest.tar.gz
 
tar -xvf vchdfe-windows-latest.tar.gz
```
Note: To be able to use wget on Windows Internet Explorer should has been launched at least once. 

3. (OPTIONAL): Add the current directory to PATH:

```
setx PATH "$env:path;$pwd\vchdfe\bin" -m
```
Note: To permanently change the path, you need to run powershell as adminstrator. 

4. Test using the sample test file:

```
vchdfe vchdfe\bin\test.csv
```



## MacOS
### Installation Instructions for the Executable

1. Open Terminal: Press COMMAND + SPACE to open spotlight search, and type terminal and hit RETURN.

2. To download the compressed file in the home directory, run the following code:

```
cd
wget -qO- https://github.com/HighDimensionalEconLab/VarianceComponentsHDFE.jl/releases/download/v0.1.5.1/vchdfe-v0.1.5.1-macos-latest.tar.gz | tar -xzv
```

One can change the version number (i.e. v0.2) in the url to download a specific version. 

3. (OPTIONAL): To add the bin folder to the PATH, you have to modify the `.bash_profile` file in your home directory. Add the following line to the bottom of `.bash_profile` and save it. 
```
export PATH="~/vchdfe/bin:$PATH"
```

To source the changes in `.bash_profile` file, type:
```
source ~/.bash_profile
```


## Ubuntu
### Installation Instructions for the Executable

1. Open Terminal: CTRL + ALT + T

2. To download the compressed file in the home directory, run the following code:

```
cd ~
wget -qO- https://github.com/HighDimensionalEconLab/VarianceComponentsHDFE.jl/releases/download/v0.1.5.1/vchdfe-v0.1.5.1-ubuntu-latest.tar.gz | tar -xzv
```
One can change the version number (i.e. v0.2) in the url to download a specific version. 

3. (OPTIONAL): To add the bin folder to the PATH, you have to modify the `.bashrc` file in your home directory. Add the following line to the bottom of `.bashrc` and save it. 
```
export PATH="~/vchdfe/bin:$PATH"
```

To source the changes in `.bashrc` file, type:
```
source ~/.bashrc
```


# Testing

After installation, run the follwoing command to test:
```
vchdfe vchdfe/bin/test.csv 
```


# The Executable

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




