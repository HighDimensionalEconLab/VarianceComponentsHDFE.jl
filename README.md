# VarianceComponentsHDFE

[![Build Status](https://github.com/HighDimensionalEconLab/VarianceComponentsHDFE.jl/workflows/CI/badge.svg)](https://github.com/HighDimensionalEconLab/VarianceComponentsHDFE.jl/actions)
[![Coverage](https://codecov.io/gh/HighDimensionalEconLab/VarianceComponentsHDFE.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/HighDimensionalEconLab/VarianceComponentsHDFE.jl)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://HighDimensionalEconLab.github.io/VarianceComponentsHDFE.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://HighDimensionalEconLab.github.io/VarianceComponentsHDFE.jl/dev)

- [Development Setup](develop.md)


The instructions below show how to download, install, and run the executable using a test sample. Follow the installation guide based on the operating system in use.

## Linux
### Installation Instructions for the Executable

Open Terminal: ctrl + alt + T

In the terminal, run:

```
wget -qO- https://github.com/HighDimensionalEconLab/VarianceComponentsHDFE.jl/releases/download/v0.2/vchdfe-windows-latest.tar.gz | tar -xzv
```
This command downloads the binary files. One can change the version number (i.e. v0.2) in the url to download a specific version. 

### Using the Executable

In the same terminal, use the command below to run the executable:

```
vchdfe\bin\vchdfe.exe vchdfe\bin\test.csv --algorithm=JLAAlgorithm
```

## Windows

### Installation Instructions for the Executable

1. Choose where you want to install the package, open up a powershell terminal (windows + R, then type "powershell"), and go to that location.

2. In the powershell, install the latest version:

```
wget -qO- https://github.com/HighDimensionalEconLab/VarianceComponentsHDFE.jl/releases/download/v0.2/vchdfe-windows-latest.tar.gz | tar -xzv
```
3. Add the current directory to PATH:

```
$env:Path += ";$pwd\vchdfe\bin"
```

4. Test using the sample test file:

```
vchdfe\bin\vchdfe.exe vchdfe\bin\test.csv --algorithm=JLAAlgorithm
```

### Using the Executable

In the same powershell, use the command below to run the executable:

```
vchdfe.exe vchdfe\bin\test.csv --algorithm=JLAAlgorithm
```



## MacOS
### Installation Instructions for the Executable

### Using the Executable



The command line arguments are as follows:
```
    PathToExecutableDir/VarianceComponentsHDFEExecutable/bin/VarianceComponentsHDFE PathToCSVDataFile/data.csv --id=1 --firmid=2 --y=4 --algorithm=JLA --simulations=1000 --write_CSV --output_path=PathToCSVOutput/output.csv
```
  - The first argument is required and it is the path the the CSV file containing the data. The options `--id`, `--firmid`, `--y` indicate which columns of the CSV file contain the data for the worker IDs, firm IDs, and wages. `--algorithm` can be set to `Exact` or `JLA` and `--simulations` is the number of simulations in the JLA algorithm. `--write_CSV` is a flag that indicates the output will be written to a CSV file at `--output_path`. Additionally, you can run `PathToExecutableDir/VarianceComponentsHDFEExecutable/bin/VarianceComponentsHDFE --help` to see the arguments, options, and flags and their descriptions, and if applicable, default values.

To install a specifc version of the package, for example v0.2, and for a specific operating system: OSX, Windows, Ubuntu, one can run the following:

```
wget -qO- https://github.com/HighDimensionalEconLab/VarianceComponentsHDFE.jl/releases/download/v0.2/vchdfe-windows-latest.tar.gz | tar -xzv
```

