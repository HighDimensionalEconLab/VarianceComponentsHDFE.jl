# VarianceComponentsHDFE

[![Build Status](https://github.com/HighDimensionalEconLab/VarianceComponentsHDFE.jl/workflows/CI/badge.svg)](https://github.com/HighDimensionalEconLab/VarianceComponentsHDFE.jl/actions)
[![Coverage](https://codecov.io/gh/HighDimensionalEconLab/VarianceComponentsHDFE.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/HighDimensionalEconLab/VarianceComponentsHDFE.jl)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://HighDimensionalEconLab.github.io/VarianceComponentsHDFE.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://HighDimensionalEconLab.github.io/VarianceComponentsHDFE.jl/dev)

- [Development Setup](develop.md)

## Using the Executable
The command line arguments are as follows:
```
    PathToExecutableDir/VarianceComponentsHDFEExecutable/bin/VarianceComponentsHDFE PathToCSVDataFile/data.csv --id=1 --firmid=2 --y=4 --algorithm=JLA --simulations=1000 --write_CSV --output_path=PathToCSVOutput/output.csv
```
  - The first argument is required and it is the path the the CSV file containing the data. The options `--id`, `--firmid`, `--y` indicate which columns of the CSV file contain the data for the worker IDs, firm IDs, and wages. `--algorithm` can be set to `Exact` or `JLA` and `--simulations` is the number of simulations in the JLA algorithm. `--write_CSV` is a flag that indicates the output will be written to a CSV file at `--output_path`. Additionally, you can run `PathToExecutableDir/VarianceComponentsHDFEExecutable/bin/VarianceComponentsHDFE --help` to see the arguments, options, and flags and their descriptions, and if applicable, default values.
  
## Installation Instructions for the Executable
To install a specifc version of the package, for example v0.2, and for a specific operating system: OSX, Windows, Ubuntu, one can run the following:

```
wget -qO- https://github.com/HighDimensionalEconLab/VarianceComponentsHDFE.jl/releases/download/v0.2/vchdfe-windows-latest.tar.gz | tar -xzv
```

