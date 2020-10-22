# VarianceComponentsHDFE

[![Build Status](https://github.com/HighDimensionalEconLab/VarianceComponentsHDFE.jl/workflows/CI/badge.svg)](https://github.com/HighDimensionalEconLab/VarianceComponentsHDFE.jl/actions)
[![Coverage](https://codecov.io/gh/HighDimensionalEconLab/VarianceComponentsHDFE.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/HighDimensionalEconLab/VarianceComponentsHDFE.jl)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://HighDimensionalEconLab.github.io/VarianceComponentsHDFE.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://HighDimensionalEconLab.github.io/VarianceComponentsHDFE.jl/dev)

- [Development Setup](develop.md)


The instructions below show how to manually install, and run the executable using a test sample. Follow the appropriate installation guide based on the operating system in use.

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
vchdfe.exe vchdfe\bin\test.csv --algorithm=JLAAlgorithm
```

### Using the Executable

Use 

```
vchdfe.exe --help
```

to see a complete list of available arguments. 



## MacOS
### Installation Instructions for the Executable

### Using the Executable


## Linux
### Installation Instructions for the Executable

1. Open Terminal: CTRL + ALT + T

2. To download the compressed file in the home directory, run the following code:

```
cd ~
wget -qO- https://github.com/HighDimensionalEconLab/VarianceComponentsHDFE.jl/releases/download/v0.2/vchdfe-windows-latest.tar.gz | tar -xzv
```
One can change the version number (i.e. v0.2) in the url to download a specific version. 

3. OPTIONAL: To add the bin folder to the PATH, you have to modify the `.bashrc` file in your home directory. Add 
```
export PATH="/HOMEDIRECTORY/vchdfe/bin:$PATH"
```
where `HOMEDIRECTORY` is the absolute address of your home directory. You can find this address by running `pwd()` in the home directory. 

To sourse the changes in `.bashrc` file, type:
```
source ~/.bashrc
```

4. Test using the sample test file:
``` 
vchdfe vchdfe\bin\test.csv --algorithm=JLAAlgorithm
```

### Using the Executable

Use 

```
vchdfe.exe --help
```

to see a complete list of available arguments. 




