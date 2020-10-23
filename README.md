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
wget https://github.com/HighDimensionalEconLab/VarianceComponentsHDFE.jl/releases/download/v0.1.3/vchdfe--windows-latest.tar.gz -O vchdfe-windows-latest.tar.gz
 
tar -xvf vchdfe-ubuntu-latest.tar.gz
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


## Linux
### Installation Instructions for the Executable

1. Open Terminal: CTRL + ALT + T

2. To download the compressed file in the home directory, run the following code:

```
cd ~
wget -qO- https://github.com/HighDimensionalEconLab/VarianceComponentsHDFE.jl/releases/download/v0.1.3/vchdfe-v0.1.3-ubuntu-latest.tar.gz | tar -xzv
```
One can change the version number (i.e. v0.2) in the url to download a specific version. 

3. (OPTIONAL): To add the bin folder to the PATH, you have to modify the `.bashrc` file in your home directory. Add the following line to the bottom of `.bashrc` and save it. 
```
export PATH="~/vchdfe/bin:$PATH"
```

To sourse the changes in `.bashrc` file, type:
```
source ~/.bashrc
```


# Testing

After installation, run the follwoing command to test:
```
vchdfe vchdfe/bin/test.csv 
```


### The Executable

Use 

```
vchdfe --help
```

to see a complete list of available arguments. 




