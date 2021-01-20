# VarianceComponentsHDFE

[![Build Status](https://github.com/HighDimensionalEconLab/VarianceComponentsHDFE.jl/workflows/CI/badge.svg)](https://github.com/HighDimensionalEconLab/VarianceComponentsHDFE.jl/actions)
[![Coverage](https://codecov.io/gh/HighDimensionalEconLab/VarianceComponentsHDFE.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/HighDimensionalEconLab/VarianceComponentsHDFE.jl)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://HighDimensionalEconLab.github.io/VarianceComponentsHDFE.jl/dev)

This Julia package implements the leave out correction of 
[Kline, Saggio and Soelvsten (2020)](https://eml.berkeley.edu/~pkline/papers/KSS2020.pdf) for estimating variance components in two-way fixed effects models. We provide the usual Julia Package as well as an executable/app that can be run in the terminal and does not require the installation of Julia or any other statistical software.

For instructions on how to install and use the executable follow this [link](https://highdimensionaleconlab.github.io/VarianceComponentsHDFE.jl/dev/Executable/).

For instructions on how to install the Julia Package and know about the exported methods follow this [link](https://highdimensionaleconlab.github.io/VarianceComponentsHDFE.jl/dev/Package/).

## Julia Package Development

Some instructions are provided at the [Development Setup](develop.md) and at the [Package docs](https://highdimensionaleconlab.github.io/VarianceComponentsHDFE.jl/dev/Package/).

## Matlab Version

The Matlab version of the package can be 
found [here](https://github.com/rsaggio87/LeaveOutTwoWay).

## Summary of Latest Version of Executable

* By default, the code runs a leave-out correction by leaving a match out as opposed to leaving an observation out. See Matlab's [vignette](https://github.com/rsaggio87/LeaveOutTwoWay/blob/master/doc/VIGNETTE.pdf) for details.
* New, optimized, random projection algorithm for calculation of the statistical leverages which scale extremely well to large datasets.
* Controls are partialled out before performing the KSS routine.
* High dimensional linear combination can be performed on the fixed effects based on the second identifier (e.g. firm effects). 

