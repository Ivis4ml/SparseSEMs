## Introduction

`SparseSEMs` is a package that ultilizes the Proximal Alternating Linearized Maximal to solve the
non-convex non-smooth joint sparse structrual equation model. Presently the models that `SparseSEMs`
has implemented include fused-sparse structural equation model.

## Installation

`SparseSEMs` package contains a lot of necessary scripts to analyze large dataset such as microarray and SNP data
from GEO database, so it has not been submitted to CRAN yet for these non-standard directory.
To install `SparseSEMs`, you need a C++ compiler such as `g++` or `clang++` with C++11 feature,
and for Windows users, the [Rtools](https://cran.r-project.org/bin/windows/Rtools/index.html)
software is needed (unless you can configure the toolchain by yourself).

The installation follows the typical way of R packages on Github:

```r
library(devtools)
install_github("Ivis4ml/SparseSEMs")
```

