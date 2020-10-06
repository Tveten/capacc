
<!-- README.md is generated from README.Rmd. Please edit that file -->

# capacc

An `R` package for detecting (collective and point) anomalies (CAPA-CC)
or changepoints (CPT-CC) in cross-correlated data. It also contains code
to reproduce the simulation study in Tveten, Eckley, Fearnhead (2020)
“Scalable changepoint and anaomly detection in cross-correlated data
with an application to condition monitoring”.

## Overview

Functionality:

  - Function for running the CAPA-CC algorithm on data, as well as
    CPT-CC for a single change.
  - Functions for visualising the detected collective and point
    anomalies.
  - Functions for estimating a precision matrix restricted to a given
    adjacency matrix by the GLASSO method.
  - Functions for generating test data and running the simulation study
    in “Scalable changepoint and anaomly detection in cross-correlated
    data with an application to condition monitoring”.

## Installation

You can install capacc from github with:

``` r
# install.packages("devtools")
devtools::install_github("Tveten/capacc")
```

## Exported and documented functions

For more information, see the documentation of the functions below
inside R.

  - capacc
  - cptcc
  - estimate\_precision\_mat
  - plot\_capa

## Example

``` r
library(capacc)

# Generate the test sets used:
x <- simulate_cor(200, 4)
```
