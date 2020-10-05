
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mvcapaCor

An `R` package for detecting (collective and point) anomalies or
changepoints in cross-correlated data. It also contains code to
reproduce the simulation study in “Scalable changepoint and anaomly
detection in cross-correlated data with an application to condition
monitoring”.

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

You can install tpcaMonitoring from github with:

``` r
# install.packages("devtools")
devtools::install_github("Tveten/mvcapaCor")
```

## Exported and documented functions

For more information, see the documentation of the functions below
inside R.

  - mvcapa\_cor
  - estimate\_precision\_mat
  - plot\_capa

## Example

``` r
library(mvcapaCor)
#> Loading required package: data.table
#> Loading required package: magrittr
#> Loading required package: tidyverse
#> ── Attaching packages ────────────────────────────────────────── tidyverse 1.3.0 ──
#> ✓ ggplot2 3.2.1     ✓ purrr   0.3.3
#> ✓ tibble  2.1.3     ✓ dplyr   0.8.3
#> ✓ tidyr   1.0.0     ✓ stringr 1.4.0
#> ✓ readr   1.3.1     ✓ forcats 0.4.0
#> ── Conflicts ───────────────────────────────────────────── tidyverse_conflicts() ──
#> x dplyr::between()   masks data.table::between()
#> x tidyr::extract()   masks magrittr::extract()
#> x dplyr::filter()    masks stats::filter()
#> x dplyr::first()     masks data.table::first()
#> x dplyr::lag()       masks stats::lag()
#> x dplyr::last()      masks data.table::last()
#> x purrr::set_names() masks magrittr::set_names()
#> x purrr::transpose() masks data.table::transpose()
#> Loading required package: lubridate
#> 
#> Attaching package: 'lubridate'
#> The following objects are masked from 'package:data.table':
#> 
#>     hour, isoweek, mday, minute, month, quarter, second, wday, week,
#>     yday, year
#> The following object is masked from 'package:base':
#> 
#>     date
#> Loading required package: ggvis
#> 
#> Attaching package: 'ggvis'
#> The following object is masked from 'package:ggplot2':
#> 
#>     resolution
# To run the entire simulation study, use run_simstudy().
# It will take several days to run, but the output from it is also contained in
# this package so that one can still get plots of results and study the results in more detail.
# The code in run_simstudy() specifies in a simple way the simulation
# setup and the necessary steps in the study.

# Generate the test sets used:
x <- simulate_cor(200, 4)
```
