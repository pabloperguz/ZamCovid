# ZamCovid <img src='man/figures/logo.png' align="right" height="150.5" />

<!-- badges: start -->
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
<!-- badges: end -->

ZamCovid is a mechanistic model of SARS-CoV-2 transmission (COVID-19) designed to be used with Zambian epidemiological surveillance data. It is implemented as stochastic compartmental models in `odin` and `dust`, and it uses `mcstate` to perfom parsimonious Bayesian evidence synthesis, requiring at a minimum seroprevalence and all deaths data. It is flexiblo to incorporate PCR-positive/LFT-positive cases data and COVID-19 confirmed deaths, when official public health guidelines match model space variable definitions (see `?ZamCovid::zamcovid_data()` for details).

## Installation

Contact `ZamCovid` developer if you require access to the package. After access is given, you can install the development version of `ZamCovid` using:

``` r
devtools::install_github("pabloperguz/ZamCovid")
```

## Important notes

Please note this package has been designed for very specific purposes, using selected datasets. Contact the developer before attempting to use `ZamCovid` for other projects or analyses.

`ZamCovid` has been developed with the use of [dust](https://mrc-ide.github.io/dust/), [mcstate](https://mrc-ide.github.io/mcstate/), [odin](https://mrc-ide.github.io/odin/), [odin.dust](https://mrc-ide.github.io/odin/). All these packages are under constant development and `ZamCovid` might not compile with versions other than the ones specified in `DESCRIPTION`.

For further information on the above packages and practical examples of their functionality see [FitzJohn et al., 2021. Reproducible parallel inference and simulation of stochastic state space models using odin, dust, and mcstate.](https://wellcomeopenresearch.org/articles/5-288/v2)

