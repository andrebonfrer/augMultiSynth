
# augMultiSynth

<!-- badges: start -->
<!-- badges: end -->

Multi-outcome augmented synthetic control for staggered adoption designs. 
This package is under development.

## Installation

You can install the development version of augMultiSynth like so:

``` r
devtools::install_github("andrebonfrer/augMultiSynth")
```

## Example

The repository includes additional simulation scripts under the `scripts/` directory.  
These are intended for reproducible examples and extended simulations that 
go beyond the package vignettes.

### `scripts/sim_example.R`

This script runs a full simulation–estimation–evaluation pipeline using
`augMultiSynth`, including:

- simulation of multi-outcome panel data under an interactive fixed-effects DGP,
- estimation via multi-outcome synthetic control,
- unit-level and avg(1..K) treatment effect evaluation,
- optional pooling / shrinkage steps.

Unlike the vignette, this script is designed to be:
- easily modified for large simulation runs,
- suitable for benchmarking and robustness checks,
- runnable non-interactively (e.g., on a server or cluster).

### How to run

From the repository root:

```r
library(augMultiSynth)

source("scripts/sim_example.R")
This is a basic example which shows you how to run staggered synthetic control,
producing individual treatment effects (and/or their weights). 

