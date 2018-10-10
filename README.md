# MASC: Mixed-effects association testing for single cells
MASC is a novel reverse single cell association strategy for testing whether a specified covariate influences the membership of single cells in any of multiple cellular subsets while accounting for technical confounds and biological variation.
## Current Version
MASC 0.1.0

## Requirements
MASC is written for **R 3.5**. It requires the following package:
* lme4

## Installation
```R
install.packages("devtools")
library(devtools)
install_github("immunogenomics/masc")
```

## Usage
MASC expects a data frame that contains, at minimum, a factor indicating cluster membership for single cells, a factor representing the covariate of interest, and other random- and fixed-effects covariates.
These latter terms should be given as character vectors that represent the name of the column with the covariate information in the input data frame.

*Work in progress*
