Rate Stick Experiment Simulations for PBX 9501 with Various Jacket
Materials
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

This repository contains data regarding simulated rate stick experiments
for PBX9501 with $14$ different jacket materials. For an extended
abstract see the document `man/pbx9501.pdf` and for a description of the
simulation setup see the manuscript Rumsey et al. (2023). This
repository also contains scripts to reproduce the example from Section 4
of the aforementioned manuscript.

To run these scripts, you will need to install the free R package
`concordance` using the following command

``` r
# install.packages("devtools")
devtools::install_github("knrumsey/concordance")
```

The scripts included in this repository include

- `RateStickAnalysis.R` - Code to conduct the “Concordance Analysis” of
  Section 4. Reproduces Figures 5 and 6.
- `RateStickAnalysis_SS_U_Ni.R` - Code to conduct the analysis of
  Sections 4.1 and 4.2, regarding Uranium, Stainless Steel 304 and
  Nickel. Reproduces Tables 3, 4, and 5 and Figure 7.
- `RateStickAnalysis_CV.R` - Code to conduct the cross validation sim
  study which reproduces Table 5.

Note that the scripts must be run in the order listed above, as they
depend on previous values.

#### Release Information

*Rate Stick Experiment Simulations*

- LA-UR-23-33239
- 24-T-0194

*Co-Active Subspace Methods*

- LA-UR-23-33047
- 23-S-3475
