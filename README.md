# regCtsem

Regularized Continuous Time Structural Equation Modeling

# Installation

If you want to install regCtsem from GitHub, use the following commands in R:

    if(!require(devtools))install.packages("devtools")

    devtools::install_github("jhorzek/regCtsem")
    
In case you also want to install the vignettes, use:


    if(!require(devtools))install.packages("devtools")

    devtools::install_github("jhorzek/regCtsem", build_vignettes = TRUE)
    
    
# Getting Started

A good place to start is the help page of the main function:

    ?regCtsem::regCtsem

If you also installed the vignettes, use ``vignette("regCtsem", package = "regCtsem")'' to get some more examples and general guidelines.
