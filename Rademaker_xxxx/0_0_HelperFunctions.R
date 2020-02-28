################################################################################
#
#   Paper: "Fit indices for composite models in structural equation modeling"
#
#   Authors: Manuel Rademaker, Florian Schuberth, Sebastian Groﬂ
#
#   Topic: Helper functions
#
#   Last modified: 21.01.2020 (by Manuel Rademaker)
#
################################################################################
### General Preparation --------------------------------------------------------

## Install if not already installed
# if(!require(MASS)) install.packages("MASS")
# if(!require(dplyr)) install.packages("dplyr")

# Helper functions for the simulation ==========================================

DataGeneration <- function(.sample_size, .Sigma, .ndraws, .empirical = F,
                           .scale_factor = c("normal.kurt", "medium.kurt","high.kurt")) {
  #----------- Arguments:
  #   .sample_size : numeric value indicating the sample size to draw
  #   .Sigma       : (symmetric) indicator correlation matrix
  #   .ndraws      : number of replications that should be drawn
  #   .empirical   : see ?mvrnorm; defaults to "FALSE" 
  #   .scale_factor: the scaling factor to be used for each vector of indicator values. 
  #                  The scaling factor is used to induce non-normality in the data. If
  #                  "normal.kurt" (default) the data has a normal kurtosis (kurtosis of 3). 
  #                  "medium.kurt" increases the kurtosis by 1.7124. 
  #                  "high.kurt" increases the kurtosis by 6.
  
  #----------- What it does:
  #   Draws .ndraws samples from a multivariate normal distribution for a given sample_size
  #   and Sigma matrix. The resulting data is either multivariate normal or non-normal depending on 
  #   the .scale_factor argument. The function creates a list of lenght .ndraw (=the number of replications) 
  #   where each list element is a matrix of simulated data based on .Sigma. 
  #   The matrix has .sample_size rows and ncol columns. 
  
  type      <- match.arg(.scale_factor)
  scale_fac <- if(type == "normal.kurt") {
    
    1
  } else if(type == "medium.kurt") {
    
    sqrt(abs(rnorm(.sample_size))*sqrt(pi/2))
  } else if(type == "high.kurt") {
    
    sqrt(3)*(rchisq(.sample_size, 5))^(-1/2)
  }
  
  replicate(.ndraws,
            scale(scale_fac*MASS::mvrnorm(.sample_size, 
                                          mu = rep(0, ncol(.Sigma)), 
                                          Sigma = .Sigma,
                                          empirical = .empirical)),
            simplify = FALSE)
}
computeRelevant <- function(.object) {
  #----------- Arguments:
  #   .object            : a cSEMResults object
  
}