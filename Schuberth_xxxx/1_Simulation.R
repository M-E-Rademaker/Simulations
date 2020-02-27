################################################################################
#
#   Paper: "Composite of composites"
#
#   Topic: Monte Carlo Simulation
#
#   Last modified: 17.02.2020 (by Manuel Rademaker)
#
#   Notes: none
#
################################################################################
rm(list = ls())

# General Preparation ==========================================================
## Install packages if necessary
#if(!require(foreach)) install.packages("foreach")
#if(!require(doParallel)) install.packages("doParallel")
#if(!require(parallel)) install.packages("parallel")
#if(!require(cSEM)) install.packages("cSEM")

## Load packages 
library(cSEM)       # (version: 0.1.0)
library(foreach)    

## Load DGPs and Models
load("DGPs/DGPs.RData")
load("Models/Models.RData")

## Source helper functions and models to estimate
source("0_0_HelperFunctions.R")

# Simulation ===================================================================
### Preparation ----------------------------------------------------------------
# Things to loop over 
sample_size      <- list(100, 300, 500, 1000)
number_of_draws  <- 200
number_boot_reps <- 1000
dgp              <- list(
  "Sigma_dgp_2ndorder" = Sigma_dgp_2ndorder, 
  "Sigma_dgp_alt_low"  = Sigma_dgp_alt_small,
  "Sigma_dgp_alt_med"  = Sigma_dgp_alt_med,
  "Sigma_dgp_alt_high" = Sigma_dgp_alt_large)
model            <- mget(c("model_2ndorder", "model_alt"))
weighting_scheme <- c("centroid", "factorial", "path") 

## Test 
#sample_size      <- list(100)
#number_of_draws  <- 3
#number_boot_reps <- 10
#weighting_scheme <- c("centroid")
number_cores     <- parallel::detectCores()
number_cores
### Monte Carlo simulation -----------------------------------------------------

## List hierachy of the resulting simulation object (bottom to top)

# - A cSEMResults_multi objects with as many cSEMResults objects as 
#   number_of_draws where used.
# - Sample size
# - Weightning scheme
# - Approach
# - DGP

## Create cluster 
cl <- parallel::makeCluster(number_cores)
doParallel::registerDoParallel(cl)

sim <- 
  foreach(j = seq_along(dgp), .packages = c("MASS", "dplyr", "cSEM")) %:%
  foreach(m = seq_along(model)) %:%
  foreach(l = weighting_scheme) %:%
  foreach(i = sample_size) %dopar%  {
    
    # Create dummy structures make data processing
    # easier
    
    out1 <- data.frame(
      "Estimate_type"  = NA,
      "Name"           = NA,
      "Estimate"       = NA
    )
    out2 <- data.frame(
      "Stage"            = NA,
      "Distance_measure" = NA,
      "99%"              = NA,
      "95%"              = NA,
      "90%"              = NA,
      check.names = FALSE
    )
    
    out_dummy <- list(
      "Parameter_estimates" = out1, 
      "Test"                = out2
    )
    
    out_dummy_status <- data.frame(
      "Status_code"   = NA,
      "Status_code_1" = NA,
      "Status_code_2" = NA
    ) 
    
    out_status <- list(
      list(
        out_dummy_status, 
        out_dummy_status, 
        out_dummy_status 
      )
    )
    
    ## Note: 
    ##  j := the index for the DGP
    ##  m := the index for the model
    out_estimation  <- list()
    counter <- 0
    repeat{
      
      counter <- counter + 1
      
      ## Generate data
      samp <- DataGeneration(
        .sample_size  = i,
        .Sigma        = dgp[[j]],
        .ndraws       = 1,
        .scale_factor = "normal.kurt",
        .empirical    = FALSE
      )[[1]]
      
      ## The (extended) repeated indicators approach requires the repeated indicators to
      ## be attached to the data set
      coln <- c(colnames(samp), paste0(colnames(samp), "_temp"))
      samp_RI <- samp[, c(1:ncol(samp), 1:ncol(samp))]
      colnames(samp_RI) <- coln
      
      ### Two-stage approach ---------------------------------------------------
      ## Do the estimation (Result is a list of length number_of_draws with
      ## elements First_stage and Second_stage)
      out_2stage <- csem(
        .data  = samp,
        .model =  model[[m]],
        .approach_2ndorder = "2stage",
        .PLS_weight_scheme_inner = l
      )
      
      status_2stage <- sum(unlist(verify(out_2stage)))
      
      ## Only compute extended repeated indicators approach and mixed approach 
      ## for the 2nd order model (model_pop_Sigma_2ndorder) and if the DGP is 
      ## the second order DGP (Sigma_2ndorder)
      status_mixed <- 0
      status_RI_extended <- 0
      if(j == 1 & m == 1) {
        
        ### Mixed two-stage approach -------------------------------------------
        ## Do the estimation (Result is a list of length number_of_draws with
        ## elements First_stage and Second_stage)
        out_mixed <- csem(
          .data  = samp,
          .model = model[[m]],
          .approach_2ndorder = "mixed",
          .PLS_weight_scheme_inner = l
        ) 
        
        status_mixed <- sum(unlist(verify(out_mixed)))
        
        ### Extended repeated indicators approach ------------------------------
        out_RI_extended <- csem(
          .data  = samp_RI,
          .model = model_RI_extended,
          .PLS_weight_scheme_inner = l
        )
        
        status_RI_extended <- sum(unlist(verify(out_RI_extended)))
        
      } else {
        ## Set up dummy structure to facilitate data processing
        out_mixed       <- out_dummy
        out_RI_extended <- out_dummy
        
      } # END if j == 1 & m == 1 
      
      if(status_2stage == 0 & status_mixed == 0 & status_RI_extended == 0) {
        
        out_2stage <- computeRelevant(out_2stage, .R = number_boot_reps)
        if(j == 1 & m == 1) {
          out_mixed       <- computeRelevant(out_mixed, .R = number_boot_reps)
          out_RI_extended <- computeRelevant(out_RI_extended, .R = number_boot_reps) 
        }
        
        out_estimation[[counter]] <- list(
          "two-stage"   = out_2stage, 
          "embedded_TS" = out_mixed, 
          "extended_RI" = out_RI_extended
        )   
      } else {# status is not ok
        ll <- list(out_2stage, out_mixed, out_RI_extended)
        
        out_status[[length(out_status) + 1]] <- lapply(ll, function(x) {
          if(inherits(x, "cSEMResults_default")) {
            data.frame(
              "Status_code"   = names(verify(x)),
              "Status_code_1" = c(verify(x)),
              "Status_code_2" = NA
            )  
          } else if(inherits(x, "cSEMResults_2ndorder")) {
            data.frame(
              "Status_code"   = names(verify(x)[[1]]),
              "Status_code_1" = c(verify(x)[[1]]),
              "Status_code_2" = c(verify(x)[[2]])
            )
          } else {
            data.frame(
              "Status_code"   = NA,
              "Status_code_1" = NA,
              "Status_code_2" = NA
            )  
          }
        })
        # Reset counter
        counter <- counter - 1
      }
      
      # Break repeat loop if .R results have been created.
      if(length(out_estimation) == number_of_draws) {
        ## Give names
        names(out_status) <- 1:length(out_status)
        break
      }
    } # END repeat
    
    out <- list(
      "out_estimation"  = out_estimation,
      "out_status"      = out_status
    )
  } # END simulation

closeAllConnections() # close connection to relase RAM

## Inspect resulting object (only if object is small. If number_of_draws is
## large this will cause the session to abort).

# listviewer::jsonedit(sim, mode = "view") 

# Save objects =================================================================

save(list = c("sim", "sample_size", "number_of_draws"),
     file = "sim_hpc_200runs2.RData")
