################################################################################
#
#   Paper: "Fit indices for composite models in structural equation modeling"

#   Authors: Manuel Rademaker, Florian Schuberth, Sebastian Gro?
#
#   Last modified: 27.02.2020 (by Manuel Rademaker)
#
#   Topic: Simulation
#
#   Notes:
#
################################################################################
rm(list = ls())

# General Preparation ==========================================================
## Install packages if necessary
if(!require(foreach)) install.packages("foreach")
if(!require(doParallel)) install.packages("doParallel")
if(!require(parallel)) install.packages("parallel")

## Load packages 
library(purrr)
library(cSEM)       # (version: 0.1.0.9000; Feb 2020)
library(foreach)    

## Load DGPs and Models
load("DGPs/dgps.RData")
load("Models/Models.RData")

## Source helper functions and models to estimate
source("0_0_HelperFunctions.R")

# Simulation ===================================================================
### Preparation ----------------------------------------------------------------
# Things to loop over (full design)
sample_size      <- seq(100, 1000, by = 50)
number_of_draws  <- 500
number_of_boot_reps <- 499
dgp_and_model    <- purrr::map(purrr::transpose(list("Dgp" = dgps, "Model" = models)), purrr::transpose)

#Small design
# sample_size     <- c(100, 200, 500)
# number_of_draws <- 3
# number_of_boot_reps <- 10
# dgp_and_model   <- dgp_and_model[c(1, 3, 6)] %>%
#   map(.f = ~ .x[c(2, 4)])
# listviewer::jsonedit(dgp_and_model, mode = "view")

### Monte Carlo simulation -----------------------------------------------------
## List hierachy of the resulting simulation object (bottom to top)
number_cores    <- parallel::detectCores() 
## Create cluster 
cl <- parallel::makeCluster(number_cores)
doParallel::registerDoParallel(cl)

sim <- 
  foreach(n = sample_size, .packages = c("MASS", "dplyr", "cSEM")) %:%
  foreach(i = seq_along(dgp_and_model)) %:% # index for number of composites
  foreach(j = seq_along(dgp_and_model[[1]])) %dopar% { # index for number of indicators
    
    # Generate data
    dat <- DataGeneration(.sample_size = n, .Sigma = dgp_and_model[[i]][[j]]$Dgp, .ndraws = number_of_draws)
    
    # Estimation 
    est <- csem(dat, dgp_and_model[[i]][[j]]$Model, .PLS_weight_scheme_inner = "centroid")
    
    # Delete estimations that did not pass verify()
    est[] <- est[unlist(lapply(verify(est), function(x) sum(x) == 0))]
    
    # Compute fit statistics using assess
    out_assess <- assess(est, .quality_criterion =  c("chi_square", "chi_square_df",
                                                      "cfi", "gfi", "ifi", "nfi", "nnfi", 
                                                      "rmsea", "rms_theta", "rms_theta_mi", 
                                                      "srmr"))
    out_test <- testOMF(est, .fit_measures = TRUE, .R = number_of_boot_reps,
                        .handle_inadmissibles = "replace", .verbose = FALSE,
                        .alpha = c(0.01, 0.05, 0.1))
    out_test <- lapply(out_test, function(x) {
      y <- data.frame(
        "Distance_measure" = names(x$Test_statistic),
        stringsAsFactors = FALSE
      )
      y <- cbind(y, x$Decision)
      rownames(y) <- NULL
      y
    })
    out <- list(
      "Assess" = out_assess,
      "Test"   = out_test
    )
  } # END simulation


closeAllConnections() # close connection to relase RAM

## Inspect resulting object (only if object is small. If number_of_draws is
## large this will cause the session to abort).

# listviewer::jsonedit(sim, mode = "view")

# Save objects =================================================================

save(list = c("sim", "sample_size"),
     file = "sim_hpc.RData")

#
# i = 3; j = 3; n = 200
# dat <- DataGeneration(.sample_size = n, .Sigma = dgp_and_model[[i]][[j]]$Dgp, .ndraws = 100, .empirical = FALSE)
# 
# est <- csem(dat, dgp_and_model[[i]][[j]]$Model, .PLS_weight_scheme_inner = "centroid")
# est[] <- est[unlist(lapply(verify(est), function(x) sum(x) == 0))]
# 
# out_assess <- assess(est, .quality_criterion =  c("chi_square", "chi_square_df",
#                                                   "cfi", "gfi", "ifi", "nfi", "nnfi", 
#                                                   "rmsea", "rms_theta", "rms_theta_mi", 
#                                                   "srmr"))
# 
# out_test <- testOMF(est, .fit_measures = TRUE, .R = 3,
#                     .handle_inadmissibles = "replace", .verbose = FALSE,
#                     .alpha = c(0.01, 0.05, 0.1))
# 
# out_test$Data_1$Decision
# 
# out_test <- lapply(out_test, function(x) {
#   y <- data.frame(
#     "Distance_measure" = names(x$Test_statistic),
#     stringsAsFactors = FALSE
#   )
#   y <- cbind(y, x$Decision)
#   rownames(y) <- NULL
#   y
# })
