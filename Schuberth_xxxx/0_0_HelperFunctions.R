################################################################################
#
#   Paper: "Composite of composites"
#
#   Authors: Florian Schuberth, Manuel Rademaker, Joerg Henseler
#
#   Topic: Helper functions
#
#   Last modified: 17.02.2020 (by Manuel Rademaker)
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

computeRelevant <- function(.object, .R) {
  
  # to do
  # - alle weights
  # - alternatice DGPs mit alternative models
  
  #----------- Arguments:
  #   .object            : a cSEMResults object
  
  ## Which columns are relevant in Path_ , Weight_ , and Effect_estimates?
  relevant <- c("Name", "Estimate")
  
  ## Which parameters are of interest for the simulation
  # All path          : gamma1, gamma2, gamma3, beta1, beta2
  # All weights       : w_y11, w_y21, w_y31, w_y41, w_y51, w_y61, w_c1, w_c2, w_c3
  # All total_effects : te_eta1_eta4, te_eta2_eta4, te_eta3_eta4
  
  selector_path <- c(
    "eta1 ~ xi", 
    "eta2 ~ xi", 
    "eta3 ~ xi", 
    "eta3 ~ eta1", 
    "eta3 ~ eta2"
    )
  selector_weights_1stage <- c(
    "c1 <~ y11", "c1 <~ y12", 
    "c2 <~ y21", "c2 <~ y22", "c2 <~ y23", "c2 <~ y24", 
    "c3 <~ y31", "c3 <~ y32", "c3 <~ y33", "c3 <~ y34", 
    "c3 <~ y35", "c3 <~ y36", "c3 <~ y37", "c3 <~ y38",
    "xi <~ y41", "xi <~ y42", "xi <~ y43",
    "eta1 <~ y51", "eta1 <~ y52", "eta1 <~ y53", 
    "eta2 <~ y61", "eta2 <~ y62", "eta2 <~ y63"
    )
  selector_weights_2ndstage <- c("eta3 <~ c1", "eta3 <~ c2", "eta3 <~ c3")
  selector_te <- c("eta3 ~ xi", "eta3 ~ eta1", "eta3 ~ eta2")
  
  if(inherits(.object, "cSEMResults_2ndorder")) {
    
    ### Parameter estimates --------------------------------------------------
    sum_res <- cSEM:::summarize(.object)
    path     <- sum_res$Second_stage$Estimates$Path_estimates[
      sum_res$Second_stage$Estimates$Path_estimates$Name %in% selector_path, relevant, drop = FALSE]
    weights_1stage <- sum_res$First_stage$Estimates$Weight_estimates[
      sum_res$First_stage$Estimates$Weight_estimates$Name %in% selector_weights_1stage, relevant, drop = FALSE]
    weights_2stage <- sum_res$Second_stage$Estimates$Weight_estimates[
      sum_res$Second_stage$Estimates$Weight_estimates$Name %in% selector_weights_2ndstage, relevant, drop = FALSE]
    total_effect <- sum_res$Second_stage$Estimates$Effect_estimates$Total_effect[
      sum_res$Second_stage$Estimates$Effect_estimates$Total_effect$Name %in% selector_te, relevant, drop = FALSE]
    
    out1 <- dplyr::bind_rows(
      list("Path" = path, 
           "Weights" = dplyr::bind_rows(weights_1stage, weights_2stage), 
           "Total_effect" = total_effect), 
      .id = "Estimate_type")
    
    ### Overall model fit test -------------------------------------------------
    ## Compute test only for two-stage approach
    if(.object$Second_stage$Information$Approach_2ndorder == "2stage") {
      
      ## Joint test
      out2a <- testOMF(.object, 
                       .handle_inadmissibles = "replace",
                       .alpha = c(0.1, 0.05, 0.01),
                       .verbose = FALSE, .R = .R)
      
      out2b <- data.frame(
        "Distance_measure" = names(out2a$Test_statistic),
        stringsAsFactors = FALSE
      )
      
      out2a <- cbind(out2b, out2a$Decision)
      rownames(out2a) <- NULL
      out2 <- cbind("Stage" = "Both", out2a) 
      
      ## Test first and second stage separately
      
      out2c <- lapply(.object, function(x) {
        x <- testOMF(x, .handle_inadmissibles = "replace", .alpha = c(0.1, 0.05, 0.01), 
                     .verbose = FALSE, .R = .R)
        xa <- data.frame(
          "Distance_measure" = names(x$Test_statistic),
          stringsAsFactors = FALSE
        )
        
        x <- cbind(xa, x$Decision)
        rownames(x) <- NULL
        x
      })
      ## Combine
      out2c[["Both"]] <- out2a
      out2 <- dplyr::bind_rows(out2c, .id = "Stage")
    } else {
      ## Set up empty data set to retain the same structure as for the second
      ## case
      out2 <- data.frame(
        "Stage"            = NA,
        "Distance_measure" = NA,
        "99%"              = NA,
        "95%"              = NA,
        "90%"              = NA,
        check.names = FALSE
      )
    }
    
    # out3 <- data.frame(
    #   "Status_code"   = names(verify(.object)[[1]]),
    #   "Status_code_1" = c(verify(.object)[[1]]),
    #   "Status_code_2" = c(verify(.object)[[2]])
    # )
  } else {
    
    selector_weights2 <- c("eta3 ~ c1", "eta3 ~ c2", "eta3 ~ c3")
    
    sum_res      <- cSEM:::summarize(.object)
    path         <- sum_res$Estimates$Path_estimates[
      sum_res$Estimates$Path_estimates$Name %in% selector_path, relevant, drop = FALSE]
    weights      <- sum_res$Estimates$Weight_estimates[
      sum_res$Estimates$Weight_estimates$Name %in% c(selector_weights_1stage, selector_weights_2ndstage), relevant, drop = FALSE]
    weights2     <- sum_res$Estimates$Path_estimates[
      sum_res$Estimates$Path_estimates$Name %in% selector_weights2, relevant, drop = FALSE]
    
    if(nrow(weights2) == 3) {
      weights2$Name <- c("eta3 <~ c1", "eta3 <~ c2", "eta3 <~ c3")
      weights <- dplyr::bind_rows(weights, weights2)
    }
    total_effect <- sum_res$Estimates$Effect_estimates$Total_effect[
      sum_res$Estimates$Effect_estimates$Total_effect$Name %in% selector_te, relevant, drop = FALSE]
    
    out1 <- dplyr::bind_rows(
      list("Path" = path, 
           "Weights" = weights,
           "Total_effect" = total_effect), 
      .id = "Estimate_type")

    ## Only compute test if not repeated indicators approach (in this case the 
    ## length of data would be 37)
    if(ncol(.object$Information$Data) == 23) {
      ## Joint test
      out2a <- testOMF(.object, 
                       .handle_inadmissibles = "replace",
                       .alpha = c(0.1, 0.05, 0.01),
                       .verbose = FALSE, .R = .R)
      
      out2b <- data.frame(
        "Distance_measure" = names(out2a$Test_statistic),
        stringsAsFactors = FALSE
      )
      
      out2a <- cbind(out2b, out2a$Decision)
      rownames(out2a) <- NULL
      out2 <- cbind("Stage" = "Both", out2a) 
      
    } else {
      ## Set up empty data set to retain the same structure as for the second-order
      ## case
      out2 <- data.frame(
        "Stage"            = NA,
        "Distance_measure" = NA,
        "99%"              = NA,
        "95%"              = NA,
        "90%"              = NA,
        check.names = FALSE
      )      
    }
  }
  
  ## Return list 
  out <- list(
    "Parameter_estimates" = out1, 
    "Test"                = out2
    # "Status_code"         = out3
  )
}

checkGuidlines <- function(.data, .model, .approach, .what) {
  
  # Estimate
  if(.approach == "emb_TS") {
    out <- csem(.data, .model, .approach_2ndorder = "mixed")
  } else {
    out <- csem(.data, .model)
  }

  if(.what %in% c("Weight_estimates", "Loading_estimates")) {
    sw <- switch (.approach,
           "TS"     = {cSEM::summarize(out)$First_stage$Estimates[[.what]]},
           "emb_TS" = {cSEM::summarize(out)$First_stage$Estimates[[.what]][1:21, ]},
           "RI"     = {cSEM::summarize(out)$Estimates[[.what]][c(1:3, 10:21, 4:9), ]}
    )
    
    sw$group <- c(rep("xi", 3), rep("c1", 2), rep("c2", 4), rep("c3", 6), rep("eta1", 3), rep("eta2", 3))
    out1 <- round(unlist(lapply(split(sw, sw$group), function(x) min(x$Estimate))), 2)[c("xi", "eta1", "eta2", "c1", "c2", "c3")]
    
    if(.approach %in% c("TS", "emb_TS")) {
      out2 <- round(min(cSEM::summarize(out)$Second_stage$Estimates[[.what]]$Estimate), 2)
    } else {
      if(.what == "Weight_estimates") {
        out2 <- round(min(cSEM::summarize(out)$Estimates$Path_estimates[15:17, ]$Estimate), 2)  
      } else {
        out2 <- round(min(cor(out$Estimates$Construct_scores)[4:6, 7]), 2)
      }
    }
  } else if(.what == "VIF") { # VIF
    if(.approach %in% c("TS", "emb_TS")) {
      out1 <- round(unlist(lapply(cSEM:::calculateVIFModeB(out$First_stage), max))[c("xi", "eta1", "eta2", "c1", "c2", "c3")], 1)
      out2 <- round(unlist(lapply(cSEM:::calculateVIFModeB(out$Second_stage), max))["eta3"], 1)
    } else {
      out1 <- round(unlist(lapply(cSEM:::calculateVIFModeB(out), max))[c("xi", "eta1", "eta2", "c1", "c2", "c3")], 1)
      out2 <- round(max(out$Estimates$VIF$eta3), 1)
    }
  }
  names(out2) <- "eta3"
  c(out1, out2)
}