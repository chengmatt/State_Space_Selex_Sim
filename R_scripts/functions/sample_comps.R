# Purpose: To sample age comps from the OM
# This function is meant to be used within the OM, when all objects are 
# loaded into the environment
# Creator: Matthew LH. Cheng
# Date: 11/3/22

#' @param alpha Dirichlet Alpha parameter
#' @param Input_N Effective sample size for the compositional dataset, corresponding to a fleet
rdirmultinom <- function(Input_N, alpha){ 
  # Make independent draws for a gamma -> Dirichlet draws
  x <- rgamma(length(alpha), alpha)
  dir_draws <- x/sum(x) # normalize proportions
  # Sample from multinomial conditional on dirichlet
  rdirmultinom_probs <- as.vector(stats::rmultinom(1, Input_N, dir_draws))
  return(rdirmultinom_probs)
}


#' @param error Error structure for sampling age comps - options currently are multinomial
#' @param Input_N Effective sample size for the compositional dataset, corresponding to a fleet
#' @param DM_Param Weight/Theta parameter governing the variability of the Dirichlet-Multinomial
#' @param prob Underlying probability (vector) of sampling

sample_comps <- function(error, Input_N, DM_Param = NULL, prob) {
  
    if(error == "multinomial") {
      # Generate comps based on the expected catch at age
      Age_Comps <- rmultinom(n = 1, size = Input_N, prob = prob)
    } # if our error is multinomial
  
  if(error == "dirichlet") {
    lambda <- Input_N/DM_Param^2 - 1
    alpha <- prob * lambda # Get concentration parameter
    Comp_Probs <- rgamma(length(alpha), alpha) # Sample from dirichlet
    Age_Comps <- Comp_Probs/sum(Comp_Probs) # Normalize!
  }
  
  if(error == "dirichlet_multinomial") { # Thorson et al. 2017
    # The Dirichlet portion of this essentially modifies the inherent sampling probability of CAA
    # Dir_Probs <- gtools::rdirichlet(n = length(prob), alpha = prob * DM_Param)
    alpha <- prob * DM_Param * Input_N # Get dirichlet alpha parameter
    # Age_Comps <- rdirmultinom(Input_N = Input_N, alpha = alpha) # Get dirichlet multinomial samples
    Age_Comps <- extraDistr::rdirmnom(n = 1, size = Input_N, alpha = alpha)
  } # if our error is dirichlet multinomial

  return(Age_Comps)

} # end function
