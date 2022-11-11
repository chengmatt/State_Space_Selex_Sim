# Purpose: To sample age comps from the OM
# This function is meant to be used within the OM, when all objects are 
# loaded into the environment
# Creator: Matthew LH. Cheng
# Date: 11/3/22

#' @param Comp_Fleet Fleet type - options are Fishery or Survey
#' @param error Error structure for sampling age comps - options currently are multinomial
#' @param N_eff Effective sample size for the compositional dataset, corresponding to a fleet

sample_comps <- function(Comp_Fleet, error, N_eff) {
  
  if(Comp_Fleet == "Fishery") {
    
    # Get proportions at a given age using the expected CPUE at age for the index in numbers - normalize to sum to 1
    Prob_Catch_at_age <- (N_at_age[y-1,,s,sim] * Fish_selex_at_age[y-1,,f,s,sim]) / sum(N_at_age[y-1,,s,sim] * Fish_selex_at_age[y-1,,f,s,sim])
    
    if(error == "multinomial") {
      # Generate comps based on the expected catch at age
      Age_Comps <- rmultinom(n = 1, size = N_eff, prob = Prob_Catch_at_age)[,1]
    } # if our error is multinomial
    
  } # if our fleet is fishery
  
  if(Comp_Fleet == "Survey") {
    
    # Calculate selected individuals based on the expected CPUE at age
    Prob_Surv_at_age <- N_at_age[y-1,,s,sim] * Surv_selex_at_age[y-1,,sf,s,sim] / sum(N_at_age[y-1,,s,sim] * Surv_selex_at_age[y-1,,sf,s,sim])
    
    if(error == "multinomial") {
      # Generate comps based on the expected CPUE at agew
      Age_Comps <- rmultinom(n = 1, size = N_eff, prob = Prob_Surv_at_age)[,1]
    } # if our error structure is multinomial
    
  } # if our fleet is survey 

  return(Age_Comps)

} # end function
