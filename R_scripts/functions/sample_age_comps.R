# Purpose: To sample age comps from the OM
# This function is meant to be used within the OM, when all objects are 
# loaded into the environment
# Creator: Matthew LH. Cheng
# Date: 11/3/22

#' @param Comp_Fleet Fleet type - options are Fishery or Survey
#' @param error Error structure for sampling age comps - options currently are multinomial

sample_age_comps <- function(Comp_Fleet, error) {
  
  if(Comp_Fleet == "Fishery") {
    
    # Get proportions at a given age using the expected CPUE at age for the index in numbers
    Prob_Catch_at_age <- (N_at_age[y-1,,sim] * Fish_selex_at_age[y-1,,sim]) /sum(N_at_age[y-1,,sim] * Fish_selex_at_age[y-1,,sim])
    
    if(error == "multinomial") {
      # Generate comps based on the expected catch at age
      Age_Comps <- rmultinom(n = 1, size = N_eff_fish, prob = Prob_Catch_at_age)[,1]
    } # if our error is multinomial
    
  } # if our fleet is fishery
  
  if(Comp_Fleet == "Survey") {
    
    # Calculate selected individuals based on the expected CPUE at age
    Prob_Surv_at_age <- N_at_age[y-1,,sim] * Surv_selex_at_age[y-1,,sim] /  sum(N_at_age[y-1,,sim] * Surv_selex_at_age[y-1,,sim])
    
    if(error == "multinomial") {
      # Generate comps based on the expected CPUE at agew
      Age_Comps <- rmultinom(n = 1, size = N_eff_surv, prob = Prob_Surv_at_age)[,1]
    } # if our error structure is multinomial
    
  } # if our fleet is survey 

  return(Age_Comps)

} # end function
