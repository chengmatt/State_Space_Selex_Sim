# Purpose: To sample age comps from the OM
# This function is meant to be used within the OM, when all objects are 
# loaded into the environment
# Creator: Matthew LH. Cheng
# Date: 11/3/22

#' @param Comp_Fleet Fleet type - options are Fishery or Survey
#' @param error Error structure for sampling age comps - options currently are multinomial
#' @param N_eff Effective sample size for the compositional dataset, corresponding to a fleet

sample_comps <- function(Comp_Fleet, error, N_eff, prob) {
  if(Comp_Fleet == "Fishery") {
    if(error == "multinomial") {
      # Generate comps based on the expected catch at age
      Age_Comps <- rmultinom(n = 1, size = N_eff, prob = prob)[,1]
    } # if our error is multinomial
    
  } # if our fleet is fishery
  
  if(Comp_Fleet == "Survey") {
    if(error == "multinomial") {
      # Generate comps based on the expected CPUE at agew
      Age_Comps <- rmultinom(n = 1, size = N_eff, prob = prob)[,1]
    } # if our error structure is multinomial
    
  } # if our fleet is survey 

  return(Age_Comps)

} # end function
