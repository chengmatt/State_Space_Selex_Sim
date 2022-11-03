# Purpose: To sample indices of abundance from the OM
# This function is meant to be used within the OM, when all objects are 
# loaded into the environment
# Creator: Matthew LH. Cheng
# Date: 11/3/22


#' @param Idx_Fleet Type of fleet we are sampling for: options are Fishery or Survey
#' @param error Error structure we want to generate our data from: options are  normal

sample_index <- function(Idx_Fleet, error = "normal") {
  
  if(Idx_Fleet == "Fishery") {
    
    # Calculate selected individuals in abundance
    selected_indv_N <- N_at_age[y-1,,sim] * Fish_selex_at_age[y-1,,sim] 
    
    # Translate that to biomass at age
    selected_indv_B <- selected_indv_N * wt_at_age[y-1,,sim]
    
    # sum across to get total biomass
    true_B <- sum(selected_indv_B)
    
    # Now, calculate the true index
    true_index <- true_B * q_Fish[y-1,sim] 
    
    # Convert CV to SD
    sd <- Fishery_CV * true_index
    
    if(error == "normal") {
      sampled_index <- rnorm(n=1, mean=true_index, sd=sd)
    } # if generating via normal

  } # if we are sampling from the fishery
  
  if(Idx_Fleet == "Survey") {
    
    # Calculate selected individuals in abundance
    selected_indv_N <- N_at_age[y-1,,sim] * Surv_selex_at_age[y-1,,sim] 
    
    # Translate that to biomass at age
    selected_indv_B <- selected_indv_N * wt_at_age[y-1,,sim]
    
    # sum across to get total biomass
    true_B <- sum(selected_indv_B)
    
    # Now, calculate the true index
    true_index <- true_B * q_Surv[y-1,sim] 
    
    # Convert CV to SD
    sd <- Survey_CV * true_index
    
    if(error == "normal") {
      sampled_index <- rnorm(n=1, mean=true_index, sd=sd)
    } # if generating via normal
    
  } # if we are sampling from the survey

  return(sampled_index) # output this out

} # end function
