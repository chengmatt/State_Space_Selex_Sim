# Purpose: To sample indices of abundance from the OM
# This function is meant to be used within the OM, when all objects are 
# loaded into the environment
# Creator: Matthew LH. Cheng
# Date: 11/3/22


#' @param Idx_Fleet Type of fleet we are sampling for: options are Fishery or Survey
#' @param error Error structure we want to generate our data from: options are  normal

sample_index <- function(Idx_Fleet, error = "normal") {
  
  if(Idx_Fleet == "Fishery") {
    
    # Calculate selected individuals in numbers
    true_N <- sum(N_at_age[y-1,,sim] * Fish_selex_at_age[y-1,,sim])

    # Now, calculate the true index
    true_index <- true_N * q_Fish[y-1,sim] 
    
    if(error == "normal") {
      
      # Convert CV to SD
      sd <- Fishery_CV * true_index
      
      sampled_index <- rnorm(n=1, mean=true_index, sd=sd)
    } # if generating via normal
    
    if(error == "log_normal") {
      
      # Convert CV to sd for log normal
      sd <- log((Fishery_CV^2) + 1)
      
      # Sample!
      sampled_index <- true_index * exp(rnorm(1, 0, sd))
      
    } # if generating via log normal

  } # if we are sampling from the fishery
  
  if(Idx_Fleet == "Survey") {
    
    # Calculate selected individuals in biomass
    true_N <- sum(N_at_age[y-1,,sim] * Surv_selex_at_age[y-1,,sim])
    
    # Now, calculate the true index
    true_index <- true_N * q_Surv[y-1,sim] 
    
    if(error == "normal") {
      
      # Convert CV to SD
      sd <- Survey_CV * true_index
      
      # Sample!
      sampled_index <- rnorm(n=1, mean=true_index, sd=sd)
      
    } # if generating via normal
    
    if(error == "log_normal") {
      
      # Convert CV to sd for log normal
      sd <- log((Survey_CV^2) + 1)
      
      # Sample!
      sampled_index <- true_index * exp(rnorm(1, 0, sd))
      
    } # if generating via log normal
    
  } # if we are sampling from the survey

  return(sampled_index) # output this out

} # end function
