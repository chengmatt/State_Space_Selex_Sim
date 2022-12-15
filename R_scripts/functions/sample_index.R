# Purpose: To sample indices of abundance from the OM
# This function is meant to be used within the OM, when all objects are 
# loaded into the environment
# Creator: Matthew LH. Cheng
# Date: 11/3/22


# Sample Index ------------------------------------------------------------

#' @param Idx_Fleet Type of fleet we are sampling for: options are Fishery or Survey

sample_index <- function(Idx_Fleet) {
  
  if(Idx_Fleet == "Fishery") {
    
    # Calculate selected individuals in numbers
    true_N <- sum(N_at_age[y-1,,s,sim] * wt_at_age[y-1,,s,sim] * Fish_selex_at_age[y-1,,f,s,sim])

    # Now, calculate the true index
    true_index <- true_N * q_Fish[y-1,f,sim] 

  } # if we are sampling from the fishery
  
  if(Idx_Fleet == "Survey") {
    
    # Calculate selected individuals in biomass
    true_N <- sum(N_at_age[y-1,,s,sim] * Surv_selex_at_age[y-1,,sf,s,sim])
    
    # Now, calculate the true index
    true_index <- true_N * q_Surv[y-1,sf,sim] 
    
  } # if we are sampling from the survey

  return(true_index) # output this out

} # end function



# Generate Obs Error ------------------------------------------------------

#' @param error Type of error we want to generate - log normal or normal
#' @param true_index Input index observed with no error
#' @param CV Amount of CV we want to apply onto the index

idx_obs_error <- function(error, true_index, CV) {
  
  if(error == "normal") {

    # Convert CV to SD
    sd <- CV * true_index

    # Sample!
    sampled_index <- rnorm(n=1, mean=true_index, sd=sd)

  } # if generating via normal

  if(error == "log_normal") {

    # Convert CV to sd for log normal
    sd <- log((CV^2) + 1)

    # Sample!
    sampled_index <- true_index * exp(rnorm(1, 0, sd))

  } # if generating via log normal
  
  return(sampled_index)
  
}
