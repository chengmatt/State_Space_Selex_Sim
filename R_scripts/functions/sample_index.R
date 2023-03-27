# Purpose: To sample indices of abundance from the OM
# This function is meant to be used within the OM, when all objects are 
# loaded into the environment
# Creator: Matthew LH. Cheng
# Date: 11/3/22

# Generate Obs Error ------------------------------------------------------

#' @param error Type of error we want to generate - log normal or normal
#' @param true_index Input index observed with no error
#' @param CV Amount of CV we want to apply onto the index

idx_obs_error <- function(error, true_index, CV) {

  if(error == "log_normal") {

    # Convert CV to sd for log normal - inner fxn = variance
    sd <- sqrt(log((CV^2) + 1))
    
    # Sample!
    sampled_index <- true_index * exp(rnorm(1, -sd^2/2, sd))
    
  } # if generating via log normal
  
  return(sampled_index)
  
}
