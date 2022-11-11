# Purpose: To specify catchability for the fishery and the survey in the OM (only time-invariant options for now)
# Creator: Matthew LH. Cheng
# Date: 10/30/22

#' @param q_Mean_Fish Mean catchability parameter for the fishery
#' @param q_Mean_Surv Mean catchability parameter for the survey

specify_q <- function(q_Mean_Fish, q_Mean_Surv) {
  
  if(n_fish_fleets != length(q_Mean_Fish)) stop("Length of q vector for fishery does not equal the number of fishery fleets specified!")
  if(n_srv_fleets != length(q_Mean_Surv)) stop("Length of q vector for syrvey does not equal the number of survey fleets specified!")
  
  for(f in 1:n_fish_fleets) {
    q_Fish[,f,] <- q_Mean_Fish[f]
  } # loop through to input q_Mean_Fish into this
  
  for(sf in 1:n_srv_fleets) {
    q_Surv[,sf,] <- q_Mean_Surv[sf]
  } # loop through for survey fleets
  
  # Output these to the environment
  q_Fish <<- q_Fish
  q_Surv <<- q_Surv
  
} # end function
