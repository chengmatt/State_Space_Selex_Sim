# Purpose: To specify catchability for the fishery and the survey in the OM
# Creator: Matthew LH. Cheng
# Date: 10/30/22


#' @param q_Fish_Time Time parameterization of catchability for fishery
#' @param q_Surv_Time Time parameterization of catchability for survey
#' @param q_Mean_Fish Mean catchability parameter for the fishery
#' @param q_Mean_Surv Mean catchability parameter for the survey
#' @param q_Fish_rho Correaltion for AR1 in q fish
#' @param q_Fish_sigma Standard deviation for AR1 in q fish

specify_q <- function(q_Fish_Time = "Constant", q_Surv_Time = "Constant",
                             q_Mean_Fish, q_Mean_Surv, 
                             q_Fish_rho, q_Fish_sigma) {
  
  if(q_Fish_Time == "Constant") {
    # Input value into array
    q_Fish[,] <- q_Mean_Fish
  } # time invariant catchability
  
  
  if(q_Surv_Time == "Constant") {
    # Input value into array
    q_Surv[,] <- q_Mean_Surv
  } # time invariant catchability
  
  # Output these to the environment
  q_Fish <<- q_Fish
  q_Surv <<- q_Surv
  
  print(paste("Catchability for the Fishery is specified as:", q_Fish_Time))
  print(paste("Catchability for the Survey is specified as:", q_Surv_Time))
  
} # end function
