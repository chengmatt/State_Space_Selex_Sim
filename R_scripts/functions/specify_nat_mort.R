# Purpose: To specify natural mortality for the OM
# Creator: Matthew LH. Cheng
# Date: 10/30/22



#' @param Mort_Time Time parameterization for natural mortality
#' @param Mean_M Mean value of natural mortality

specify_nat_mort <- function(Mort_Time = "Constant", Mean_M) {
  
  if(Mort_Time == "Constant") {
    
    Mort_at_age[,,] <- Mean_M
    
  } # if this is time invariant natural mortality
  
  # Output this into our environment
  Mort_at_age <<- Mort_at_age
  
  print(paste("Natural Mortality is specified as:", Mort_Time))
  
} # end function
