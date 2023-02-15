# Purpose: To predict new recruits following a beverton holt stock
# recruit relationship
# Creator: Matthew LH. Cheng
# Date: 10/29/22

#' @param ssb Spawning stock biomass 
#' @param h steepness (how steep you want the BH curve to be)
#' @param r0 Virgin recruitment (where recruitment asymptotes at)

beverton_holt_recruit <- function(ssb, h, r0) {
  
  # Get SPR
  SPR_SSB0 <- vector()
  
  for(a in 1:length(ages)) {
    SPR_SSB0[a] = exp(-Mort_at_age[1,a,1] * (a - 1)) * wt_at_age[1, a, 1, 1] * mat_at_age[1, a, 1, 1]
  }
  
  # Now, get SPR rate
  SPR0 <- sum(SPR_SSB0)
  ssb0 <- SPR0 * r0
  
  # Output to environemnt
  SPR_SSB0 <<- SPR_SSB0
  
  # Get alpha and beta for BH
  alpha <- ((1 - h) * ssb0 )/ (4 * h * r0 )
  beta <- ((5 * h - 1) ) / (4 * h * r0)
  
  # Now calculate BH parameterization
  rec <- ssb / (alpha + (beta * ssb))

  return( rec )
}
