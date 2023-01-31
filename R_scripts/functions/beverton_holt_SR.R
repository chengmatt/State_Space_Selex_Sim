# Purpose: To predict new recruits following a beverton holt stock
# recruit relationship
# Creator: Matthew LH. Cheng
# Date: 10/29/22

#' @param ssb Spawning stock biomass 
#' @param h steepness (how steep you want the BH curve to be)
#' @param r0 Virgin recruitment (where recruitment asymptotes at)
#' @param ssb0 Spawning stock biomass at unfished biomass

beverton_holt_recruit <- function(ssb, h, r0, ssb0) {
  return( (ssb*4*h*r0)/(ssb0*(1-h) + ssb*(5*h-1)) )
}

