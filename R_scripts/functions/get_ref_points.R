# Purpose: Compute Reference Points from OM
# Creator: Matthew LH. Cheng
# Date 3/22/23

#' Title Function to compute Fx% via optimization
#'
#' @param Trial_F Trial "parameter" values to be fed into an optimizer
#' @param F_x SPR% we want to acheive
#' @param fratio Ratio of Fs if > 1 fleet
#' @param SelexAA Selectivity at Age
#' @param WAA Weight at Age
#' @param MatAA Maturity at Age
#' @param MortAA Mortality at Age
#' @param N_Init Inintial Abundance for SPR (doesn't really matter)
#' @param SSB_Unfished Unfished SSB at equilibrium
#'
#' @return
#' @export
#'
#' @examples
Fx_optim <- function(Trial_F,
                   F_x = F_x, 
                   fratio = fratio,
                   SelexAA = SelexAA,
                   WAA = WAA,
                   MatAA = MatAA,
                   MortAA = MortAA,
                   N_Init, 
                   SSB_Unfished) {
  
  # One fleet
  if(length(fratio) == 1)  FAA <- Trial_F * (fratio * SelexAA)  # Calculate fishing mortality
  # Multiple fleets
  if(length(fratio) > 1)  FAA <- Trial_F * colSums((fratio * SelexAA))  # Calculate fishing mortality

  Z <- c(0, FAA[-length(FAA)] + MortAA[-length(MortAA)]) # Get Total Mortality w/o + group
  N <- N_Init * exp(-cumsum(Z)) # Calculate Numbers over the lifespan of cohort
  SSB_Fished <- sum(N * WAA * MatAA) # Get SSB in biomass units
  
  # Evaluating our objective function here 
  # F rate that results in x% of fished/unfished SSB
  diff <- abs((SSB_Fished / SSB_Unfished) - F_x)  
  return(diff)
  
} # end get F_x

#' Title Get F-Based spawning Biomass per recruit reference points
#'
#' @param ages vector of ages
#' @param MortAA Natural Mortality at Age Matrix
#' @param SelexAA Selectivity at Age Matrix fleet x ages
#' @param MatAA Maturity at Age Matrix
#' @param WAA Weight at age Matrix
#' @param n_fleets Number of Fleets
#' @param Terminal_F Terminal Fishing Mortality Rates
#' @param F_x F rate that reduces Fished SSB/Unfished SSB to x%
#'
#' @return
#' @export
#'
#' @examples
get_Fx_refpt <- function(ages, 
                           MortAA, 
                           SelexAA,
                           MatAA,
                           WAA,
                           Terminal_F,
                           F_x = 0.4
                           ) {
  # Pre-Processing
  N_Init <- 1 # Starting Numbers 
  # Get ratio of Fs - our assumption of a given "fishing regime"
  fratio <- Terminal_F/sum(Terminal_F)

  # Get Unfished SSB
  Z <- c(0, MortAA[-length(MortAA)]) # Get Natural mortality
  N <- N_Init * exp(-cumsum(Z)) # Calculate Numbers over the lifespan of cohort
  SSB_Unfished <- sum(N * WAA * MatAA) # Get SSB in biomass units
  
  # Optimize to get F that would result in SBPR being x% of unfished
  F_x_optim <- bbmle::mle2(Fx_optim, 
             data = list(fratio = fratio, 
                         SelexAA = SelexAA,
                         WAA = WAA, 
                         MatAA = MatAA,
                         MortAA = MortAA, 
                         N_Init = N_Init, 
                         SSB_Unfished = SSB_Unfished,
                         F_x = F_x), 
             start = list(Trial_F = 0.01), 
             optimizer = "nlminb",
             method="Nelder-Mead",
             control = list(maxit = 1e3))
  
  # Get estimated Fx% out
  F_x_val <- F_x_optim@coef

  return(F_x_val)

} # end function
