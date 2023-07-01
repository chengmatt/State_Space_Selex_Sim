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

#' Title get_trialF_spr
#'
#' @param MortAA Vector of mortality at age
#' @param SelexAA Vector of selectivity at age
#' @param MatAA Vector of maturity at age
#' @param WAA Vector of weight at age
#' @param trial_F Vector of trial F valeus
#' @param F_x F_x SPR rate
#'
#' @return
#' @export
#'
#' @examples
get_trialF_spr <- function(MortAA, 
                           SelexAA,
                           MatAA,
                           WAA,
                           trial_F,
                           F_x = 0.4
                           ) {
  # Pre-Processing
  N_Init <- 1 # Starting Numbers 

  # Get Unfished SSB
  Z <- c(0, MortAA[-length(MortAA)]) # Get Natural mortality
  N <- N_Init * exp(-cumsum(Z)) # Calculate Numbers over the lifespan of cohort
  SSB_Unfished <- sum(N * WAA * MatAA) # Get SSB in biomass units
  
  # Get SPR rates
  spr_values <- vector()
  spr_diff <- vector()
  
  # Loop through to look at SPR across Fs
  for(i in 1:length(trial_F)) {
    
    # Get FAA
    FAA = trial_F[i] * SelexAA
    Z <- c(0, FAA[-length(FAA)] + MortAA[-length(MortAA)]) # Get Total Mortality w/o + group
    N <- N_Init * exp(-cumsum(Z)) # Calculate Numbers over the lifespan of cohort
    SSB_Fished <- sum(N * WAA * MatAA) # Get SSB in biomass units
    
    # F rate that results in x% of fished/unfished SSB
    spr_values[i] <- abs((SSB_Fished / SSB_Unfished)) 
    spr_diff[i] <- abs((SSB_Fished / SSB_Unfished) - F_x)  
  } # end i loop

  return(data.frame(SPR = spr_values, Diff = spr_diff, trial_F = trial_F))

} # end function


#' Title Get ABC reference points using SPR Fx%
#'
#' @param Fx_proj Fx% reference point from SPR calculations
#' @param n_ages Number of ages
#' @param n_sex Number of sexes
#' @param n_fleets Number of fishery fleets
#' @param mean_rec Mean recruitment
#' @param term_NAA Array of terminal year NAA
#' @param term_F_Slx Array of terminal Fishery selex
#' @param term_F Vector of terminal year Fs
#' @param MortAA Array of natural mortality by age
#' @param WAA Array of weight at age values
#'
#' @return Numeric value of ABC
#' @export
#'
#' @examples
get_ABC_refpt <- function(Fx_proj,
                          terminal_yr = 50,
                          sex_ratio = 0.5,
                          n_ages = 30,
                          n_sex = 2,
                          n_fleets = 2,
                          mean_rec = 10,
                          term_NAA,
                          term_F_Slx,
                          term_F,
                          MortAA,
                          WAA) {
  
  # TESTING
  # term_NAA = oms$N_at_age[50,,,1]
  # term_F = oms$fish_mort[50,2,1]
  # MortAA = oms$Mort_at_age[50,,1]
  # term_F_Slx = array(oms$Fish_selex_at_age[50,,2,,1],
  #                    dim = c(30,1,2))
  # Fx_proj = 0.35
  # WAA = oms$wt_at_age[50,,,1]
  
  # Empty matrices to store projections in
  N_proj <- array(data = 0, dim = c(n_ages, n_sex))
  F40_proj <- array(data = 0, dim = c(n_ages, n_sex))
  Z_proj <- array(data = 0, dim = c(n_ages, n_sex))
  CAA_proj <- array(data = 0, dim = c(n_ages, n_sex))
  Catch_proj <- array(data = 0, dim = c(n_ages, n_sex))
  
  # Multiply mean recruitment by sex ratio
  mean_rec = mean_rec * sex_ratio
  # Input into N_proj
  N_proj[1,] = mean_rec
  
  for(s in 1:n_sex) {
    for(a in 1:n_ages) {
      # Get terminal total mortality here (constant M)
      term_ZAA = sum(term_F * term_F_Slx[a,,s], na.rm = TRUE) + MortAA[a]
      if(a < n_ages) { # not plus group
        N_proj[a + 1, s] = term_NAA[a, s] * exp(-term_ZAA)
      } else{ 
        N_proj[a, s] = N_proj[a, s] + (term_NAA[a, s] * exp(-term_ZAA))
      } # else = plus group
      
      # Get F40 projections here

      # get F40 projections using fratio
      fratio <- term_F / sum(term_F) # Calculate Fratio here
      for(f in 1:n_fleets) {
        F40_proj[a, s] = F40_proj[a, s] + (fratio[f] * Fx_proj * term_F_Slx[a, f, s])
      } # end f (fleets) loop

      # Get mortality and catch projections now
      Z_proj[a,s] <- F40_proj[a,s] + MortAA[a]
      # Now, get catch projections
      CAA_proj[a,s] <- N_proj[a,s] * (1 - exp(-Z_proj[a,s])) * (F40_proj[a,s] / Z_proj[a,s]) 
      Catch_proj[a,s] <- CAA_proj[a,s] * WAA[a,s] # Turn to biomass units
      ABC <- sum(Catch_proj) # sum to get abc
      
    } # end a (age) loop
  } # end s (sex) loop
  
  return(ABC)
  
} # end function