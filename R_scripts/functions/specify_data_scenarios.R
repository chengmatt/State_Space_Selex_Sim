# Purpose: To specify data quality and quantity scenarios for the fishery
# and the survey, start year for the fishery, and start year for the survey
# Creator: Matthew LH. Cheng
# Date: 10/30/22

#' @param Fish_Start_yr Fishery start year
#' @param Surv_Start_yr Survey start year
#' @param fish_Neff_max Fishery Effective Sample sizes
#' @param srv_Neff_max Survey effective sample sizes
#' @param fish_CV Fishery CVs for indices
#' @param srv_CV Survey CVs for indices
#' @param Neff_Fish_Time Whether Neff for the fisheries remain constant, or vary as a function of time (Constant or F_Vary)
#' @param fish_mort Input fishing mortality array if we want to allow fish_Neff to vary as a function of that
#' @param fixed_Neff Value for which Neff is fixed at - user-specified so that F_vary's N_eff determination doesn't
#' drop below a certain value. 

fish_surv_data_scenarios <- function(Fish_Start_yr, Surv_Start_yr, fish_Neff_max, srv_Neff_max,
                                     fish_CV, srv_CV, Neff_Fish_Time = "Constant", fish_mort = NULL,
                                     fixed_Neff) {
  
  # Warning/Error Messages
  if(n_fish_fleets != length(Fish_Start_yr) & 
     n_fish_fleets != length(fish_Neff) &
     n_fish_fleets != length(fish_CV) ) stop("Vector lengths of fishery start years, Neff, or CV are not the same as the number of fishery fleets specified")
  
  if(n_srv_fleets != length(Surv_Start_yr) &
     n_srv_fleets != length(srv_Neff) &
     n_srv_fleets != length(srv_CV) ) stop("Vector lengths of survey start years, Neff, or CV are not the same as the number of fishery fleets specified")
  
  if(!Neff_Fish_Time %in% c("Constant", "F_Vary")) stop("Neff_Fish_Time is not specified correctly")

  # Create objects to output into the environment for use
  Fish_Start_yr <<- Fish_Start_yr # Fishery start year
  Surv_Start_yr <<- Surv_Start_yr # Survey Start Year
  fish_CV <<- fish_CV # Fishery CV
  srv_CV <<- srv_CV # Survey CV
  
  if(Neff_Fish_Time == "Constant") { # If our effective sample sizes stay constant throughout the fishery
    fish_Neff_mat <- matrix(data = fish_Neff, nrow = n_years, ncol = n_fish_fleets, byrow = TRUE) # create matrix
  } # if constant
  
  if(Neff_Fish_Time == "F_Vary") {
    
    # Basically, we are going to calculate a scalar based on the desired effective sample
    # size at the peak for a given fishery
    F_scalar <- vector() # Create an empty vector to store the F scalar

    # Create empty matrix
    fish_Neff_mat <- matrix(data = NA, nrow = n_years, ncol = n_fish_fleets, byrow = TRUE)
  
    for(f in 1:n_fish_fleets) {
      
      # Caclulate scalar such that the max fish mort results in the desired Neff
      F_scalar[f] <- fish_Neff_max[f] / (fish_Neff_max[f] * max(fish_mort[,f,])) 
      
      for(y in 1:n_years) {

          # Scale N_effs by relative fishing mortality rates - fix sim index at 1,
          # because these will remain constant across simulations
          fish_Neff_mat[y,f] <- round(fish_Neff_max[f] * F_scalar[f] * fish_mort[y,f,1])

      } # end fleet loop
      
    } # end year loop
    
    # Fix Neff at a specified value for those that are < said specified value
    fish_Neff_mat[fish_Neff_mat < fixed_Neff & fish_Neff_mat != 0] <- fixed_Neff
    
  } # effective sample sizes vary as a function of fishery specific Fs
  
  # create empty matrix for survey Neff
  srv_Neff_mat <- matrix(data = NA, nrow = n_years, ncol = n_srv_fleets, byrow = TRUE)
  
  # specify matrix for survey neff
  for(sf in 1:n_srv_fleets) {
    
    for(y in 1:n_years) {
      
      if(y < Surv_Start_yr[sf]) {
        srv_Neff_mat[y,sf] <- 0 # force neff to be 0 when fishing is not occuring
      } else{
        srv_Neff_mat[y,sf] <- srv_Neff_max[sf]
      }
    }
    
  } # end sf loop
  
  # output into environment
  fish_Neff <<- fish_Neff_mat
  srv_Neff <<- srv_Neff_mat
  
} # end function
