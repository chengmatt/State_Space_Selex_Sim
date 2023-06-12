# Purpose: To specify data quality and quantity scenarios for the fishery
# and the survey, start year for the fishery, and start year for the survey
# Creator: Matthew LH. Cheng
# Date: 10/30/22

#' @param Fish_Start_yr Fishery start year
#' @param Surv_Start_yr Survey start year
#' @param Input_Fish_N_Max Fishery Effective Sample sizes
#' @param Input_Srv_N_Max Survey effective sample sizes
#' @param fish_CV Fishery CVs for indices
#' @param srv_CV Survey CVs for indices
#' @param Input_N_Fish_Time Whether Neff for the fisheries remain constant, or vary as a function of time (Constant or F_Vary)
#' @param fish_mort Input fishing mortality array if we want to allow fish_Neff to vary as a function of that
#' @param Input_N_Fish_Fixed Value for which Neff is fixed at - user-specified so that F_vary's N_eff determination doesn't
#' drop below a certain value. 

fish_surv_data_scenarios <- function(Fish_Start_yr, Surv_Start_yr, Input_Fish_N_Max, Input_Srv_N_Max,
                                     fish_CV, srv_CV, Input_N_Fish_Time = "Constant", fish_mort = NULL,
                                     Input_N_Fish_Fixed) {
  
  # Warning/Error Messages
  if(n_fish_fleets != length(Fish_Start_yr) & 
     n_fish_fleets != length(Input_Fish_N_Max) &
     n_fish_fleets != length(fish_CV) ) stop("Vector lengths of fishery start years, Neff, or CV are not the same as the number of fishery fleets specified")
  
  if(n_srv_fleets != length(Surv_Start_yr) &
     n_srv_fleets != length(Input_Srv_N_Max) &
     n_srv_fleets != length(srv_CV) ) stop("Vector lengths of survey start years, Neff, or CV are not the same as the number of fishery fleets specified")
  
  if(!Input_N_Fish_Time %in% c("Constant", "F_Vary")) stop("Input_N_Fish_Time is not specified correctly")

  # Create objects to output into the environment for use
  Fish_Start_yr <<- Fish_Start_yr # Fishery start year
  Surv_Start_yr <<- Surv_Start_yr # Survey Start Year
  fish_CV <<- fish_CV # Fishery CV
  srv_CV <<- srv_CV # Survey CV
  
  if(Input_N_Fish_Time == "Constant") { # If our effective sample sizes stay constant throughout the fishery
    fish_Input_N_mat <- matrix(data = Input_Fish_N_Max, nrow = n_years, ncol = n_fish_fleets, byrow = TRUE) # create matrix
  } # if constant
  
  if(Input_N_Fish_Time == "F_Vary") {
    
    # Basically, we are going to calculate a scalar based on the desired effective sample
    # size at the peak for a given fishery
    F_scalar <- vector() # Create an empty vector to store the F scalar

    # Create empty matrix
    fish_Input_N_mat <- matrix(data = NA, nrow = n_years, ncol = n_fish_fleets, byrow = TRUE)
  
    for(f in 1:n_fish_fleets) {
      
      # Scale this up so that there is a linear ramp increase with bounds of minimum and maximum
      scaled_vector <- floor(
        ((fish_mort[,f,1] - min(fish_mort[,f,1], na.rm =T)) / 
           (max(fish_mort[,f,1], na.rm =T) - min(fish_mort[,f,1], na.rm =T))) * 
          (Input_Fish_N_Max[f]  - Input_N_Fish_Fixed[f]) + Input_N_Fish_Fixed[f]
      )
      
      # Input into matrix
      fish_Input_N_mat[,f] <- scaled_vector
      
    } # end fleet loop
  } # effective sample sizes vary as a function of fishery specific Fs
  
  # create empty matrix for survey Neff
  srv_Input_N_mat <- matrix(data = NA, nrow = n_years, ncol = n_srv_fleets, byrow = TRUE)
  
  # specify matrix for survey neff
  for(sf in 1:n_srv_fleets) {
    for(y in 1:(n_years-1)) {
      
      if(y < Surv_Start_yr[sf]) {
        srv_Input_N_mat[y,sf] <- 0 # force neff to be 0 when fishing is not occuring
      } else{
        srv_Input_N_mat[y,sf] <- Input_Srv_N_Max[sf]
      }
      
    } # end y loop
  } # end sf loop
  
  # output into environment
  Input_N_Fish <<- fish_Input_N_mat
  Input_N_Srv <<- srv_Input_N_mat
  
} # end function
