# Purpose: To specify F pattern scenarios for the fishery
# Creator: Matthew LH. Cheng
# Date: 10/31/22


# Specify Fs --------------------------------------------------------------

#' @param Start_F Starting F vector
#' @param Fish_Start_yr Starting Fishing years for fleets
#' @param n_years Number of years in simulation
#' @param max_rel_F_M Maximum F value relative to the M
#' @param desc_rel_F_M Descending F value relative to the M
#' @param F_type Specify type of fishing mortality pattern
#' @param mean_nat_mort Mean natural mortality
#' @param yr_chng If we are simulating a fleet structure change, when do we want this change to occur at 
#' Note that this is only applicable for (Const_Inc, Contrast, Const_Ramp_Const, and Contrast_Const)
#' @param yr_chng_end If we are simulating a fleet structure change, when do we want this fleet structure change
#' to end. Applicable for (Const_Ramp_Const, Contrast_Const)

specify_F_pattern <- function(Start_F, Fish_Start_yr, F_type, n_years, max_rel_F_M, desc_rel_F_M,
                              mean_nat_mort, yr_chng = NULL, yr_chng_end = NULL) {
  
  if(F_type == "Contrast") { # Contrast F type
    
    if(desc_rel_F_M > max_rel_F_M) stop("Descending limb of fishing mortality contrast is larger than ascending limb! Specify another pattern or adjust the relative mortality rates to reflect a contrast pattern!")

    # Change point here
    chngpoint <- yr_chng
    
    # Increase fishing mortality (ramp up); here we are specifying fishing mortality relative to natural mortality.
    # We will allow it go over natural mortality first by 1.5x
    fish_mort_1 <- seq(Start_F, 
                       (max_rel_F_M * mean_nat_mort), 
                       length.out = length(Fish_Start_yr:chngpoint)) 
    
    # Now, decrease fihsing mortality (ramp down) - decrease to 0.75 of the mean natural mortality
    fish_mort_2 <- seq(((max_rel_F_M * mean_nat_mort)), 
                       (desc_rel_F_M * mean_nat_mort), 
                       length.out = length(chngpoint:(n_years-1)))
    
    # Now bind all of these fish morts into a vector
    F_vec <- c(fish_mort_1, fish_mort_2[-1]) # Remove first element so it connects smoothly
    
  } # if statement for F pattern that is contrasting
  
  if(F_type == "Constant") { # Constant F type
    
    # Constant fishing mortality - 0.75 of natural mortality
    F_vec <- rep(Start_F, times = length(Fish_Start_yr:(n_years-1))) 
    
  } # if statement for F pattern that is Constant
  
  if(F_type == "Increase") { # Increasing F type
    
    F_vec <- seq(Start_F, (max_rel_F_M * mean_nat_mort), length.out = length(Fish_Start_yr:(n_years-1)))
    
  } # if statement for F pattern that is increasing
  
  if(F_type == "Decrease") { # Decreasing F type
    
    F_vec <- seq(Start_F, (max_rel_F_M * mean_nat_mort), length.out = length(Fish_Start_yr:(n_years-1)))
    
  } # if statement for F pattern that is decreasing
  

  if(F_type == "Increase_Plat") { # if increases and plateaus out
    
    # Calculate plateau year
    plat_year <- round(0.95 * (n_years-1))  # make it 0.5 of the number of years in the simulation
    
    # Calculate increasing f mort
    increasing <- seq(Start_F, (max_rel_F_M * mean_nat_mort), length.out = length(Fish_Start_yr:(plat_year-1))) 
    
    # Get vector of plateau fmorts
    plat <- rep((max_rel_F_M * mean_nat_mort), length.out = length((plat_year):(n_years-1)))
    
    F_vec <- c(increasing, plat)
    
  } # if statement for F pattern that is increasing
  
  if(F_type == "Const_Inc_or_Dec") { # Constant and Increasing or Decreasing F type
    
    # Change point for fmort
    chngpoint <- yr_chng
    # F stays constant initially
    constant_F <- seq(Start_F, Start_F, length.out = length(Fish_Start_yr:chngpoint))
    # Change point for F
    change_F <- seq(Start_F, (max_rel_F_M * mean_nat_mort), length.out = length(chngpoint:(n_years-1)))
    # F vector 
    F_vec <- c(constant_F, change_F[-1]) # Removing the first year from this sequence
    
  } # if fleet structure is constant and increases after
  
  if(F_type == "Const_Inc_or_Dec_Plat") { # Constant and Increasing or Decreasing F type followed by a plateau
    
    # Change point for fmort
    chngpoint <- yr_chng
    # F stays constant initially
    constant_F <- seq(Start_F, Start_F, length.out = length(Fish_Start_yr:chngpoint))
    # Change point for F
    change_F <- seq(Start_F, (max_rel_F_M * mean_nat_mort), length.out = length(chngpoint:(yr_chng_end-1)))
    # F vector 
    F_vec <- c(constant_F, 
               change_F[-1],  # Removing the first year from this sequence in change F
               rep(change_F[length(chngpoint:yr_chng_end)-1], length(yr_chng_end:(n_years-1))))
    
  } # if fleet structure is constant and increases after
  
  if(F_type == "Const_Ramp_Const") { # Constant, Ramp (can be increase or decrease), 
    # and then remains constant
    
    chngpoint <- yr_chng     # Change point for fmort
    chngpoint_end <- yr_chng_end # when does change point end
    
    # F stays constant initially
    constant_F <- seq(Start_F, Start_F, length.out = length(Fish_Start_yr:chngpoint))
    # Change point for F (Ramp)
    change_F <- seq(Start_F, (max_rel_F_M * mean_nat_mort), length.out = length(chngpoint:chngpoint_end))
    # Change point end (plateaus/stays constant now)
    change_F_const <- seq((max_rel_F_M * mean_nat_mort), (max_rel_F_M * mean_nat_mort),
                          length.out = length(chngpoint_end:(n_years-1)))
    
    # Make F vector here
    F_vec <- c(constant_F, change_F[-1], change_F_const[-1])
  }
  
  if(F_type == "Contrast_Const") { # Constrast but declines to a constant rate
    
    chngpoint <- yr_chng    # Change point here
    chngpoint_end <- yr_chng_end # when does change point end

    # Increase fishing mortality (ramp up); here we are specifying fishing mortality relative to natural mortality.
    # We will allow it go over natural mortality first by 1.5x
    fish_mort_1 <- seq(Start_F, 
                       (max_rel_F_M * mean_nat_mort), 
                       length.out = length(Fish_Start_yr:chngpoint)) 
    
    # Now, decrease fihsing mortality (ramp down) - decrease to 0.75 of the mean natural mortality
    fish_mort_2 <- seq(((max_rel_F_M * mean_nat_mort)), 
                       (desc_rel_F_M * mean_nat_mort), 
                       length.out = length(chngpoint:chngpoint_end))
    
    # Finally, make this stay constant at a given rate
    fish_mort_3 <- seq(((desc_rel_F_M * mean_nat_mort)), 
                       (desc_rel_F_M * mean_nat_mort), 
                       length.out = length(chngpoint_end:(n_years-1)))
    
    # Now bind all of these fish morts into a vector
    F_vec <- c(fish_mort_1, fish_mort_2[-1], fish_mort_3[-3]) # Remove first element so it connects smoothly
    
  }
  
  return(F_vec)
  
} # end function


# Get fishing mortality rates ---------------------------------------------

#' @param Start_F Starting F vector
#' @param Fish_Start_yr Starting Fishing years for fleets
#' @param n_years Number of years in simulation
#' @param max_rel_F_M Maximum F value relative to the M
#' @param desc_rel_F_M Descending F value relative to the M
#' @param F_type Specify type of fishing mortality pattern
#' @param mean_nat_mort Mean natural mortality

get_Fs <- function(Start_F, Fish_Start_yr, F_type, n_years, max_rel_F_M, desc_rel_F_M,
                   mean_nat_mort, yr_chng = NULL, yr_chng_end = NULL) {
  
  if(n_fish_fleets != length(Start_F) & n_fish_fleets != length(Fish_Start_yr) &
     n_fish_fleets != length(F_type) & n_fish_fleets != length(max_rel_F_M) &
     n_fish_fleets != length(desc_rel_F_M)
     ) stop("Vector length of Start_F, Fish_Start_yr, F_type, max_rel_F_M, desc_rel_F_M does not equal the number of fleets specified.")
  
  # Loop through fleets to specify Fs
  for(f in 1:n_fish_fleets) {
    
    # Set Fs to 0 before fishery starts
    fish_mort[1:(Fish_Start_yr-1)[f],f,] <- 0
    
    # Set F pattern to fishery - feed in inputs to arguments above
    fish_mort[Fish_Start_yr[f]:(n_years-1),f,] <- specify_F_pattern(Start_F = Start_F[f], # Starting F
                                              Fish_Start_yr = Fish_Start_yr[f], # When fishery starts
                                              F_type = F_type[f], # What type of F pattern
                                              n_years = n_years, # Number of model years
                                              max_rel_F_M = max_rel_F_M[f], # Apical F relative to M
                                              desc_rel_F_M = desc_rel_F_M[f], # Descending F after apical relative to M
                                              mean_nat_mort = mean_nat_mort, # Mean Natural Mortality
                                              yr_chng = yr_chng, # what year fleet structure begins to change
                                              yr_chng_end = yr_chng_end) # year fleet structure chagne ends
  } # end f loop
  
  # Output into environment
  fish_mort <<- fish_mort
  Fish_Start_yr <<- Fish_Start_yr 
  
} # end function



