# Purpose: To feed the outputs of our OM into WHAM format!
# Creator: Matthew LH. Cheng
# Date: 10/31/22

#' Title
#'
#' @param n_fleets number of fleets
#' @param n_indices number of indices
#' @param Catch_CV_Val CV for catch - needs to be > 0
#' @param catch_error Whether we want some error in catch - uses CV val - boolean
#' @param n_sims Which simulation we are at - for use in continuous loop
#' @param bias_process Wheter we want to bias correct process error
#' @param bias_obs Whether we want to bias correct observations
#' @param units_indices Units for abundance indices (1 = biomass, 2 = numbers) needs to be specified as a vector
#' @param units_index_paa Units for proportions at age (1 = biomass, 2 = numbers) needs to be specified as a vector
#' @param single_fleet Whether or not we want to prepare our data according to a single fleet, when
#' we actually have multiple fleets - Boolean
#' @param time_block Whether or not we want to prepare our data such that a time-block is possible.
#' @param block_period_sel When we want the time block to occur for fleet selex
#' @param block_period_idx When we want the time block to occur for idx selex
#' @param noFish_Idx Whether or not we want to use a fishery dependent index of abundance


make_wham_input <- function(n_fleets, n_indices, Catch_CV_Val = 0, catch_error = TRUE,
                       n_sims, bias_process = FALSE, bias_obs = FALSE, units_indices, units_index_paa,
                       single_fleet = FALSE, time_block = FALSE, block_period_sel = NULL,
                       block_period_idx = NULL, noFish_Idx = FALSE) {
  
  if(n_indices != length(units_indices) & n_indices != length(units_index_paa)) stop("units_indices and units_index_paa vectors need to be the same length as the number of indices")
  if(length(Catch_CV_Val) != n_fleets) stop("Catch CV values are not the same length as the number of fishery fleets specified")
  if(single_fleet == TRUE & n_fleets > 1) stop("Arguments n_fleets and single_fleet contradict each other. Make sure to have n_fleets as 1 if single_fleet = TRUE")
  require(tidyverse)
  
  # The structure of this is taken from https://timjmiller.github.io/wham/reference/prepare_wham_input.html (input_info section)
  
# Set up ------------------------------------------------------------------
  
  # Create a list to store objects 
  input <- list() 
  
  # Vector of ages - the last being the plus group
  input$ages <- as.integer(ages) 
  
  # Assessment years - starts at the first fishery start year
  input$years <- min(Fish_Start_yr):(n_years-1) 
  
  # Number of fishing fleets
  input$n_fleets <- n_fleets  
  
  # Specify number of indices we want to model in the assessment
  input$n_indices <- n_indices  

  # Create some shortcut variables
  n_y <- length(input$years) # Years
  n_a <- length(input$ages) # ages
  
  # Ages and years
  input$n_ages <- n_a
  input$n_years_model <- n_y
  
  # Bias correction
  input$bias_correct_process <- bias_process
  input$bias_correct_observation <- bias_obs
  
# Catch related quantities (Fishery) ------------------------------------------------
  
  # Catch CV matrix(length(years), n fleets)
  input$catch_cv <- matrix(Catch_CV_Val, n_y, input$n_fleets) # Catch CV
  
  if(single_fleet == TRUE) {
    # Prepare effective sample sizes as a single fleet. i.e., sum up multinomial ESS
    input$catch_Neff <- matrix(rowSums(Input_N_Fish[min(Fish_Start_yr):(min(Fish_Start_yr)+n_y-1),,drop = FALSE]))
  } else{
    # Catch effective sample size for catch_paa (comps for fishery) matrix(length(years), n_fleets)
    input$catch_Neff <- Input_N_Fish[min(Fish_Start_yr):(min(Fish_Start_yr)+n_y-1),,drop = FALSE ] # Effective sample size for catch
  } # if one single fleet = TRUE
  
  # Specify selectivity block pointers for the fishery fleet (length(years), n_fleets) 
  # Still not really sure what this does - presuambly its to make sure stuff is getting indexed properly? (i.e.,
  # is it an independent fleet, are there are time-blocks, etc?)
  if(time_block == TRUE) {

    # Create matrix to store point values in
    sel_block_mat <- matrix(ncol = input$n_indices, nrow = n_y)
    
    # First figure out how much we need to sequence to
    map_blk_sel <- sum(!is.na(unlist(block_period_sel))) + length(block_period_sel)
    
    # Create sequence here
    seq_blk_vals <- 1:map_blk_sel
    
    for(i in 1:ncol(sel_block_mat)) {
      
      if(is.na(block_period_sel[i])) { # if this is not a time block
        blk_val <- seq_blk_vals[1] # Get the first element for a block value
        sel_block_mat[,i] <- blk_val # Input block value into the column here
        # Remove the inputted unique values
        seq_blk_vals <- seq_blk_vals[(!seq_blk_vals %in% c(unique(sel_block_mat[,i])))] 
      } else{
        
        # Get block sequences and add last year into sequence
        seq_blks <- c(block_period_sel[[i]], n_y)
        # Get block values from the generated sequence earlier
        blk_val <- seq_blk_vals[1:(length(seq_blks))]
        
        counter <- 1 # Create counter
        for(j in 1:length(seq_blks)) {
          
          # Index the block value here
          idx_blk_val <- blk_val[j]
          
          # Put the indexed block value into said location
          sel_block_mat[c(counter:seq_blks[j]),i] <- idx_blk_val
          
          counter <- 1 + seq_blks[j] # Add counter to update sequence
          
        } # end j loop
        
        # Remove these unique values from our unique block values
        seq_blk_vals <- seq_blk_vals[(!seq_blk_vals %in% c(unique(sel_block_mat[,i])))] 
      }
      
    } # end i
  
    input$selblock_pointer_fleets <- sel_block_mat
    
  } else{
    input$selblock_pointer_fleets <- t(matrix(1:input$n_fleets, input$n_fleets, n_y)) 
  } # not a time block

  ### Catch -------------------------------------------------------------------

  # Aggregated catch - matrix(length(years) x  n_fleets) of aggregate catch (biomass)
  agg_catch <- matrix(nrow = n_y, ncol = input$n_fleets)
  
  # Need to munge to get into matrix format
    catch_df <- melt(Catch_at_age) %>%
    drop_na() %>% # drop nas in the last year
    rename(Year = Var1, Age = Var2, # Rename varialbes
           Fleet = Var3, Sex = Var4, Sim = Var5, Catch = value) %>%
    mutate(Year = parse_number(as.character(Year)),
           Sim = parse_number(as.character(Sim))) %>%  # Parse number for year and simulation
      filter(Year >= Fish_Start_yr[1],
             Sim == n_sims)
      
    if(single_fleet == TRUE) {
      # Prepare catch by summing catch across fleets
      catch_mat <- matrix(with(catch_df, tapply(Catch, list(Year), FUN = sum)))
    } else{
      # Get catch and sum by years and fleets - turn into matrix 
      catch_mat <- matrix(with(catch_df, tapply(Catch, list(Year, Fleet), FUN = sum))[,1:input$n_fleets], nrow = n_y)
    }
    
  # Now stick that into the matrix and put it into the input list - add some lognormal
  # noise to this
    
    if(catch_error == TRUE) {
      for(f in 1:input$n_fleets) {
        agg_catch[, f] <- catch_mat[,f] * exp(rnorm(length(catch_mat[,f]), 0, log((Catch_CV_Val[f]^2) +1)))
      } # end f loop
    } else{
      for(f in 1:input$n_fleets) agg_catch[, f] <- catch_mat[,f]
    } # else = no error
    
    # make 0 catches into NAs
    agg_catch[agg_catch == 0] <- 0
    
  input$agg_catch <- agg_catch
  

  ### Proportions at age ------------------------------------------------------

  # Create array to store values in catch_paa[fleet, year, age]
  catch_paa <- array(dim = c(input$n_fleets, n_y, n_a))
  
  # Munge the comp dataframe into the correct format. 
  fish_comps <- melt(Fish_Age_Comps) %>%
    drop_na() %>% # drop na values
    rename(Year = Var1, Age = Var2,
           Fleet = Var3, Sex = Var4, Sim = Var5,
           Count = value) %>% # Parse number for year, age, and simulation
    mutate(Year = parse_number(as.character(Year)),
           Sim = parse_number(as.character(Sim)),
           Age = parse_number(as.character(Age))) %>% 
    filter(Sim == n_sims) # filter to the simulation number
  
  if(single_fleet == TRUE) {
    # Group by and summarize to one single fleet
    fish_comps <- fish_comps %>% 
      group_by(Year, Age, Sex, Sim) %>% 
      summarize(Count = sum(Count)) %>% 
      mutate(Fleet = "Fish_Fleet_1") %>% 
      filter(Sim == n_sims)
  } # if this is a single fleet, condense the above if it is multi fleet into a single fleet
  
  for(f in 1:input$n_fleets) {
    
    # Filter to fishery fleet f
    fish_comps_fleet <- fish_comps[fish_comps$Fleet == paste("Fish_Fleet",f, sep = "_"),]
    
    # Get fish comps and proportions
    fishcomps_df <- with(fish_comps_fleet, tapply(Count,list(Year), FUN=function(x) { x/sum(x) }))
    
    # Get length of unique years to loop through
    unique_yrs <- length(unique(as.numeric(names(fishcomps_df))))
    
    # Get number of years for use in y loop as a counter, so we can idnex properply for uneven years
    yr_counter <- input$n_years_model - length(unique(as.numeric(names(fishcomps_df))))
    
    # Now, feed it into the catch_paa array
    for(y in 1:unique_yrs) {
      
      # Filter the dataset to a given year and age and stick it into the dataframe
      catch_paa[f, y + yr_counter,] <- fishcomps_df[[y]]
      
    } # end year loop
    
  } # end f for fishing fleets

  # Put this into the list
  input$catch_paa <- catch_paa

  

# Indices and related quantities -----------------------------------------------------------------

  # Note: Indices are numbers based
  
  # CV for indices of abundance - rows = years, indices = cols (first col = fishery, second = survey)
  index_cv <- matrix(nrow = n_y, ncol = input$n_indices)
  
  # CVs need to be specified in the order of the fishery fleets, followed by the survey fleets
  for(f in 1:input$n_fleets) {
    index_cv[,f] <- fish_CV[f]
  } # end f loop for fishery
  
  # Do the same for the survey fleets
  for(sf in 1:n_srv_fleets) {
    index_cv[,(input$n_fleets + sf)] <- srv_CV[sf]
  } # end sf loop
  
  # Put this into our list
  input$index_cv <- index_cv

  # Specify effective sample sizes for the fishery and survey (first cols = fishery, followed by the survey)
  index_Neff <- matrix(nrow = n_y, ncol = input$n_indices)
  
  if(single_fleet == TRUE) {
    # Prepare effective sample sizes as a single fleet. i.e., sum up multinomial ESS
    Neff_fish_mat <- matrix(rowSums(Input_N_Fish[min(Fish_Start_yr):(min(Fish_Start_yr)+n_y-1),,drop = FALSE]))
  } else{
    # Filter to the specified years in the model
    Neff_fish_mat <- Input_N_Fish[min(Fish_Start_yr):(min(Fish_Start_yr)+n_y-1),,drop = FALSE]
  } # if one single fleet = TRUE

  for(f in 1:input$n_fleets) {
    index_Neff[,f] <- Neff_fish_mat[,f]
  } # end f fishery loop
  
  # do the same above, but for the survey - keeping the min as fish start year because
  # we want to retain the 0s
  Neff_srv_mat <- Input_N_Srv[min(Fish_Start_yr):(min(Fish_Start_yr)+n_y-1),,drop = FALSE ]
  
  for(sf in 1:n_srv_fleets) {
    index_Neff[,(input$n_fleets + sf)] <- Neff_srv_mat[,sf]
  } # end f fishery loop
  
  # put this into our list
  input$index_Neff <- index_Neff
  
  # Set units for the indices 1 = biomass, 2 = abundance
  input$units_indices <- units_indices
  # Set units for index proprotions at age - abundnace e.g., numbers for both age comps 1 = biomass, 2 = abundance
  input$units_index_paa <- units_index_paa
  
  # Specify whether to use these indices comps or not (length(years) x n_indices) 
  # row = years, col = indices (fishery followed by survey). Given that we turned comps off
  # in our catch_paa, we will turn this on to make sure these comps read into WHAM and are 
  # associated with the fishery fleet 
  use_index_paa <- matrix(nrow = n_y, ncol = input$n_indices) # Note (0s = don't use, 1 = use)
  
  # Create fishery and survey start years vector relative to the model start year for use in loops below
  fish_use_yrs <- Fish_Start_yr - min(Fish_Start_yr)
  srv_use_yrs <- Surv_Start_yr - min(Fish_Start_yr)
  
  for(f in 1:input$n_fleets) { # loop through
    # For fishery comps - use all when avaliable
    for(y in 1:n_y){
      if(y <= fish_use_yrs[f]) {
        use_index_paa[y,f] <- 0
      } else{
        use_index_paa[y,f] <- 1
      } # if else for when to use indices
    } # end y loop
  } # end f fishery fleet loop
  
  # For survey fleets, do the same, but leave in 0s during periods where no sruvey is occuring
  for(sf in 1:n_srv_fleets) {
    for(y in 1:n_y) {
      if(y <= srv_use_yrs[sf]) {
        use_index_paa[y,(input$n_fleets + sf)] <- 0
    } else{
      use_index_paa[y,(input$n_fleets + sf)] <- 1
    } # end if else statement for specify when to use survey 
  } # end year loop
} # end sf loop
  
  # Put these into our list
  input$use_index_paa <- use_index_paa
  input$use_indices <- use_index_paa 
  
  # Specify selectivity block options so the fishery fleets get linked to the index fleet for 
  # the fishery index and comps. We need to make sure to specifty the selblock pointer fleets to
  # be the same as selblock pointer indices for the fishery fleet and index of abundance (e.g., both
  # need to be 1; Tim Miller Personal Communication Nov 1 2022)
  if(time_block == TRUE) {

    # Create matrix to store point values in
    sel_block_idx_mat <- matrix(ncol = input$n_indices, nrow = n_y)
    
    # First figure out how much we need to sequence to
    map_blk_idx <- sum(!is.na(unlist(block_period_idx))) + length(block_period_idx)
    
    # Create sequence here
    seq_blk_vals <- 1:map_blk_idx
    
    for(i in 1:ncol(sel_block_idx_mat)) {
      
      if(is.na(block_period_idx[i])) { # if this is not a time block
        blk_val <- seq_blk_vals[1] # Get the first element for a block value
        sel_block_idx_mat[,i] <- blk_val # Input block value into the column here
        # Remove the inputted unique values
        seq_blk_vals <- seq_blk_vals[(!seq_blk_vals %in% c(unique(sel_block_idx_mat[,i])))] 
      } else{
        
        # Get block sequences and add last year into sequence
        seq_blks <- c(block_period_idx[[i]], n_y)
        # Get block values from the generated sequence earlier
        blk_val <- seq_blk_vals[1:(length(seq_blks))]

        counter <- 1
        for(j in 1:length(seq_blks)) {
          
          # Index the block value here
          idx_blk_val <- blk_val[j]
          
          # Put the indexed block value into said location
          sel_block_idx_mat[c(counter:seq_blks[j]),i] <- idx_blk_val
          
          counter <- 1 + seq_blks[j] # Add counter to update sequence
          
        } # end j loop
        
        # Remove these unique values from our unique block values
        seq_blk_vals <- seq_blk_vals[(!seq_blk_vals %in% c(unique(sel_block_idx_mat[,i])))] 
        }
      
    } # end i
    
    input$selblock_pointer_indices <- sel_block_idx_mat
    
  } else{
    input$selblock_pointer_indices <- t(matrix(1:input$n_indices, input$n_indices, n_y))
  } # not a time block
  
  # Specify fracion of year the indices are observed - time elapsed since observation
  input$fracyr_indices <- matrix(0, n_y, input$n_indices)
  
  
  ### Indices of abundance (Fishery + Survey) ---------------------------------------------------
  
  # Create aggregate indices matrix to store values in col 1 = fishery, col 2 = survey
  agg_indices <- matrix(nrow = n_y, ncol = input$n_indices)
  
  # Munge fishery first
  fish_idx_df <- melt(Fishery_Index_Agg) %>% 
    drop_na() %>% # drop nas in the last year
    rename(Year = Var1, Fleet = Var2, Sim = Var3, Index = value) %>%  # Rename varialbes
    mutate(Year = parse_number(as.character(Year)),
           Sim = parse_number(as.character(Sim))) %>%  # get numbers from character strings
    filter(Sim == n_sims)
  
  if(single_fleet == TRUE) {
    
    # Group by and summarize across fleets by taking the average of the two
    fish_idx_df <- fish_idx_df %>% 
      group_by(Year, Sim) %>% 
      summarize(Index = mean(Index)) %>% 
      mutate(Fleet = "Fish_Fleet_1") %>% 
      filter(Sim == n_sims)
    
  } # if we are using a single fleet structure
  
  # Munge survey index
  surv_idx_df <- melt(Survey_Index_Agg) %>%
    drop_na() %>% # drop nas in the last year
    rename(Year = Var1, Fleet = Var2, Sim = Var3, Index = value) %>%  # Rename varialbes
    mutate(Year = parse_number(as.character(Year)),
           Sim = parse_number(as.character(Sim))) %>%  # Parse number for year and simulation
    filter(Sim == n_sims)
  
  # Loop through to put the fishery abundance indices into our matrix
  for(f in 1:input$n_fleets) {
    
    # Filter to fishery fleet f and only fishery years - note that although we are inputting data
    # from fisheries during the same years, those data will not get used because we specified so in
    # use_index and use_index_paa
    fish_idx_fleet <- fish_idx_df[fish_idx_df$Fleet == paste("Fish_Fleet",f, sep = "_") &
                                  fish_idx_df$Year %in% c(min(Fish_Start_yr):(n_years-1)),]
    
    # Get length of unique years to loop through
    unique_yrs <- length(unique(fish_idx_fleet$Year))
    
    # Get number of years for use in y loop as a counter, so we can idnex properply for uneven years
    yr_counter <- input$n_years_model - unique_yrs
    
    for(y in 1:unique_yrs) {
      # Put into our matrix
      agg_indices[(y + yr_counter),f] <- fish_idx_fleet$Index[y]
    }
    
  } # end f fishey loop

  # Do the above but for survey indices
  for(sf in 1:n_srv_fleets) {
    
    # Filter to survey fleet sf 
    srv_idx_fleet <- surv_idx_df[surv_idx_df$Fleet == paste("Srv_Fleet",sf, sep = "_") &
                                   surv_idx_df$Year %in% c(min(Fish_Start_yr):(n_years-1)),]
    
    # Get length of unique years to loop through
    unique_yrs <- length(unique(srv_idx_fleet$Year))
    
    # Get number of years for use in y loop as a counter, so we can idnex properply for uneven years
    yr_counter <- input$n_years_model - unique_yrs
    
    for(y in 1:unique_yrs) {
      # Put into our matrix
      agg_indices[(y + yr_counter),(input$n_fleets + sf)] <- srv_idx_fleet$Index[y]
    } # end year loop
    
  } # end sf loop
  
  # Put this into our list
  input$agg_indices <- agg_indices
  

  ### Proportions at age (Fishery and Survey) ---------------------------------
  
  # Given that we have an index of abundance, we will specify to use the age comps within the indices 
  # realm for WHAM purposes so it doesn't get counted twice in the obj fxn
  # Array dimensioned by [n_indices, year, and age] - first row should be for the fishery, followed by the survey
  index_paa <- array(dim = c(input$n_indices, n_y, n_a))
  
  for(f in 1:input$n_fleets) {
    
    # Filter to fishery fleet f
    fish_comps_fleet <- fish_comps[fish_comps$Fleet == paste("Fish_Fleet",f, sep = "_"),]
    
    # Get fish comps and proportions
    fishcomps_df <- with(fish_comps_fleet, tapply(Count,list(Year), FUN=function(x) { x/sum(x) }))
    
    # Get length of unique years to loop through
    unique_yrs <- length(unique(as.numeric(names(fishcomps_df))))
    
    # Get number of years for use in y loop as a counter, so we can idnex properply for uneven years
    yr_counter <- input$n_years_model - unique_yrs
    
    # Now, feed it into the catch_paa array
    for(y in 1:unique_yrs) {
      
      # Filter the dataset to a given year and age and stick it into the dataframe
      index_paa[f, (y + yr_counter),] <- fishcomps_df[[y]]
      
    } # end year loop
    
  } # end f for fishing fleets
  
  # Munge survey data and repeat the above
  surv_comps <- melt(Survey_Age_Comps) %>%
    drop_na() %>% # drop na values
    rename(Year = Var1, Age = Var2, Fleet = Var3,
           Sex = Var4, Sim = Var5, Count = value) %>% # Parse number for year, age, and simulation
    mutate(Year = parse_number(as.character(Year)),
           Sim = parse_number(as.character(Sim)),
           Age = parse_number(as.character(Age))) %>% 
    filter(Sim == n_sims)
  
  for(sf in 1:n_srv_fleets) {
    
    # Filter to fishery fleet f
    srv_comps_fleet <- surv_comps[surv_comps$Fleet == paste("Srv_Fleet",sf, sep = "_"),]
    
    # Get fish comps and proportions
    srv_comps_df <- with(srv_comps_fleet, tapply(Count,list(Year), FUN=function(x) { x/sum(x) }))
    
    # Get length of unique years to loop through
    unique_yrs <- length(unique(as.numeric(names(srv_comps_df))))
    
    # Get number of years for use in y loop as a counter, so we can idnex properply for uneven years
    yr_counter <- input$n_years_model - unique_yrs
    
    # Now, feed it into the catch_paa array
    for(y in 1:unique_yrs) {
      
        # Filter the dataset to a given year and age and stick it into the dataframe
        index_paa[(input$n_fleets + sf), (y + yr_counter),] <- srv_comps_df[[y]]

    } # end year loop
  } # end f for fishing fleets
  
  # Put this into our list
  input$index_paa <- index_paa

  
# Biologicals (Weight at age and Maturity at age) --------------------------------------
  
  # Specify fracyr_ssb - when to calcualte ssb
  input$fracyr_SSB <- rep(0, n_y)

  # Maturity at age - matrix needs to be row = years, cols = ages
  maturity <- matrix(nrow = n_y, ncol = n_a)
  # Fill this in with our matrix - filter matrix to dimensions of assessment
  maturity[,] <- mat_at_age[min(input$years):max(input$years),,1,n_sims] # filter to only sex 1 because wham is not sex specific
  # Put into our list
  input$maturity <- maturity
  
  # Input weight at age (may need to tweak this potentially?)
  nwaa <- input$n_indices + input$n_fleets + 2
  # use for each index and fleet if applicable
  
  # Create array to store values
  waa <- array(dim = c(nwaa, n_y, n_a))
  
  # Loop through to fill these in
  for(i in 1:nwaa) {
    
    # Filter to sex 1 because wham is not sex specific
    waa[i,,] <- wt_at_age[min(input$years):max(input$years),,1,n_sims]
    
  } # end nwaa loop
  
  # Now put this into our input list
  input$waa <- waa
  
  if(noFish_Idx == TRUE) { # removing the first column for a fishery dependent index
    input$n_indices <- 1
    # input$use_catch_paa[] <- 1
    input$index_cv <-  matrix(input$index_cv[,-1])
    input$index_Neff <-  matrix(input$index_Neff[,-1])
    input$units_indices <- input$units_indices[-1]
    input$units_index_paa <- input$units_index_paa[-1]
    input$use_index_paa <- matrix(input$use_index_paa[,-1])
    input$use_indices <- matrix(input$use_indices[,-1])
    input$selblock_pointer_indices <- matrix( input$selblock_pointer_indices[,-1])
    input$fracyr_indices <- matrix(input$fracyr_indices[,-1])
    input$agg_indices <- matrix(input$agg_indices[,-1])
    input$index_paa <- input$index_paa[-1,,]
  }
  
  return(input)
  
} # end function
