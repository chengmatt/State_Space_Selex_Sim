# Purpose: To feed the outputs of our OM into WHAM format!
# Creator: Matthew LH. Cheng
# Date: 10/31/22

make_input <- function(n_fleets = 1, n_indices = 2, Catch_CV_Val = 0) {
  
  require(tidyverse)
  
  # The structure of this is taken from https://timjmiller.github.io/wham/reference/prepare_wham_input.html (basic_info section)
  
  # Note: A lot of these objects are getting called in from the environment
  # when we first specify the OM!
  
# Set up ------------------------------------------------------------------
  
  # Create a list to store objects 
  input <- list() 
  
  # Vector of ages - the last being the plus group
  input$ages <- as.integer(ages) 
  
  # Years (Fish start - end year because our OM lags by a year)
  input$years <- Fish_Start_yr:(n_years-1) 
  
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
  

# Catch related quantities (Fishery) ------------------------------------------------
  
  # Catch CV matrix(length(years), n fleets)
  input$catch_cv <- matrix(Catch_CV_Val, n_y, input$n_fleets) # Catch CV
  
  # Catch effective sample size for catch_paa (comps for fishery) matrix(length(years), n_fleets)
  input$catch_Neff <- matrix(N_eff_fish, n_y, input$n_fleets) # Effective sample size for catch
  
  # Specify whether to use fleet age comps - we will specify 0 as default since we will have
  # a fishery dependent index of abundance (length(years), n_fleets)
  input$use_catch_paa <- matrix(0, n_y, input$n_fleets)
  
  # Specify selectivity block pointers for the fishery fleet (length(years), n_fleets) 
  # Still not really sure what this does - presuambly its to make sure stuff is getting indexed properly? (i.e.,
  # is it an independent fleet, are there are time-blocks, etc?)
  input$selblock_pointer_fleets <- t(matrix(1:input$n_fleets, input$n_fleets, n_y)) 
  
  # Fishing mortality rates to initialize the model
  F_mat <- matrix(nrow = n_y, ncol = input$n_fleets)
  
  # Fill this in with random normal draws and stick into our list
  F_mat[,] <- abs(rnorm(n = length(input$years), mean = 0, sd = 0.1))
  input$`F` <- F_mat # Random starting values
  

  ### Catch -------------------------------------------------------------------

  # Aggregated catch - matrix(length(years) x  n_fleets) of aggregate catch (biomass)
  agg_catch <- matrix(nrow = n_y, ncol = input$n_fleets)
  
  # Need to munge to get into matrix format
    catch_df <- melt(Catch_at_age) %>%
    drop_na() %>% # drop nas in the last year
    rename(Year = Var1, Age = Var2, # Rename varialbes
           Sim = Var3, Catch = value) %>%
    mutate(Year = parse_number(as.character(Year)),
           Sim = parse_number(as.character(Sim))) %>% # Parse number for year and simulation
    group_by(Year, Sim) %>% # summarize and aggregate catch
    summarize(Catch = sum(Catch, na.rm = TRUE)) %>%
    filter(Year >= Fish_Start_yr) # Filter out to the start of the fishing year
  
  # Now stick that into the matrix and put it into the input list
  agg_catch[, 1] <- catch_df$Catch
  input$agg_catch <- agg_catch
  

  ### Proportions at age ------------------------------------------------------

  # Create array to store values in catch_paa[fleet, year, age]
  catch_paa <- array(dim = c(input$n_fleets, n_y, n_a))
  
  # Munge the comp dataframe into the correct format
  fish_comps <- melt(Fish_Age_Comps) %>%
    drop_na() %>% # drop na values
    rename(Year = Var1, Age = Var2,
           Sim = Var3, Count = value) %>% # Parse number for year, age, and simulation
    mutate(Year = parse_number(as.character(Year)),
           Sim = parse_number(as.character(Sim)),
           Age = parse_number(as.character(Age))) %>% 
    group_by(Year) %>%
    mutate(Sum = sum(Count),
           Count = Count/Sum) # get proportions here
  
  # # Create a vector of years and ages we can loop through
  year_vec <- unique(fish_comps$Year)
  age_vec <- unique(fish_comps$Age)

  # Now, feed it into the catch_paa array
  for(y in 1:ncol(catch_paa)) {

    for(a in 1:n_a) {

      # Filter the dataset to a given year and age and stick it into the dataframe
      catch_paa[, y,a] <- fish_comps[fish_comps$Year == year_vec[y] &  fish_comps$Age == age_vec[a], ]$Count

      # print(a)
    } # end age loop

  } # end year loop

  # Put this into the list
  input$catch_paa <- catch_paa

  

# Indices and related quantities -----------------------------------------------------------------

  # Note: Indices are biomass based
  
  # CV for indices of abundance - rows = years, indices = cols (first col = fishery, second = survey)
  index_cv <- matrix(nrow = n_y, ncol = input$n_indices)
  index_cv[,1] <- Fishery_CV # Fishery CV
  # Specify survey cv - make sure to leave the NAs in 
  index_cv[length(Fish_Start_yr:Surv_Start_yr):n_y,2] <- Survey_CV
  # Put this into our list
  input$index_cv <- index_cv
  
  
  # Specify effective sample sizes for the fishery and survey (first col = fishery, second = survey)
  index_Neff <- matrix(nrow = n_y, ncol = input$n_indices)
  index_Neff[,1] <- N_eff_fish
  # Specify effective N or input for survey, but leave the NAs in
  index_Neff[length(Fish_Start_yr:Surv_Start_yr):n_y,2] <- N_eff_surv
  # put this into our list
  input$index_Neff <- index_Neff
  
  
  # Set units for the indices - biomass based (should be 1)
  input$units_indices <- rep(1, input$n_indices) # Biomass
  # Set units for index proprotions at age - in abundance (should be 2)
  input$units_index_paa <- rep(1, input$n_indices) # Abundance
  
  
  # Specify whether to use these indices comps or not (length(years) x n_indices) 
  # row = years, col = indices (fishery followed by survey). Given that we turned comps off
  # in our catch_paa, we will turn this on to make sure these comps read into WHAM and are 
  # associated with the fishery fleet 
  use_index_paa <- matrix(nrow = n_y, ncol = input$n_indices) # Note (0s = don't use, 1 = use)
  # For fishery comps - use all
  use_index_paa[,1] <- 1
  
  # For survey comps - make sure to leave the 0s in there
  use_index_paa[1:(length(Fish_Start_yr:Surv_Start_yr)-1),2] <- 0
  
  # Next fill in the 1s for when the survey start
  use_index_paa[length(Fish_Start_yr:Surv_Start_yr):nrow(use_index_paa),2] <- 1
  
  # Put this into our list
  input$use_index_paa <- use_index_paa
  input$use_indices <- use_index_paa # also stick this into here as well
  
  # Specify selectivity block options so the fishery fleets get linked to the index fleet for 
  # the fishery index and comps. We need to make sure to specifty the selblock pointer fleets to
  # be the same as sleblock pointer indices for the fishery fleet and index of abundance (e.g., both
  # need to be 1; Tim Miller Personal Communication Nov 1 2022)
  input$selblock_pointer_indices <- t(matrix(1:input$n_indices, input$n_indices, n_y))
  
  # Specify fracion of year the indices are observed - time elapsed since observation
  input$fracyr_indices <- matrix(1, n_y, input$n_indices)
  
  
  ### Indices of abundance (Fishery + Survey) ---------------------------------------------------
  
  # Create aggregate indices matrix to store values in col 1 = fishery, col 2 = survey
  agg_indices <- matrix(nrow = n_y, ncol = input$n_indices)
  
  # Munge fishery first
  fish_idx_df <- melt(Fishery_Index) %>% 
    drop_na() %>% # drop nas in the last year
    rename(Year = Var1, Sim = Var2, Index = value) %>%  # Rename varialbes
    mutate(Year = parse_number(as.character(Year)),
           Sim = parse_number(as.character(Sim))) # get numbers from character strings
  
  # Munge survey index
  surv_idx_df <- melt(Survey_Index) %>%
    drop_na() %>% # drop nas in the last year
    rename(Year = Var1, Sim = Var2, Index = value) %>% # Rename varialbes
    mutate(Year = parse_number(as.character(Year)),
           Sim = parse_number(as.character(Sim))) # Parse number for year and simulation

  # Put fishery index into matrix
  agg_indices[,1] <- fish_idx_df$Index
  
  # Put survey index into matrix - for missing values, we will specify as NAs
  # length(Fish_Start_yr:Surv_Start_yr) - when the survey first starts relative to the fishery (fishery always
  # starts at year 1)
  agg_indices[length(Fish_Start_yr:Surv_Start_yr):n_y,2] <- surv_idx_df$Index
  
  # Put this into our list
  input$agg_indices <- agg_indices
  

  ### Proportions at age (Fishery and Survey) ---------------------------------
  
  # Given that we have an index of abundance, we will specify to use the age comps within the indices 
  # realm for WHAM purposes so it doesn't get counted twice in the obj fxn
  # Array dimensioned by [n_indices, year, and age] - first row should be for the fishery, followed by the survey
  index_paa <- array(dim = c(input$n_indices, n_y, n_a))
  
  # Respecify year and age vecs for the fishery
  year_vec <- unique(fish_comps$Year)
  age_vec <- unique(fish_comps$Age)
  
  # Use the fishery data munged from above
  for(y in 1:ncol(index_paa)) {
    
    for(a in 1:n_a) {
      
      # Filter the dataset to a given year and age and stick it into the dataframe
      index_paa[1, y, a] <- fish_comps[fish_comps$Year == year_vec[y] &  fish_comps$Age == age_vec[a], ]$Count
      
    } # end age loop
    
  } # end year loop
  
  # Munge survey data and repeat the above
  surv_comps <- melt(Survey_Age_Comps) %>%
    drop_na() %>% # drop na values
    rename(Year = Var1, Age = Var2,
           Sim = Var3, Count = value) %>% # Parse number for year, age, and simulation
    mutate(Year = parse_number(as.character(Year)),
           Sim = parse_number(as.character(Sim)),
           Age = parse_number(as.character(Age))) %>% 
    group_by(Year) %>%
    mutate(Sum = sum(Count),
           Count = Count/Sum) # get proportions here
  
  # Create a vector of years and ages we can loop through for the survey
  year_vec <- unique(surv_comps$Year)
  age_vec <- unique(surv_comps$Age)
  
  # Now, feed it into the index_paa array
  for(y in 1:length(year_vec)) {  # length year vec - because needs to correspond to # of unique survey years
    
    for(a in 1:length(input$ages)) {
      
      # Filter the dataset to a given year and age and stick it into the dataframe
      # The (length(Fish_Start_yr:Surv_Start_yr)-2) there is to make sure we are indexing stuff correctly 
      # We are just adding the y index iteratively to 24 (when survey first starts)
      index_paa[2, y + (length(Fish_Start_yr:Surv_Start_yr)-1) ,a] <- surv_comps[surv_comps$Year == year_vec[y] & 
                                                                                   surv_comps$Age == age_vec[a], ]$Count
      
    } # end age loop
    
  } # end year loop
  
  # Put this into our list
  input$index_paa <- index_paa

  
# Biologicals (Weight at age and Maturity at age) --------------------------------------
  
  # Specify fracyr_ssb - when to calcualte ssb
  input$fracyr_SSB <- rep(1, n_y)

  # Maturity at age - matrix needs to be row = years, cols = ages
  maturity <- matrix(nrow = n_y, ncol = n_a)
  # Fill this in with our matrix - filter matrix to dimensions of assessment
  maturity[,] <- mat_at_age[min(input$years):max(input$years),,]
  # Put into our list
  input$maturity <- maturity
  
  
  # Input weight at age (may need to tweak this potentially?)
  nwaa <- input$n_indices + input$n_fleets + 2
  # use for each index and fleet if applicable
  
  # Create array to store values
  waa <- array(dim = c(nwaa, n_y, n_a))
  
  # Loop through to fill these in
  for(i in 1:nwaa) {
    
    waa[i,,] <- wt_at_age[min(input$years):max(input$years),,]
    
  } # end nwaa loop
  
  # Now put this into our input list
  input$waa <- waa

  return(input)
  
} # end function
