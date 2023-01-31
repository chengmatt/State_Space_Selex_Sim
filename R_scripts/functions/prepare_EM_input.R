# Creator: Matthew LH. Cheng
# Purpose: To preapre inputs for the estimation model in TMB
# Date 1/10/23

#' Title Prepares data inputs for the EM 
#'
#' @param years vector of years
#' @param catch_cv vector of catch CVs
#' @param F_Slx_Blocks_Input matrix of when you want to block
#' @param S_Slx_Blocks_Input same as above
#' @param use_catch Boolean (TRUE OR FALSE) TRUE = Use, FALSE = don't use
#' @param use_fish_index same as above
#' @param use_srv_index same as above
#' @param use_fish_comps Boolean (TRUE OR FALSE) TRUE = USE, FALSE = don't use
#' @param use_srv_comps same above abovce
#' @param rec_model recruitment model == 0 (mean recruitment)
#' @param F_Slx_Model_Input Fishery selectivity model character (logistic, gamma, double_logistic)
#' @param S_Slx_Model_Input same as above
#' @param sim simulation indexing
#'
#' @return
#' @export
#'
#' @examples
#' Note that n_years here refers to the all the rows of the array
prepare_EM_input <- function(years,
                          n_fleets, 
                          catch_cv,
                          F_Slx_Blocks_Input,
                          S_Slx_Blocks_Input,
                          use_catch,
                          use_fish_index,
                          use_srv_index, 
                          use_fish_comps,
                          use_srv_comps,
                          F_Slx_Model_Input,
                          S_Slx_Model_Input,
                          rec_model = 0,
                          sim) {
  
  # Make data into a list 
  data <- list()
  # Make parameters into a list
  pars <- list()
  
  # Specify dimensions 
  ages <- ages # ages
  Fish_Start_yr <- Fish_Start_yr # Fishery start year
  n_sexes <- n_sex # sexes
  n_years <- nrow(Catch_agg) # number of years
  n_fleets <- n_fleets # number of fleets
  n_fish_comps <- n_fleets # number of fleets = number of fish comps
  n_fish_indices <- n_fleets # number of fish indices = number of fleets
  n_srv_comps <- dim(Surv_selex_at_age)[3] # 3rd dimension of this array = number of survey fleets
  n_srv_indices <- dim(Survey_Index_Agg)[2] # 2nd dimension of this array = number of survey fleets
  
# Catch -------------------------------------------------------------------

  if(n_fleets == 1) { # if single fleet - sum across or leave as is
    obs_catches <- as.matrix(apply( as.matrix(Catch_agg[(Fish_Start_yr[1]:(n_years-1)),,sim]), 1, FUN = sum),
                             ncol = 1)
    
    # Get catch weighting if there is more than one fleet
    if(dim(Catch_agg)[2] > 1) {
      # sum across rows
      catch_sum <- rowSums(as.matrix(Catch_agg[(Fish_Start_yr[1]:(n_years-1)),,sim]))
      catch_weight <- as.matrix(Catch_agg[(Fish_Start_yr[1]:(n_years-1)),,sim]) / catch_sum
    } # if catch weighting

  } else{ # multi fleet - leave as is
    obs_catches <- matrix(Catch_agg[(Fish_Start_yr[1]:(n_years-1)),,sim], nrow = length(years), ncol = n_fleets)
  } # end else
  

# Fishery Age Comps -------------------------------------------------------

  obs_fish_age_comps <- array(data = 0, dim = c(length(years), length(ages), n_fleets, n_sexes))
  
  if(n_fleets == 1) { # one fleet = sum across
    
    # Loop through to sum - index the 3rd dimension to get fleets
    for(s in 1:n_sexes) {
      for(f in 1:dim(Fish_Age_Comps)[3]) {
        # Filter to save as an object
        fish_age_comps <-  Fish_Age_Comps[Fish_Start_yr[1]:(n_years - 1),,f,s,sim] * catch_weight[,f]
        # Increment comps - fixing fleet index to 1 here
        obs_fish_age_comps[,,1,s] <- obs_fish_age_comps[,,1,s] + fish_age_comps
      } # end f loop
    } # end s loop
    
    # Empty array to store values in
    obs_fish_age_Neff <- array(data = NA, dim = c(length(years), n_fleets, n_sexes))
    
      # Effective sample size
    for(s in 1:n_sexes) {
      obs_fish_age_Neff[,1,s] <- rowSums(floor(obs_fish_age_comps[,,,s]))
    } # s loop
    
    # Now, apply the proportion function over a single fleet
    for(s in 1:n_sexes) {
      obs_fish_age_comps[,,,s] <- t(apply(obs_fish_age_comps[,,,s], MARGIN = 1, 
                                          FUN=function(x) { x/sum(x) }))
    } # s loop
    
  } else{ # more than one fleet here
    
    for(f in 1:n_fleets) { # needs to loop through transpose and apply for easy retnetion of array dimensions
      for(s in 1:n_sexes) {
        obs_fish_age_comps[,,f,s] <- t(
          apply(Fish_Age_Comps[Fish_Start_yr[1]:(n_years - 1),,f,s,sim], MARGIN = 1, 
                FUN=function(x) { x/sum(x) })
        )
      } # end s loop
    } # end f loop
    
    # Effective Sample Sizes
    obs_fish_age_Neff <- matrix(fish_Neff[Fish_Start_yr[1]:(n_years - 1),], 
                                    nrow = length(years), ncol = n_fleets)
    
  } # end else
  

# Survey Age Comps --------------------------------------------------------
# Right now, this is only able to accomdate one single survey fleet
  
  obs_srv_age_comps <- array(data = NA, dim = c(length(years), length(ages), n_srv_comps, n_sexes))

  # Apply function to get proportions and munge into matrix
  # Now, apply the proportion function over a single fleet
  for(s in 1:n_sexes) {
    obs_srv_age_comps[,,,s] <- t(apply(Survey_Age_Comps[Fish_Start_yr[1]:(n_years - 1),,,s,sim],
                                        MARGIN = 1, FUN=function(x) { x/sum(x) }))
  } # s loop

  # Get survey age neffs
  obs_srv_age_Neff <- array(srv_Neff[Fish_Start_yr[1]:(n_years - 1),], 
                             dim = c(length(years), n_srv_comps, n_sexes))
  

# Abundance Indices -------------------------------------------------------

  # Fishery index
  if(n_fish_indices == 1) { # if single fleet index, we need to do some munging
    # basically, we are going to take the average of the two fleets, if there are multiple 
    # fleets, but this can be genearlizable to a single fleet index, where it just retains
    # the original matrix structure
    obs_fish_indices <- matrix(nrow = length(years), ncol = 1)
    
    # Loop through matrix to average across rows
    for(i in 1:length(years))
      obs_fish_indices[i, 1] <- mean(Fishery_Index_Agg[(Fish_Start_yr[1]-1) + i,,sim])

  } else{
    obs_fish_indices <-as.matrix( Fishery_Index_Agg[Fish_Start_yr[1]:(n_years - 1),,sim], 
                                  nrow = length(years), ncol = n_fish_indices)
  } # multi fleet index = leave as is
  
  # Survey index
  obs_srv_indices <- as.matrix(Survey_Index_Agg[Fish_Start_yr[1]:(n_years - 1),,sim], 
                               nrow = length(years), ncol = n_srv_indices)
  
  
# Biological inputs -------------------------------------------------------

  # Weight at age
  WAA <- array(wt_at_age[Fish_Start_yr[1]:(n_years),,1,sim], 
               dim = c(length(Fish_Start_yr[1]:(n_years)), length(ages), n_sexes))
  
  # Maturity at age
  MatAA <- array(mat_at_age[Fish_Start_yr[1]:(n_years),,1,sim]  ,
                 dim = c(length(Fish_Start_yr[1]:(n_years)), length(ages), n_sexes))
  
  # Sex Ratios
  Sex_Ratio <- sex_ratio[1,]
  
# CV inputs ---------------------------------------------------------------
  
  # Fishery, survey, and catch CV
  fish_cv <- as.vector(fish_CV)
  srv_cv <- as.vector(srv_CV)
  catch_cv <- as.vector(catch_cv)

# Selectivity options ------------------------------------------------------
  
  F_Slx_model <- vector()
  S_Slx_model <- vector()
  
  # Fishery selex model input
  for(i in 1:n_fleets) {
    if(F_Slx_Model_Input[i] == "logistic") F_Slx_model[i] <- 0
    if(F_Slx_Model_Input[i] == "gamma") F_Slx_model[i] <- 1
    if(F_Slx_Model_Input[i] == "double_logistic") F_Slx_model[i] <- 2
    if(F_Slx_Model_Input[i] == "exp_logistic") F_Slx_model[i] <- 3
  } # end loop for selectivity model input for fishery
  
  # Survey selex model input
  for(i in 1:n_srv_comps) {
    if(S_Slx_Model_Input[i] == "logistic") S_Slx_model[i] <- 0
    if(S_Slx_Model_Input[i] == "gamma") S_Slx_model[i] <- 1
    if(S_Slx_Model_Input[i] == "double_logistic") S_Slx_model[i] <- 2
    if(S_Slx_Model_Input[i] == "exp_logistic") S_Slx_model[i] <- 3
  } # end loop for selectivity model input for survey
  
  # Specify selectivity blocks here
  F_Slx_Blocks <- F_Slx_Blocks_Input # fishery blocks
  S_Slx_Blocks <- S_Slx_Blocks_Input # survey blocks

# Data Indicators ---------------------------------------------------------
  
  # Catch data
  if(use_catch == TRUE) {
    use_catch <- matrix(1, nrow = length(years), ncol = n_fleets)
  } else{
    use_catch <- matrix(0, nrow = length(years), ncol = n_fleets)
  }
  
  # Fishery index of abundance
  if(use_fish_index == TRUE) {
    use_fish_index <- matrix(1, nrow = length(years), ncol = n_fish_indices)
  } else{
    use_fish_index <- matrix(0, nrow = length(years), ncol = n_fish_indices)
  }
  
  # Survey index of abundance
  if(use_srv_index == TRUE) {
    use_srv_index <- matrix(1, nrow = length(years), ncol = n_srv_indices)
  } else{
    use_srv_index <- matrix(0, nrow = length(years), ncol = n_srv_indices)
  }
  
  # Fishery comps
  if(use_fish_comps == TRUE) {
    use_fish_comps <- array(1, dim = c(length(years), n_fish_comps, n_sexes))
  } else{
    use_fish_comps <- array(0, dim = c(length(years), n_fish_comps, n_sexes))
  }
  
  # Survey comps
  if(use_srv_comps == TRUE) {
    use_srv_comps <- array(1, dim = c(length(years), n_srv_comps, n_sexes))
  } else{
    use_srv_comps <- array(0, dim = c(length(years), n_srv_comps, n_sexes))
  }

  # Input these data into a list object
  data$ages <- ages
  data$years <- years
  data$n_sexes <- n_sexes
  data$n_fleets = n_fleets
  data$n_fish_comps = n_fish_comps
  data$n_srv_comps = n_srv_comps
  data$n_fish_indices = n_fish_indices
  data$n_srv_indices = n_srv_indices
  data$obs_catches <- obs_catches
  data$obs_fish_age_comps <- obs_fish_age_comps
  data$obs_fish_age_Neff <- obs_fish_age_Neff
  data$obs_srv_age_comps <- obs_srv_age_comps
  data$obs_srv_age_Neff <- obs_srv_age_Neff
  data$obs_fish_indices <- obs_fish_indices
  data$obs_srv_indices <- obs_srv_indices
  data$WAA <- WAA
  data$MatAA <- MatAA
  data$Sex_Ratio <- Sex_Ratio
  data$fish_cv <- fish_cv
  data$srv_cv <- srv_cv
  data$catch_cv <- catch_cv
  data$F_Slx_Blocks <- F_Slx_Blocks
  data$S_Slx_Blocks <- S_Slx_Blocks
  data$use_catch <- use_catch
  data$use_fish_index <- use_fish_index
  data$use_srv_index  <- use_srv_index 
  data$use_fish_comps <- use_fish_comps
  data$use_srv_comps  <- use_srv_comps 
  data$rec_model <- rec_model
  data$S_Slx_model <- S_Slx_model
  data$F_Slx_model <- F_Slx_model
  

# Parameter specifications ------------------------------------------------

  # Set up parameters
  pars$ln_SigmaRec <- sigma_rec # recruitment variability
  pars$ln_RecDevs <- rnorm(length(years)-1, -1, 0.05) # rec devs
  pars$ln_N1_Devs <- rnorm(length(ages)-1, -1, 0.05) # intial recruitment deviaates
  pars$ln_M <- rnorm(1, 0, 0.1) # natural mortality
  pars$ln_Fy <- matrix(rnorm(n_fleets * length(years), -3, 0.05), 
                          ncol = n_fleets, nrow = length(years)) # fishing mortality
  pars$ln_q_fish <- rnorm(n_fish_indices, -1, 0.05) # catchability for fishery
  pars$ln_q_srv <- rnorm(n_srv_indices, -1, 0.05) # catchability for survey
  
  if(rec_model == 0) pars$ln_RecPars <- as.vector(rnorm(1, 0.1)) # Mean Recruitment (1 parameter)
  if(rec_model == 1) pars$ln_RecPars <- as.vector(rnorm(2, 0.1)) # BH Recruitment (2 parameters)
  
  # Selectivity parameters
  
  # Survey
  n_srv_blocks <- length(unique(as.vector(S_Slx_Blocks_Input))) # unique numbers (max surv blocks)
  
  # Now, figure out how many parameters we need to dimension the array
  # Create a vector to hold number of parameter values, and then take the max of that vector
  # to set up our array
  n_srv_pars <- vector()
  for(i in 1:length(S_Slx_Model_Input)) {
    if(S_Slx_Model_Input[i] == "logistic") n_srv_pars[i] <- 2
    if(S_Slx_Model_Input[i] == "gamma") n_srv_pars[i] <- 2
    if(S_Slx_Model_Input[i] == "exp_logistic") n_srv_pars[i] <- 3
    if(S_Slx_Model_Input[i] == "double_logistic") n_srv_pars[i] <- 4
  } # end i loop
  # Put array into our list
  pars$ln_srv_selpars <- array(rnorm(1, 0, 1), dim = c(n_srv_comps, n_sexes, n_srv_blocks, max(n_srv_pars)))
  
  # Do the same, but for the fishery
  n_fish_blocks <- length(unique(as.vector(F_Slx_Blocks_Input))) # unique numbers (max fish blocks)
  # Loop through to figure out how to dimension the last element of the array
  n_fish_pars <- vector()
  for(i in 1:length(F_Slx_Model_Input)) {
    if(F_Slx_Model_Input[i] == "logistic") n_fish_pars[i] <- 2
    if(F_Slx_Model_Input[i] == "gamma") n_fish_pars[i] <- 2
    if(F_Slx_Model_Input[i] == "exp_logistic") n_fish_pars[i] <- 3
    if(F_Slx_Model_Input[i] == "double_logistic") n_fish_pars[i] <- 4
  } # end i loop
  
  # put array into our parameter list
  pars$ln_fish_selpars <- array(log(0.5), dim = c(n_fish_comps, n_sexes, n_fish_blocks, max(n_fish_pars)))
  
  return(list(data = data, parameters = pars))
  
}
