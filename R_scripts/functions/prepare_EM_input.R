# Creator: Matthew LH. Cheng
# Purpose: To preapre inputs for the estimation model in TMB
# Date 1/10/23

#' Title Prepares data inputs for the EM 
#'
#' @param years vector of years
#' @param ages vector of ages
#' @param catch_cv vector of catch CVs
#' @param F_Slx_Blocks_Input matrix of when you want to block
#' @param S_Slx_Blocks_Input same as above
#' @param use_catch Boolean (TRUE OR FALSE) TRUE = Use, FALSE = don't use
#' @param use_fish_index same as above
#' @param use_srv_index same as above
#' @param use_fish_comps Boolean (TRUE OR FALSE) TRUE = USE, FALSE = don't use
#' @param use_srv_comps same above abovce
#' @param Fish_Comp_Like_Model Fishery composition likelihoods ("multinomial" or "dirichlet_multinomial")
#' @param Srv_Comp_Like_Model Survey composition likelihoods ("multinomial" or "dirichlet_multinomial")
#' @param rec_model recruitment model == 0 (mean recruitment)
#' @param F_Slx_Model_Input Fishery selectivity model character (logistic, gamma, double_logistic, exp_logistic). Needs to be a matrix dimensioned by
#' bn_years x n_fleets
#' @param S_Slx_Model_Input same as above
#' @param time_selex Type of time-varying to do (Options are: None, RW)
#' @param fix_pars What parameters we want to fix
#' @param n_time_selex_pars Number of time-varying selectivity parameters we want to estimate on the parametric form
#' @param sim simulation indexing
#'
#' @return
#' @export
#'
#' @examples
#' Note that n_years here refers to all the rows of the array
prepare_EM_input <- function(years,
                             ages,
                             Fish_Start_yr = 1,
                             n_fleets, 
                             catch_cv,
                             F_Slx_Blocks_Input,
                             S_Slx_Blocks_Input,
                             use_catch = TRUE,
                             use_fish_index = TRUE,
                             use_srv_index = TRUE, 
                             use_fish_comps = TRUE,
                             use_srv_comps = TRUE,
                             Fish_Comp_Like_Model,
                             Srv_Comp_Like_Model,
                             F_Slx_Model_Input,
                             S_Slx_Model_Input,
                             rec_model = 0,
                             time_selex = "None",
                             n_time_selex_pars = 1,
                             fix_pars = NA,
                             share_ages = 1,
                             sim) {
  
  # Make data into a list 
  data <- list()
  # Make parameters into a list
  pars <- list()
  # Mapping to fix certain parameters
  map <- list()
  
  # Specify dimensions 
  ages <- ages # ages
  Fish_Start_yr <- Fish_Start_yr # Fishery start year
  n_sexes <- n_sex # sexes
  n_years <- length(years) # number of years
  n_fleets <- n_fleets # number of fleets
  n_fish_comps <- n_fleets # number of fleets = number of fish comps
  n_fish_indices <- n_fleets # number of fish indices = number of fleets
  n_srv_comps <- dim(Surv_selex_at_age)[3] # 3rd dimension of this array = number of survey fleets
  n_srv_indices <- dim(Survey_Index_Agg)[2] # 2nd dimension of this array = number of survey fleets
  
  
  # Catch -------------------------------------------------------------------
  
  if(n_fleets == 1) { # if single fleet - sum across or leave as is
    obs_catches <- as.matrix(apply( as.matrix(Catch_agg[(Fish_Start_yr[1]:(n_years)),,sim]), 1, FUN = sum),
                             ncol = 1)
    
    # Catch weighting here
    catch_sum <- rowSums(as.matrix(Catch_agg[(Fish_Start_yr[1]:(n_years)),,sim]))
    catch_weight <- as.matrix(Catch_agg[(Fish_Start_yr[1]:(n_years)),,sim]) / catch_sum
    
  } else{ # multi fleet - leave as is
    obs_catches <- matrix(Catch_agg[(Fish_Start_yr[1]:(n_years)),,sim], nrow = length(years), ncol = n_fleets)
  } # end else
  
  # Fishery Age Comps -------------------------------------------------------
  
  obs_fish_age_comps <- array(data = 0, dim = c(length(years), length(ages), n_fleets, n_sexes))
  
  if(n_fleets == 1) { # one fleet = sum across
    
    # Loop through to sum - index the 3rd dimension to get fleets
    for(s in 1:n_sexes) {
      for(f in 1:dim(Fish_Age_Comps)[3]) {
        # Filter to save as an object
        fish_age_comps <-  Fish_Age_Comps[Fish_Start_yr[1]:length(years),,f,s,sim] * catch_weight[,f]
        # Increment comps - fixing fleet index to 1 here
        obs_fish_age_comps[,,1,s] <- obs_fish_age_comps[,,1,s] + fish_age_comps
      } # end f loop
    } # end s loop
    
    # Empty array to store values in
    obs_fish_age_Input_N <- array(data = NA, dim = c(length(years), 1, n_sexes))
    
    if(Fish_Comp_Like_Model %in% c("dirichlet_multionimal", "multinomial")) { # if mutlinomial or DM
      # Effective sample size
      for(s in 1:n_sexes) {
        # obs_fish_age_Input_N[,1,s] <- floor(rowSums(obs_fish_age_comps[,,,s]))
        obs_fish_age_Input_N[,1,s] <- round(rowSums(Input_N_Fish[Fish_Start_yr[1]:length(years),] * catch_weight))
      } # s loop
    } else{ # Dirichlet input sample sizes - weigh by the catch and input N (bc simulation sums to 1)
      # Effective sample size
      for(s in 1:n_sexes) {
        obs_fish_age_Input_N[,1,s] <- rowSums(matrix(Input_N_Fish[Fish_Start_yr[1]:length(years),] *
                                                       catch_weight, ncol = dim(Fish_Age_Comps)[3]))
      } # s loop
    } # end else statement
    
    # Now, apply the proportion function over a single fleet
    for(s in 1:n_sexes) {
      obs_fish_age_comps[,,,s] <- t(apply(obs_fish_age_comps[,,,s], MARGIN = 1, 
                                          FUN=function(x) { x/sum(x) }))
    } # s loop
    
  } else{ # more than one fleet here
    
    for(f in 1:n_fleets) { # needs to loop through transpose and apply for easy retnetion of array dimensions
      for(s in 1:n_sexes) {
        obs_fish_age_comps[,,f,s] <- t(
          apply(Fish_Age_Comps[Fish_Start_yr[1]:length(years),,f,s,sim], MARGIN = 1, 
                FUN=function(x) { x/sum(x) })
        )
      } # end s loop
    } # end f loop
    
    # Effective Sample Sizes
    obs_fish_age_Input_N <- array(Input_N_Fish[Fish_Start_yr[1]:length(years),], 
                                  dim = c(length(years), n_fish_comps, n_sexes))
    
  } # end else
  
  
  # Survey Age Comps --------------------------------------------------------
  # Right now, this is only able to accomdate one single survey fleet
  
  obs_srv_age_comps <- array(data = NA, dim = c(length(years), length(ages), n_srv_comps, n_sexes))
  
  # Apply function to get proportions and munge into matrix
  # Now, apply the proportion function over a single fleet
  for(s in 1:n_sexes) {
    obs_srv_age_comps[,,,s] <- t(apply(Survey_Age_Comps[Fish_Start_yr[1]:length(years),,,s,sim],
                                       MARGIN = 1, FUN=function(x) { x/sum(x) }))
  } # s loop
  
  # Get survey age input sample size
  obs_srv_age_Input_N <- array(Input_N_Srv[Fish_Start_yr[1]:length(years),], 
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
    obs_fish_indices <-as.matrix( Fishery_Index_Agg[Fish_Start_yr[1]:length(years),,sim], 
                                  nrow = length(years), ncol = n_fish_indices)
  } # multi fleet index = leave as is
  
  # Survey index
  obs_srv_indices <- as.matrix(Survey_Index_Agg[Fish_Start_yr[1]:length(years),,sim], 
                               nrow = length(years), ncol = n_srv_indices)
  
  
  # Biological inputs -------------------------------------------------------
  
  WAA <- array(data = NA, dim = c(length(Fish_Start_yr[1]:(n_years)), 
                                  length(ages), n_sexes))
  # Weight at age
  for(s in 1:n_sexes) {
    WAA[,,s] <- wt_at_age[Fish_Start_yr[1]:(n_years),,s,sim]
  }
  
  # Maturity at age
  MatAA <- array(mat_at_age[Fish_Start_yr[1]:(n_years),,1,sim],
                 dim = c(length(Fish_Start_yr[1]:(n_years)), 
                         length(ages), n_sexes))
  
  # Sex Ratios
  Sex_Ratio <- sex_ratio[1,]
  
  # CV inputs ---------------------------------------------------------------
  
  # Fishery, survey, and catch CV
  fish_cv <- as.vector(fish_CV)
  srv_cv <- as.vector(srv_CV)
  catch_cv <- as.vector(catch_cv)
  
  # Selectivity options ------------------------------------------------------
  
  # For the fishery selectivity models, we want this as a matrix so we can time-block and switch
  # over to a different selectivity form. 
  F_Slx_model <- matrix(data = NA, nrow = length(years), ncol = n_fleets)
  S_Slx_model <- vector()
  
  # Fishery selex model input
  for(y in 1:n_years) {
    for(i in 1:n_fleets) {
      if(F_Slx_Model_Input[y,i] == "logistic") F_Slx_model[y, i] <- 0
      if(F_Slx_Model_Input[y,i] == "gamma") F_Slx_model[y, i] <- 1
      if(F_Slx_Model_Input[y,i] == "double_logistic") F_Slx_model[y, i] <- 2
      if(F_Slx_Model_Input[y,i] == "exp_logistic") F_Slx_model[y, i] <- 3
    } # end i loop (fleets)
  } # end y loop (years)
  
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
  
  # Compositional Likelihoods
  # Fishery
  fish_comp_likelihoods <- matrix(ncol = n_sexes, nrow = n_fish_comps)
  ln_DM_Fish_Param <- matrix(ncol = n_sexes, nrow = n_fish_comps) # set up ln fishery parameters
  map_ln_DM_Fish_Param <- matrix(ncol = n_sexes, nrow = n_fish_comps) # set up mapping for DM 
  
  for(i in 1:length(Fish_Comp_Like_Model)) {
    
    if(Fish_Comp_Like_Model[i] == "multinomial") {
      fish_comp_likelihoods[i,1:n_sexes] <- rep(0, n_sexes)
      map_ln_DM_Fish_Param[i,1:n_sexes] <- rep(NA, n_sexes)
    } # if Multinomial
    
    if(Fish_Comp_Like_Model[i] == "dirichlet_multinomial") {
      fish_comp_likelihoods[i, 1:n_sexes] <- rep(1, n_sexes)
      map_ln_DM_Fish_Param[i,1:n_sexes] <- rep(1, n_sexes)
    } # if DM 
    
    if(Fish_Comp_Like_Model[i] == "dirichlet") {
      fish_comp_likelihoods[i, 1:n_sexes] <- rep(2, n_sexes)
      map_ln_DM_Fish_Param[i,1:n_sexes] <- rep(1, n_sexes)
    } # if Dirichlet 
    
    # Random starting values for DM theta parameter
    ln_DM_Fish_Param[i,1:n_sexes] <- log(3) 
    
  } # end i fishery
  
  # Survey
  srv_comp_likelihoods <- matrix(ncol = n_sexes, nrow = n_srv_comps)
  ln_DM_Srv_Param <- matrix(ncol = n_sexes, nrow = n_srv_comps)
  map_ln_DM_Srv_Param <- matrix(ncol = n_sexes, nrow = n_srv_comps)
  for(i in 1:length(Srv_Comp_Like_Model)) {
    
    if(Srv_Comp_Like_Model[i] == "multinomial") {
      srv_comp_likelihoods[i,1:n_sexes] <- rep(0, n_sexes)
      map_ln_DM_Srv_Param[i,1:n_sexes] <- rep(NA, n_sexes)
    } # if multinomial
    
    if(Srv_Comp_Like_Model[i] == "dirichlet_multinomial") {
      srv_comp_likelihoods[i,1:n_sexes] <- rep(1, n_sexes)
      map_ln_DM_Srv_Param[i,1:n_sexes] <- rep(1, n_sexes)
    } # if DM
    
    if(Srv_Comp_Like_Model[i] == "dirichlet") {
      srv_comp_likelihoods[i,1:n_sexes] <- rep(2, n_sexes)
      map_ln_DM_Srv_Param[i,1:n_sexes] <- rep(1, n_sexes)
    } # if dirichlet
    
    # Random starting values for DM theta parameter
    ln_DM_Srv_Param[i,1:n_sexes] <- log(3) 
    
  } # end i survey
  
  # Map factor for fishery and survey comps if there is a multinomial likelihood
  if(Fish_Comp_Like_Model %in% "multinomial") map$ln_DM_Fish_Param <- factor(map_ln_DM_Fish_Param)
  if(Srv_Comp_Like_Model %in% "multinomial")  map$ln_DM_Srv_Param <- factor(map_ln_DM_Srv_Param)
  
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
  data$obs_fish_age_Input_N <- obs_fish_age_Input_N
  data$obs_srv_age_comps <- obs_srv_age_comps
  data$obs_srv_age_Input_N <- obs_srv_age_Input_N
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
  data$fish_comp_likelihoods <- fish_comp_likelihoods
  data$srv_comp_likelihoods <- srv_comp_likelihoods
  data$S_Slx_model <- S_Slx_model
  data$F_Slx_model <- F_Slx_model
  
  
  # Parameter specifications ------------------------------------------------
  
  # Set up parameters
  pars$ln_DM_Fish_Param <- ln_DM_Fish_Param # DM theta for fishery
  pars$ln_DM_Srv_Param <- ln_DM_Srv_Param # DM theta for survey
  pars$ln_SigmaRec <- log(sigma_rec) # recruitment variability
  pars$ln_RecDevs <- rec_devs[years[-1],sim] * 0.3 # rec devs
  pars$ln_N1Devs <- rnorm(length(ages)-1, 0.1, 0.0) # rec devs
  if(sum(fix_pars %in% c("ln_M")) == 1) {
    pars$ln_M <- log(mean(Mort_at_age)) # natural mortality
  } else{
    pars$ln_M <- log(mean(Mort_at_age)) * 0.3 # natural mortality
  }
  
  # row sums for fish mort
  if(n_fleets == 1 & dim(Fish_Age_Comps)[3] == 1) { # 1 fleet model and truth = multiple
    pars$ln_Fy <- matrix(log(fish_mort[Fish_Start_yr[1]:length(years),,sim] * 0.3),
                         ncol = n_fleets, nrow = length(years))
  }
  
  if(n_fleets == 1 & dim(Fish_Age_Comps)[3] > 1) { # 1 fleet model, but truth = multiple
    pars$ln_Fy <- matrix(log(rowSums(fish_mort[Fish_Start_yr[1]:length(years),,sim]) * 0.3),
                         ncol = n_fleets, nrow = length(years))
  }
  
  if(n_fleets > 1 & dim(Fish_Age_Comps)[3] > 1) { # more tha one fleet for both OM and EM
    pars$ln_Fy <- matrix(log(fish_mort[Fish_Start_yr[1]:length(years),,sim] * 0.3),
                         ncol = n_fleets, nrow = length(years))
  }
  
  pars$ln_q_fish <- log(rnorm(n_fish_indices, 0.1, 0.0)) # catchability for fishery
  pars$ln_q_srv <- log(rnorm(n_srv_indices, 0.1, 0.0)) # catchability for survey
  
  if(rec_model == "mean_rec") {
    pars$ln_RecPars <- as.vector(c(mu_rec) * 0.3) # Mean Recruitment (1 parameter)
    data$rec_model <- 0
  }
  if(rec_model == "BH") {
    pars$ln_RecPars <- as.vector(c(log(r0) * 0.3, log(h))) # R0 and steepness
    data$rec_model <- 1
  }
  
  # Survey
  n_srv_blocks <- length(unique(as.vector(S_Slx_Blocks_Input))) # unique numbers (max surv blocks)
  
  # Now, figure out how many parameters we need to dimension the array
  # Create a vector to hold number of parameter values, and then take the max of that vector
  # to set up our array
  n_srv_pars <- vector()
  for(i in 1:length(S_Slx_Model_Input)) {
    if(S_Slx_Model_Input[i] == "logistic") n_srv_pars[i] <- 2
    if(S_Slx_Model_Input[i] == "gamma") n_srv_pars[i] <- 2
    if(S_Slx_Model_Input[i] == "normal") n_srv_pars[i] <- 2
    if(S_Slx_Model_Input[i] == "exp_logistic") n_srv_pars[i] <- 3
    if(S_Slx_Model_Input[i] == "double_logistic") n_srv_pars[i] <- 4
  } # end i loop
  # Put array into our list
  pars$ln_srv_selpars <- array(log(rnorm(1, 4, 0)), dim = c(n_srv_comps, n_sexes, n_srv_blocks, max(n_srv_pars)))
  
  # Do the same, but for the fishery
  n_fish_blocks <- length(unique(as.vector(F_Slx_Blocks_Input))) # unique numbers (max fish blocks)
  # Loop through to figure out how to dimension the last element of the array
  n_fish_pars <- vector()
  for(i in 1:length(F_Slx_Model_Input)) {
    if(F_Slx_Model_Input[i] == "logistic") n_fish_pars[i] <- 2
    if(F_Slx_Model_Input[i] == "gamma") n_fish_pars[i] <- 2
    if(F_Slx_Model_Input[i] == "normal") n_fish_pars[i] <- 2
    if(F_Slx_Model_Input[i] == "exp_logistic") n_fish_pars[i] <- 3
    if(F_Slx_Model_Input[i] == "double_logistic") n_fish_pars[i] <- 4
  } # end i loop
  
  # put array into our parameter list
  # need to do slightly different starting values here, b/c exp_logistic is a lil finicky
  if(sum(reshape2::melt(F_Slx_Model_Input)$value == "exp_logistic") == 0) { 
    # if we don't have exponential logistic or double logistic
    pars$ln_fish_selpars <- array(log(rnorm(1, 3.5, 0)), 
                                  dim = c(n_fish_comps, n_sexes, 
                                          n_fish_blocks, max(n_fish_pars)))
  } else{
    pars$ln_fish_selpars <- array(0, dim = c(n_fish_comps, n_sexes, 
                                             n_fish_blocks, max(n_fish_pars)))
    
    pars$ln_fish_selpars[,,,1] <- log(0.3) # gamma parameter - degree of doming
    pars$ln_fish_selpars[,,,2] <- log(7) # beta parameter - controls peak of ascending limb
    pars$ln_fish_selpars[,,,3] <- log(0.3) # alpha parameter - moves ascending limb back (steepness of ascending limb)
  }  # if we have an exponential logistic
  
  # Time-Varying Selectivity Options (Fishery)
  if(time_selex == "None") { # No time-varying
    data$F_Slx_re_model <- matrix(1e5, nrow = n_fish_comps, ncol = n_sexes)
    pars$ln_fish_selpars_re <- log(array(rnorm(1, 0.05, 0),
                                         dim = c((length(years) - 1), n_time_selex_pars, n_fish_comps, n_sexes)))
    pars$ln_fixed_sel_re_fish <- array(rnorm(1,0.3,0), dim = c(1, n_fish_comps, n_sexes))
    
    # Set mapping here for random effects (no random effects)
    map$ln_fish_selpars_re <- factor(rep(NA, length(pars$ln_fish_selpars_re))) # Fix random effects
    map$ln_fixed_sel_re_fish <- factor(rep(NA, length(pars$ln_fixed_sel_re_fish))) # fixed deivation parameters/don't estimate
    
  } # none if statement
  
  if(time_selex == "RW") { # Random Walk
    data$F_Slx_re_model <- matrix(0, nrow = n_fish_comps, ncol = n_sexes)
    pars$ln_fish_selpars_re <- array(0,
                                     dim = c((length(years)-1), n_time_selex_pars, n_fish_comps, n_sexes))
    pars$ln_fixed_sel_re_fish <- array(log(rnorm(1,1,0)), dim = c(n_time_selex_pars, n_fish_comps, n_sexes))
  } # random walk if statement
  
  if(time_selex == "Semi-Parametric") {
    data$F_Slx_re_model = matrix(1, nrow = n_fish_comps, ncol = n_sexes)
    pars$ln_fish_selpars_re = array(0, dim = c(length(years)-1, length(ages), n_fish_comps, n_sexes))
    pars$ln_fixed_sel_re_fish = array(log(rnorm(1,1,0)), dim = c(n_fish_comps, n_sexes))

    # Share random effects for adjacent ages
    map_ln_fish_selpars_re = pars$ln_fish_selpars_re # create mapping variable here
    unique_idx = (prod(dim(pars$ln_fish_selpars_re)) / share_ages) # get number of unique elements to mape
    idx = rep(seq(1:unique_idx), each = share_ages) # create index variable to loop through
    counter = 1 # counter for looping
    # loop through indexing now
    for(y in 1:(length(years)-1)) {
      for(f in 1:n_fish_comps) {
        for(s in 1:n_sexes) {
          for(a in 1:length(ages)) {
            map_ln_fish_selpars_re[y,a,f,s] = idx[counter]
            counter = counter + 1 # update counter
          } # end a loop
        } # end s loop
      } # end f loop
    } # end y loop
    
    # update mapping
    map_ln_fish_selpars_re = factor(as.numeric(map_ln_fish_selpars_re))
    map$ln_fish_selpars_re = map_ln_fish_selpars_re
  } # end if semi-parametric
  
  # Parameter mapping -------------------------------------------------------
  if(sum(fix_pars %in% c("ln_h")) == 1) {
    map$ln_RecPars <- factor(c(1, NA)) # fixing steepness
    # Remove steepness from fix_pars vector so it goes through the next loop properly
    fix_pars <- fix_pars[fix_pars != 'ln_h']
  }
  # Loop through to map parameters
  for(i in 1:length(fix_pars)) {
    
    # Get parameter length here
    par_length <- length(unlist(pars[names(pars) == fix_pars[i]]))
    
    # Now, stick the map parameter into a list
    map_par <- list( factor(rep(NA, par_length)) )
    names(map_par) <- fix_pars[i] # name the list
    
    # Now, append this to our map list
    map <- c(map_par, map)
    
  } # end i
  
  # Parameter mapping for time-blocking purposes in fishery selectivity here
  if(length(unique(n_fish_pars)) > 1) { # only do this if there is more than 2 unique selex forms w/ different params
    
    # Get dummy selex array we can manipulate
    fish_selpars_array <- pars$ln_fish_selpars
    # Get unique number of fishery selex pars
    unique_fish_pars <- unique(n_fish_pars)
    
    # Figure out which block is the missing parameter here
    blk_missing <- which(unique_fish_pars != max(n_fish_pars))
    # Now replace this with a unique identifier so we know to
    # put an NA in the parameter mapping here.
    fish_selpars_array[,,blk_missing, max(n_fish_pars)] <- NA
    
    # now, vectorize this 
    fish_selpars_vec <- as.vector(fish_selpars_array)
    
    # now, map this out
    map_fish_selpars <- vector()
    counter <- 1
    for(i in 1:length(fish_selpars_vec)) {
      # If this is an NA, input this as an NA
      if(is.na(fish_selpars_vec[i])) map_fish_selpars[i] <- NA
      else {
        # input counter into vector here to identify parameter to estimate
        map_fish_selpars[i] <- counter
        counter <- counter + 1
      } # else statement
    } # end i loop
    
    # Turn into factor and input into mapping list
    ln_fish_selpars <- factor(map_fish_selpars)
    map <- rlist::list.append(map, ln_fish_selpars = ln_fish_selpars)
    
  } # time block mapping
  
  return(list(data = data, parameters = pars, map = map))
  
} # end function
