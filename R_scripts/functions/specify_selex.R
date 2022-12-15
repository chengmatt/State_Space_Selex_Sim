# Purpose: To specify selectivity parameterizations for the fishery and the survey
# Creator: Matthew LH. Cheng
# Date: 10/30/22

# Selectivity options -----------------------------------------------------

#' @param selex_type Selectivity type (Options are: uniform, logistic, gamma, double_logistic, double_normal)
#' @param bins vector of age or length bins
#' @param par_values parameter values specific for each selectivity type
#' # This fxn is meant to be used within the specify_selex fxn.

selex_opts <- function(selex_type, bins, par_values) {
  
  if(selex_type != "uniform") {
    
    # List out all the different parameters there are
    par_vec <- c("a50", "k", "amax", "delta", "slope_1", "slope_2", "infl_1", "infl_2",
                 "p1", "p2", "p3", "p4", "p5", "p6")
    
    # Create a matrix to assign values to - input parameter names in the 2nd col and nas (values) in the 1st
    par_mat <- matrix(data = c(rep(NA, length(par_vec)), par_vec), ncol = 2, nrow = length(par_vec))
    
    # Get parameter names specified in the array
    par_names <- names(par_values)
    
    # Put parameter values in the corresponding area
    par_mat[,1][par_mat[,2] %in% c(par_names)] <- par_values
    
  } # only look for par values if selectivity is not uniform

  if(selex_type == "logistic") { # 2 parameters
    
    # Get a50 value
    a50 <- as.numeric(par_mat[,1][par_mat[,2] == "a50"])
    # Get k value
    k <- as.numeric(par_mat[,1][par_mat[,2] == "k"])
    # Compute selex and scale to max of 1
    selex <- 1 / (1 + exp(-1 * (bins - a50) / k)) / max( 1 / (1 + exp(-1 * (bins - a50) / k)) ) / max(1 / (1 + exp(-1 * (bins - a50) / k)) / max( 1 / (1 + exp(-1 * (bins - a50) / k)) ))
    
  } # logistic selectivity
  
  if(selex_type == "uniform") { # 0 parameters 
    selex <- rep(1, length(bins))
  } # uniform selectivity
  
  if(selex_type == "gamma") { # 2 parameters 
    
    # Get delta value
    delta <- as.numeric(par_mat[,1][par_mat[,2] == "delta"])
    # Get amax value
    amax <- as.numeric(par_mat[,1][par_mat[,2] == "amax"])
    
    # Base Fucntional Form here - use Punt et al. 1994 parameterization (amax and delta)
    p <- (0.5 * (sqrt(amax^2 + 4*delta^2) - amax)) # Derive power parameter here first
    
    # Now calculate Selex for dome-shaped - scale to a max of 1
    selex <- (bins/amax) ^ (amax/p) * exp((amax - bins) / p) / 
                  max( (bins/amax) ^ (amax/p) * exp((amax - bins) / p) )
    
  } # gamma dome-shaped selectivity
  
  if(selex_type == "double_logistic") { # 4 parameters
    
    # TESTING 
    # slope_1 <- 3 # ascending limb
    # slope_2 <- 0.2 # descending limb
    # infl_1 <- 5 # inflection point for ascending limb
    # infl_2 <- 20 # inflection point for descending limb
    
    # Get slope1 value
    slope_1 <- as.numeric(par_mat[,1][par_mat[,2] == "slope_1"])
    # Get slope2 value
    slope_2 <- as.numeric(par_mat[,1][par_mat[,2] == "slope_2"])
    # Get inflection point 1 value
    infl_1 <- as.numeric(par_mat[,1][par_mat[,2] == "infl_1"])
    # Get inflection point 2 value
    infl_2 <- as.numeric(par_mat[,1][par_mat[,2] == "infl_2"])
    
    # Calculate logistic curve 1 - mediates the ascending limb
    logistic_1 <- 1/(1 + exp(-slope_1 * (bins - infl_1))) 
    
    # Calculate logistic curve 2 - mediates the descending limb
    logistic_2 <- 1 - (1/(1 + exp(-slope_2 * (bins - infl_2))))
    
    # multiply the logistic curves and scale to a max of 1
    selex <- logistic_1 * logistic_2 / max(logistic_1 * logistic_2)
    
    # plot(selex)
    
  } # double logistic selectivity
  
  if(selex_type == "double_normal") { # 6 parameters
    
    # p1 = when age at selex = 1 starts (should probably bound by min and max age)
    # p2 = partially controls where selex ends - should be in log space
    # p3 = slope of ascending section (higher values = less steep ascending)
    # p4 = slope of descending section (higher values = less steep decline)
    # p5 = selex at age 0/young ages (needs to be in inv logit space or between 0 and 1)
    # p6 = if set at 1 - asymptotic (where the right half of the curve limb ends, should also be in inv logit space)
    
    # TESTING
    bins <- 2:30
    p1 <- 5
    p2 <- 0.01
    p3 <- 3
    p4 <- 10
    p5 <- 0.2
    p6 <- 0.1
    
    # Get input parameter values
    p1 <- as.numeric(par_mat[,1][par_mat[,2] == "p1"]) # p1
    p2 <- as.numeric(par_mat[,1][par_mat[,2] == "p2"]) # p2
    p3 <- as.numeric(par_mat[,1][par_mat[,2] == "p3"]) # p3
    p4 <- as.numeric(par_mat[,1][par_mat[,2] == "p4"]) # p4
    p5 <- as.numeric(par_mat[,1][par_mat[,2] == "p5"]) # p5
    p6 <- as.numeric(par_mat[,1][par_mat[,2] == "p6"]) # p6
 
    # Derive gamma (age at which selex = 1 ends)
    gamma <- p1 + 1 + (0.99 * max(bins) - p1 - 1) / (1 + exp(-p2))
    
    # Derive joiner function 1 (ascending) mimics a logistic-like function but w/ knife-edged dynamics
    j1 <- 1 / (1 + exp( -1 * (20 / (1 + abs(bins - p1))) * (bins - p1) ))
    
    # Derive joiner function 2 (descending)
    j2 <- 1 / (1 + exp( -1 * (20 / (1 + abs(bins - gamma))) * (bins - gamma) ))

    # Derive alpha parameter vector
    alpha <- p5 + (1 - p5) * (exp( (-1*((bins - p1)^2)) / exp(p3) ) - exp(-1*((1 - p1)^2)/ exp(p3) )) / (1 - exp(-1*((1 - p1)^2)/ exp(p3)))
    
    # Derive beta parameter vector
    beta <- 1 + (p6 - 1) * (exp( (-((bins - gamma)^2)) / exp(p4) ) -1 )/ (exp(-((max(bins) - gamma)^2)/exp(p4)) - 1)
    
    # Now, get selectivity values - scale to a max of 1
    selex <- alpha * (1 - j1) + (j1 * ( (1-j2) + j2 * beta)) 
    
    plot(selex)
    
  } # double normal selectivity
  
  return(selex)
  
} # end selectivity options function




# Specify selectivity -----------------------------------------------------
  
#' @param fish_selex selectivity types for fishery selectivity (should be a vector if > 1 fleet) 
#' # (Options are: uniform, logistic, gamma, double_logistic, double_normal)
#' @param srv_selex selectivity types for survey selectivity (should be a vector if > 1 fleet) 
#' # (Options are: uniform, logistic, gamma, double_logistic, double_normal) 
#' @param fish_pars parameters for the respective fishery fleets and selectivity forms (needs to be a list)
#' if we are specifying uniform selectivity, the matrix can be specified just as matrix()
#' @param srv_pars parameters for the respective survey fleets and selectivity forms (needs to be a list of matrices with col length equal to
#' the number of parameters and row length equal to the number of rows)
#' @param bins age of length bins

#' @examples fish_pars <- list(Fleet_1_l <- matrix(data = c(2,3,
#'                             4,3), nrow = n_sex, byrow = TRUE),
#' Fleet_2_dn <- matrix(data = c(5,0.01,3,10,0.2,0.1,
#'                              10,0.01,3,15,0.2,0.05), nrow = n_sex, byrow = TRUE))

#'srv_pars <- list(Fleet_1_l <- matrix(data = c(2,3,
#'                               4,3), nrow = n_sex, byrow = TRUE),
#'                               Fleet_2_u <- matrix())

specify_selex <- function(fish_selex, srv_selex, fish_pars, srv_pars, bins) {
  
  if(length(fish_selex) != n_fish_fleets) stop("Fishery selex options are not of equal length to number of fishery fleets specified!")
  if(length(srv_selex) != n_srv_fleets) stop("Survey selex options are not of equal length to number of survey fleets specified!")
  
  # Create objects to loop through
  fish_par_list <- list()
  srv_par_list <- list()

# Fishery Selex -----------------------------------------------------------

  # Loop through to construct objects to use within the function
  for(f in 1:n_fish_fleets) {
    
    if(nrow(fish_pars[[f]]) != n_sex) stop("The number of sexes in the list of matrices (fishery) does not equal the number of sexes specified")
    
      # Specify number of cols needed for constructing our matrix using npars
      if(fish_selex[f] == "uniform") {n_pars <- 0; par_names <- NULL}
      if(fish_selex[f] == "logistic") {n_pars <- 2; par_names <- c("a50", "k")
      if(ncol(fish_pars[[f]]) != n_pars) stop("Number of parameters specified does not equal number of parameters required for a logistic (2)!")}
      if(fish_selex[f] == "gamma") {n_pars <- 2; par_names <- c("amax", "delta")
      if(ncol(fish_pars[[f]]) != n_pars) stop("Number of parameters specified does not equal number of parameters required for a gamma (2)!")}
      if(fish_selex[f] == "double_logistic") {n_pars <- 4; par_names <- c("slope_1", "slope_2", "infl_1", "infl_2")
      if(ncol(fish_pars[[f]]) != n_pars) stop("Number of parameters specified does not equal number of parameters required for a double logistic (4)!")}
      if(fish_selex[f] == "double_normal") {n_pars <- 6; par_names <- c("p1", 'p2', "p3", "p4", "p5", "p6") 
      if(ncol(fish_pars[[f]]) != n_pars) stop("Number of parameters specified does not equal number of parameters required for a double normal (6)!")}
      
      # Create an empty named matrix to store specified values in
      fish_mat <- matrix(data = fish_pars[[f]], nrow = n_sex, ncol = n_pars) # fIrst row = female, second row = male if 2 sexes
      # Specify column names based on the selectivity type
      colnames(fish_mat) <- par_names
      
      fish_par_list[[f]] <- fish_mat 
      
  } # end f loop for fishery fleets
  
  # Loop through to assign values
  for(y in 1:n_years) {
    
    for(f in 1:n_fish_fleets) {
      
      for(s in 1:n_sex) {
        
      # Selectivity type within a fleet across sexes remain the same
      Fish_selex_at_age[y,,f,s,] <- selex_opts(selex_type = fish_selex[f], bins = bins, 
                                           par_values = fish_par_list[[f]][s,]) 
      # Note that fish_par_list is a list so we need to do some weird indexing 
        
      } # end s loop with number of sexes
      
    } # end f loop - number of fishery fleets
    
  } # end y loop years
  

# Survey Selex ------------------------------------------------------------

  
  # Loop through to construct objects to use within the function
  for(sf in 1:n_srv_fleets) {
    
    if(nrow(srv_pars[[sf]]) != n_sex) stop("The number of sexes in the list of matrices (survey) does not equal the number of sexes specified")
    
    # Specify number of cols needed for constructing our matrix using npars - adding some stops to make sure things are being specified correctly.
    if(srv_selex[sf] == "uniform") {n_pars <- 0; par_names <- NULL}
    if(srv_selex[sf] == "logistic") {n_pars <- 2; par_names <- c("a50", "k")
    if(ncol(srv_pars[[sf]]) != n_pars) stop("Number of parameters specified does not equal number of parameters required for a logistic (2)!")}
    if(srv_selex[sf] == "gamma") {n_pars <- 2; par_names <- c("amax", "delta")
    if(ncol(srv_pars[[sf]]) != n_pars) stop("Number of parameters specified does not equal number of parameters required for a gamma (2)!")}
    if(srv_selex[sf] == "double_logistic") {n_pars <- 4; par_names <- c("slope_1", "slope_2", "infl_1", "infl_2")
    if(ncol(srv_pars[[sf]]) != n_pars) stop("Number of parameters specified does not equal number of parameters required for a double logistic (4)!")}
    if(srv_selex[sf] == "double_normal") {n_pars <- 6; par_names <- c("p1", 'p2', "p3", "p4", "p5", "p6") 
    if(ncol(srv_pars[[sf]]) != n_pars) stop("Number of parameters specified does not equal number of parameters required for a double normal (6)!")}
    
    # Create an empty named matrix to store specified values in
    srv_mat <- matrix(data = srv_pars[[sf]], nrow = n_sex, ncol = n_pars) # fIrst row = female, second row = male if 2 sexes
    # Specify column names based on the selectivity type
    colnames(srv_mat) <- par_names
    
    srv_par_list[[sf]] <- srv_mat # stick into list
    
  } # end f loop for fishery fleets
  
  # Loop through to assign values
  for(y in 1:n_years) {
    
    for(sf in 1:n_srv_fleets) {
      
      for(s in 1:n_sex) {
        
        # Selectivity type within a fleet across sexes remain the same
        Surv_selex_at_age[y,,sf,s,] <- selex_opts(selex_type = srv_selex[sf], bins = bins, 
                                                 par_values = srv_par_list[[sf]][s,]) 
        # Note that fish_par_list is a list so we need to do some weird indexing 
        
      } # end s loop with number of sexes
      
    } # end sf loop - number of survey fleets
    
  } # end y loop years  

  Surv_selex_at_age <<- Surv_selex_at_age # Output into environment
  Fish_selex_at_age <<- Fish_selex_at_age # Output into environment
  
} # end function here
