# Purpose: Utility functions for running TMB models
# Creator: Matthew LH. Cheng
# Date: 12/30/22

#' Title Compile TMB Model
#'
#' @param wd Working directory of TMB template
#' @param cpp .cpp file name of TMB template (character)
#'
#' @return
#' @export
#'
#' @examples

compile_tmb <- function(wd, cpp) {
  setwd(wd) # Set working directory
  TMB::compile(cpp) # Compile TMB 
  dyn.load(TMB::dynlib(unlist(strsplit(cpp, split = ".cpp")))) # load TMB template in
} 


#' Title Take additional newton steps with TMB model
#'
#' @param n.newton number of additional newton steps we want to take
#' @param ad_model MakeADFUN model object
#' @param mle_optim Optimized model object
#'
#' @return
#' @export
#'
#' @examples

add_newton <- function(n.newton, ad_model, mle_optim) {
  
  tryCatch(expr = for(i in 1:n.newton) {
    g = as.numeric(ad_model$gr(mle_optim$par))
    h = optimHess(mle_optim$par, fn = ad_model$fn, gr = ad_model$gr)
    mle_optim$par = mle_optim$par - solve(h,g)
    mle_optim$objective = ad_model$fn(mle_optim$par)
  }, error = function(e){e})
  
}


#' Title Run Estimation Model
#'
#' @param data list of data inputs
#' @param parameters list of parameter random starting values
#' @param map map to fix parameters
#' @param n.newton number of additional newton steps to take
#' @param iter.max number of iterations for nlminb to run (default = 1e5)
#' @param eval.max number of nlminb evaluations (default = 1e5)
#' @param silent whether or not to print stuff out
#' @param getsdrep whether to return a sdreport
#'
#' @return
#' @export
#'
#' @examples

run_EM <- function(data, parameters, map, n.newton, random = NULL,
                   iter.max = 1e6, eval.max = 1e6, 
                   silent = TRUE, getsdrep = TRUE) {
  
  # Make AD Function here
  model_fxn <- TMB::MakeADFun(data, parameters, map, random = random,
                              DLL="EM", silent = silent, 
                              checkParameterOrder = TRUE, tracepar = TRUE)
  
  # Optimize model here w/ nlminb
  mle_optim <- stats::nlminb(model_fxn$par, model_fxn$fn, model_fxn$gr, 
                             control = list(iter.max = iter.max, eval.max = eval.max))
  
  # Take additional newton steps
  add_newton(n.newton = n.newton, ad_model = model_fxn, mle_optim = mle_optim)
  
  if(getsdrep == TRUE) {
    # Get report with mle optimized parameters
    model_fxn$rep <- model_fxn$report(model_fxn$env$last.par.best) # Need to pass both fixed and random effects!!!
    # Get sd report here from TMB
    sd_rep <- TMB::sdreport(model_fxn)
  } # if get sdrep = TRUE
  
  return(list(model_fxn = model_fxn, mle_optim = mle_optim, sd_rep = sd_rep))
  
}


#' Title Check TMB Model Convergence (PD Hessian and parameter gradients)
#'
#' @param mle_optim MLE optimized object by nlminb or optim
#' @param sd_rep sd report object from TMB 
#' @param min_grad Minimum gradient we want to determine convergence
#' @param mod_rep Model object with report stored in the list
#'
#' @return List of objects (dim 1 = Convergence status; dim 2 = maximum gradient of the model,
#' dim 3 = parameter with the maximum gradient)
#' @export 
#'
#' @examples

check_model_convergence <- function(mle_optim, sd_rep, mod_rep, min_grad = 0.001) {
  
  # Maximum gradient of the model
  max_grad_val <- max(abs(sd_rep$gradient.fixed))
  
  # Parameter with maximum gradient
  max_grad_par <- names(sd_rep$par.fixed)[which.max(abs(sd_rep$gradient.fixed))]
  
  if(sd_rep$pdHess == TRUE & 
     max_grad_val < min_grad &
     !is.nan(sum(model$sd_rep$sd )) &
     !is.nan(mod_rep$rep$jnLL)) convergence = "Converged"
  else convergence = "Not Converged"
  
  return(list(Convergence = convergence, 
              Max_Grad = max_grad_val, 
              Max_Grad_Par = max_grad_par,
              jnLL = mod_rep$rep$jnLL))
}

#' Title Extract MLE estimates and standard errors from sd_rep object for derived variables
#'
#' @param sd_rep sd_rep object from TMB
#' @param par parameter name we want to extract
#'
#' @return Dataframe with MLE values and error estimates, with 95% confidence intervals
#' @export
#'
#' @examples

extract_ADREP_vals <- function(sd_rep, par) {
  
  # Get sd rep MLE estimate
  mle_val <- sd_rep$value[str_detect(names(sd_rep$value), par)]
  
  # Get sd rep predicted standard errors via the delta method
  mle_se <- sd_rep$sd[str_detect(names(sd_rep$value), par)]
  
  # Put this into a dataframe
  mle_df <- data.frame(mle_val = mle_val, mle_se = mle_se,
                       lwr_95 = mle_val - (1.96 * mle_se ),
                       upr_95 = mle_val + (1.96 * mle_se ))
  
  return(mle_df)
  
}

#' Title Extract MLE estimates and standard errors from sd_rep object for parameter estimates
#'
#' @param sd_rep sd_rep object from TMB
#' @param par parameter name we want to extract
#' @param trans parameter space a given parameter is estimated in ("none", "log", "logit")
#' @param logit_bounds if estimated in logit space, what are the specified bounds
#' @return Dataframe with MLE values and error estimates, with 95% confidence intervals
#' @export
#'
#' @examples

extract_parameter_vals <- function(sd_rep, par, trans = NA, logit_bounds = NA) {
  
  # Get parameter MLE estimates
  par_vals <- sd_rep$par.fixed[names(sd_rep$par.fixed) == par]
  
  # Get parameter variance estimates
  par_var_vals <- diag(sd_rep$cov.fixed)[names(diag(sd_rep$cov.fixed)) == par]
  
  # Put these values into a dataframe
  if(trans == "log") { # log transformation
    mle_df <- data.frame(trans_mle_val = par_vals, mle_val = exp(par_vals), 
                         mle_var = par_var_vals, mle_sd = sqrt(par_var_vals), 
                         lwr_95 = exp(par_vals - (1.96 * sqrt(par_var_vals))),
                         upr_95 = exp(par_vals + (1.96 * sqrt(par_var_vals))))
  } 
  
  if(trans == "logit"){ # log transformation
    
    # Inverse logit parameter values
    inv_par_vals <- logit_bounds[1] + (logit_bounds[2] - logit_bounds[1])/
      (1 + exp(-par_vals))
    
    # Get confidence intervals in logit space
    logit_q_lwr_95 <- par_vals - qnorm(0.975)*sqrt(par_var_vals)
    logit_q_upr_95 = par_vals + qnorm(0.975)*sqrt(par_var_vals)
    
    # Transform into normal space
    lwr_95 <- logit_bounds[1] + (logit_bounds[2] - logit_bounds[1])/  (1 + exp(-logit_q_lwr_95))
    upr_95 <- logit_bounds[1] + (logit_bounds[2] - logit_bounds[1])/  (1 + exp(-logit_q_upr_95))
    
    # Now construct dataframe
    mle_df <- data.frame(trans_mle_val = par_vals, mle_val = inv_par_vals, 
                         mle_var = par_var_vals, mle_sd = sqrt(par_var_vals), 
                         lwr_95 = lwr_95, upr_95 = upr_95)
  } 
  
  if(trans == "none") { # No transofmration needed
    mle_df <- data.frame(trans_mle_val = par_vals, mle_var = par_var_vals, mle_sd = sqrt(par_var_vals), 
                         lwr_95 = par_vals - (1.96 * sqrt(par_var_vals)),
                         upr_95 = par_vals + (1.96 * sqrt(par_var_vals)))
  }
  
  return(mle_df)
}

#' Title Extract true and predicted mean age values for a given fleet and composition type
#'
#' @param mod_rep Model object created by ADFUN
#' @param comp_name Composition name within the model_rep$rep object
#' @param bins Number of age or length bins
#' @param comp_start_yr When composition data collection begins 
#' @param sim Simulation number
#' @param n_fish_true_fleets Number of true fishery fleets within the OM - default = NULL for survey comps
#'
#' @return
#' @export
#'
#' @examples
#' 
extract_mean_age_vals <- function(mod_rep, comp_name, bins, comp_start_yr,
                                  sim, n_fish_true_fleets = NULL) {
  
  # Extract predicted age comps
  pred_age_comps <- mod_rep$rep[names(mod_rep$rep) == comp_name]
  
  # Munge to get mean age for a given year
  molten_pred_age <- reshape2::melt(pred_age_comps)
  colnames(molten_pred_age) <- c("year", "age", "fleet", "sex", "value", "type")
  
  # Summarize to get mean age by year, fleet, and sex
  pred_mean_ages <- molten_pred_age %>% 
    dplyr::group_by(year, fleet, sex) %>% 
    dplyr::summarize(pred_mean_age = sum(value * bins)) %>% 
    mutate(sim = sim)

  # Create global objects for use in loops
  n_years <- length(Fish_Start_yr[1]:(Fish_Start_yr[1] + max(pred_mean_ages$year) - 1))
  
  # Get number of modelled fleets
  n_mod_fleets <- length(unique(pred_mean_ages$fleet))
  
  # Get number of sexes
  n_sexes <- length(unique(pred_mean_ages$sex))
  
  # Create a matrix to store values in
  true_ages <- array(NA, dim = c(n_years, n_mod_fleets, n_sexes))
  
  # Get true mean ages if sampling without error
  # If the composition type = survey comps
  if(comp_name == "pred_srv_age_comps") { 
    
    for(i in 1:n_years) {
      for(f in 1:n_mod_fleets) {
        for(s in 1:n_sexes) {
          # Get true age proportions
          prop_ages <- N_at_age[(Fish_Start_yr[1] - 1 + i),,s,sim] * Surv_selex_at_age[(Fish_Start_yr[1] - 1 + i),,f,s,sim]
          # Get true mean age
          true_ages[i, f, s] <- sum((prop_ages / sum(prop_ages)) * bins)
        } # s loop
      } # f loop
    } # i loop
    
  } else{   # If the composition type = fishery comps
    
    if(is.null(n_fish_true_fleets)) stop("Please specify the true number of fishery fleets")
    
    # Single Fleet OM and EM
    if(n_mod_fleets == 1 & n_fish_true_fleets == 1) { 
      for(i in 1:n_years) {
        for(s in 1:n_sexes) {
          # Get true age proportions - sum across fleets
          prop_ages <- Catch_at_age[(Fish_Start_yr[1] - 1 + i),,1,s,sim]
          # Get true mean age
          true_ages[i, 1, s] <- sum((prop_ages / sum(prop_ages)) * bins)
        } # s loop
      } # i loop
    } # end if this is a single fleet
    
    # Single Fleet in EM, but multi fleet in OM
    if(n_mod_fleets == 1 & n_fish_true_fleets > 1) { 
      for(i in 1:n_years) {
        for(s in 1:n_sexes) {
          # Get true age proportions - sum across fleets for catcha t age
          prop_ages <- rowSums(Catch_at_age[(Fish_Start_yr[1]  - 1 + i),,,s,sim])
          # Get true mean age
          true_ages[i, 1, s] <- sum((prop_ages / sum(prop_ages)) * bins)
        } # s loop
      } # i loop
    } # end if this is modeled as a single fleet when there are > 1 fleet
    
    # Multi fleet in OM and EM
    if(n_mod_fleets > 1 & n_fish_true_fleets > 1) {
      
      if(n_mod_fleets != n_fish_true_fleets) {
        stop("Number of modelled fleets != number of true fleets specified")
      } # warning stop
      
      # if there are separate fleets in the true OM and the EM
      for(i in 1:n_years) {
        for(f in 1:n_fish_true_fleets) {
          for(s in 1:n_sexes) {
            # Get true age proportions - sum across fleets
            prop_ages <- Catch_at_age[(Fish_Start_yr[f]  - 1 + i),,f,s,sim] 
            # Get true mean age
            true_ages[i, f, s] <- sum((prop_ages / sum(prop_ages)) * bins)
          } # s loop
        } # f loop 
      } # i loop
      
    } # if statement
    
  } # else statement
  
  # Rename and bind to the predicted dataframe
  molten_true_ages <- reshape2::melt(true_ages) 
  colnames(molten_true_ages) <- c("year", "fleet", "sex", 'true_mean_ages')
  
  # Left join to complete dataframe
  mean_ages_df <- pred_mean_ages %>% 
    dplyr::left_join(molten_true_ages, by  = c("year", "fleet", "sex"))
  
  return(mean_ages_df)
  
} # end function


#' Title Return RE percentiles
#'
#' @param df dataframe we want to get percentiles from
#' @param est_val_col column number of estimated mle values
#' @param true_val_col column number of true values
#' @param par_name character name that we want to input into the returned dataframe
#' @param group_vars group_by variable names as strings
#'
#' @return dataframe of a range of different percentiles
#' @export
#'
#' @examples
get_RE_precentiles <- function(df, est_val_col = 1, true_val_col = 5, par_name = NULL,
                               group_vars) {
  
  # Get relative error based on indexed columns
  RE <- (df[,est_val_col] - df[,true_val_col]) / df[,true_val_col]
  df <- cbind(df, RE) # Cbind RE
  names(df)[ncol(df)] <- "RE" # Rename variable
  
  df <- df %>% 
    group_by(!!! syms(group_vars)) %>% 
    summarize(median = median(RE), 
              lwr_100 = quantile(RE, 0),
              upr_100 = quantile(RE, 1),
              lwr_95 = quantile(RE, 0.025),
              upr_95 = quantile(RE, 0.975),
              lwr_80 = quantile(RE, 0.1),
              upr_80 = quantile(RE, 0.9),
              lwr_75 = quantile(RE, 0.125),
              upr_75 = quantile(RE, 0.875)) %>% 
    mutate(par_name = par_name) 
  
  if(str_detect(par_name, "Age")) { # if this is an age variable
    
    # Make par_name different such that we can facet wrap later on
    df$par_name <- paste(df$par_name, "Fleet", df$fleet, "Sex", df$sex)
    
    # Drop grouped columns (2nd and 3rd column - fleet and sex)
    df <- df[,-c(2:3)]
  }
  
  return(df)
}


#' Title To extract estiamted and true values from model runs
#'
#' @param sd_rep sd report from TMB
#' @param model_fxn model_fxn from tmb
#' @param sim simulation number
#' @param conv convergence statistic
#' @param n_years number of years in model
#' @param ages number of ages in model
#' @param F_x x% of SPR
#' @param fish_selex_opt Selex fishery options in EM
#' @param ntrue_fish_fleets Number of true fishery fleets 
#' @param nmod_fish_fleets Number of modelled fishery fleets
#'
#' @return
#' @export
#'
#' @examples

get_quants <- function(sd_rep, 
                       model_fxn,
                       sim, 
                       conv, 
                       n_years, 
                       ages,
                       F_x,
                       fish_selex_opt,
                       ntrue_fish_fleets,
                       nmod_fish_fleets) {
  
  # Define parameters we want to get
  par_name <- c("ln_q_srv", 
                "ln_RecPars", 
                "ln_srv_selpars", 
                "ln_fish_selpars",
                "fixed_sel_re_fish")
  
  # True quantities
  t <- c(mean(q_Surv), # Catchability
         r0, # Virgin Recruitment
         NA,  # Srv selex
         NA,  # Fish selex
         NA) # Fish selex variance pars
  
  # Empty dataframe
  par_all <- data.frame()
  # Define years
  years <- 1:n_years
  
  # Loop through to extract values
  for(par in 1:length(par_name)) {
    par_df <- extract_parameter_vals(sd_rep = sd_rep, par = par_name[par], trans = "log") %>%
              dplyr::mutate(type = par_name[par], t = t[par])
    
    # Organizing parameter names for fishery selectivity
    if(par_name[par] == "ln_fish_selpars") {
      n_sexes <- dim(model_fxn$rep$F_Slx)[4] # number of sexes
      # 1 fleet options
      if(sum(fish_selex_opt == "logistic") == 1) {
        if(n_sexes == 1) {
          paste_names <- c("F1_a50", "F1_k")
        } # sex 1
        if(n_sexes == 2) {
          paste_names <- c("F1_a50_m", "F1_a50_f", 
                           "F1_k_m", "F1_k_f")
        }
      } # end if only logistic
      if(sum(fish_selex_opt == "gamma") == 1) {
        if(n_sexes == 1) {
          paste_names <- c("F1_delta", "F1_amax")
        } # sex 1
        if(n_sexes == 2) {
          paste_names <- c("F1_delta_m", "F1_delta_f", 
                           "F1_amax_m", "F1_amax_f")
        }
      } # end gamma
      
      if(sum(fish_selex_opt == "exp_logistic") == 1) {
        if(n_sexes == 1) {
          paste_names <- c("F1_gamma", "F1_alpha", "F1_beta")
        } # sex 1
        if(n_sexes == 2) {
          paste_names <- c("F1_gamma_m", "F1_gamma_f",
                           "F1_alpha_m", "F1_alpha_f",
                           "F1_beta_m", "F1_beta_f")
        }
      } # end exponential logistic
      # 2 fleet options
      if(sum(fish_selex_opt == c("logistic", "logistic")) == 2) {
        if(n_sexes == 1) {
          paste_names <- c("F1_a50", "F2_a50", "F1_k", 'F2_k')
        } # sex 1
        if(n_sexes == 2) {
          paste_names <- c("F1_a50_m", "F2_a50_m", 
                           "F1_a50_f", "F2_a50_f",
                           "F1_k_m", "F2_k_m", 
                           "F1_k_f", "F2_k_f")
        } # 2 sexes
      } # end if logistic and logistic combination
      if(sum(fish_selex_opt == c("logistic", "gamma")) == 2) {
        if(n_sexes == 1) {
          paste_names <- c("F1_a50", "F2_delta",  
                           "F1_k", "F2_amax")
        } # 1 sex
        if(n_sexes == 2) {
          paste_names <- c("F1_a50_m", "F2_delta_m", 
                           "F1_a50_f", "F2_delta_f",
                           "F1_k_m", "F2_amax_m", 
                           "F1_k_f", "F2_amax_f")
        } # 2 sexes
      } # end if logistic and gamma combination
      
      if(sum(fish_selex_opt == c("gamma", "logistic")) == 2) { # if gamma and logistic fleet 1 and 2 
        if(n_sexes == 1) {
          paste_names <- c("F1_delta", "F2_a50",  
                           "F1_amax", "F2_k")
        } # 1 sex
        if(n_sexes == 2) {
          paste_names <- c("F1_delta_m", "F2_a50_m", 
                           "F1_delta_f", "F2_a50_f", 
                           "F1_amax_m", "F2_k_m", 
                           "F1_amax_f", "F2_k_f")
        } # 2 sexes
      } # end if gamma and logistic combination
      
      par_df$type <- paste_names
    } # end if fishery selex parameters
    
    if(par_name[par] == "ln_srv_selpars") {
      n_sexes <- dim(model_fxn$rep$S_Slx)[4] # number of sexes
      if(n_sexes == 1) {
        paste_names <- c("S1_a50", "S1_k")
      }  # if nsex = 1
      if(n_sexes == 2) {
        paste_names <- c("S1_a50_m", "S1_a50_f",
                         "S1_k_m", "S1_k_f")
      } # if nsex = 2
      par_df$type <- paste_names
    } # end if survey slx pars
    
    par_all <- rbind(par_all, par_df)
  } # end par loop
  
  # Get SSB0
  ssb0_df <- extract_ADREP_vals(sd_rep = sd_rep, par = "ssb0") %>% 
    dplyr::mutate(type = "SSB0", trans_mle_val = NA,
           mle_var = NA, t = ssb0) %>%  dplyr::rename(mle_sd = mle_se)
  par_all <- rbind(par_all, ssb0_df)
  
  # Get Fx here
  # Pre-Processing first
  # Get Fs here
  if(nmod_fish_fleets > 1) {
    # Vector of fleets (define to multiply)
    fleet_num <- seq(1, nmod_fish_fleets)
    # Extract Fs here
    est_TermF <- exp(sd_rep$par.fixed[names(model$sd_rep$par.fixed) == "ln_Fy"])[c(length(years) * fleet_num)]
  } else{ # only 1 fleet
    # Extract Fs here
    est_TermF <- exp(sd_rep$par.fixed[names(model$sd_rep$par.fixed) == "ln_Fy"])[c(length(years))]
  } # else statement

  est_Selex <- model_fxn$rep$F_Slx[years[length(years)],,,1] # selectivity only for females
  
  # Get estaimted Fx% value
  Fx_val <- get_Fx_refpt(ages = ages,
                          MortAA = Mort_at_age[length(years),,sim], 
                          SelexAA = t(est_Selex), 
                          MatAA = mat_at_age[length(years),,1,sim], # females
                          WAA = wt_at_age[length(years),,1,sim],  # females
                          Terminal_F = est_TermF, 
                          F_x = F_x)
  
  # Get true f40
  true_f40 <- get_Fx_refpt(ages = ages,
                           MortAA = Mort_at_age[length(years),,sim], 
                           SelexAA = t(Fish_selex_at_age[length(years),,,1,sim]),  # females selex
                           MatAA = mat_at_age[length(years),,1,sim],  # females
                           WAA = wt_at_age[length(years),,1,sim],  # females
                           Terminal_F = fish_mort[length(years),,sim],
                           F_x = F_x)

  # Create empty df for Fx% to rbind to par all dataframe
  F_df <- data.frame(matrix(ncol = length(names(par_all))))
  colnames(F_df) <- names(par_all) # replace colnames in empty df
  
  # Now, input into empty dataframe
  F_df <- F_df %>% dplyr::mutate(mle_val = Fx_val, 
                                 type = paste("F", F_x, sep = "_"),
                                 t = true_f40)
  
  # Now, bind all quantities of interest into one df
  par_all <- rbind(F_df, par_all)
  
  # Input idenitfier and converence stats
  par_all <- par_all %>% dplyr::mutate(conv = conv, sim = sim)
  row.names(par_all) <- NULL # Get rid of pesky row names
  
  # Extract time series quantities
  # Recruitment
  rec_df <- extract_ADREP_vals(sd_rep = sd_rep, par = "Total_Rec") %>% 
    dplyr::mutate(t = rec_total[years,sim], sim = sim, conv = conv, year = years,
           type = "Total Recruitment")

  # SSB
  ssb_df <- extract_ADREP_vals(sd_rep = sd_rep, par = "SSB") %>% 
    dplyr::mutate(t = SSB[years, sim], sim = sim, conv = conv, year = years,
           type = "Spawning Stock Biomass")
  
  # Get total fishing mortality
  f_df <- extract_ADREP_vals(sd_rep = sd_rep, par = "Total_Fy") %>% 
    dplyr::mutate(t = rowSums(matrix(fish_mort[years,,sim], ncol = ntrue_fish_fleets)), 
           sim = sim, conv = conv, year = years, type = "Total Fishing Mortality")
  
  # Get depletion
  depletion_df <- extract_ADREP_vals(sd_rep = sd_rep, par = "Depletion") %>%
    dplyr::mutate(t = (SSB[years,sim]/ssb0), sim = sim, conv = conv, year = years,
           type = "Depletion")
  
  ts_all <- rbind(rec_df, ssb_df, f_df, depletion_df)
  
  return(list(TS_df = ts_all, Par_df = par_all))

} # end function

#' Title Get Results from simulations (time series and parameters)
#'
#' @param om_scenario_path Path to Opearting Model Scenarios
#'
#' @return Returns a list of dataframes
#' @export
#'
#' @examples
get_results <- function(om_scenario_path) {
  
  require(here)
  require(tidyverse)
  
  # Get OM Scenarios from OM path
  om_scenarios <- list.files(om_scenario_path)
  # Set up empty dataframes
  param_all <- data.frame()
  ts_all <- data.frame()
  
  # Loop thorugh to extract metrics
  for(i in 1:length(om_scenarios)) {
    
    em_scenarios <- list.files(here(om_scenario_path, om_scenarios[i]))
    # Remove .Rdata and .pdf from em_scenarios
    em_scenarios <- em_scenarios[str_detect(em_scenarios, ".RData") == FALSE]
    em_scenarios <- em_scenarios[str_detect(em_scenarios, ".pdf") == FALSE]
    
    for(j in 1:length(em_scenarios)) {
      # Read in parameters and time series csvs
      params <- read.csv(here(om_scenario_path, om_scenarios[i], 
                              em_scenarios[j], "Param_Results.csv"))
      ts <- read.csv(here(om_scenario_path, om_scenarios[i], 
                          em_scenarios[j], "TimeSeries_Results.csv"))
      
      param_all <- rbind(param_all, params)
      ts_all <- rbind(ts_all, ts)
      
    } # j loop
  } # i loop
  
  return(list(Parameter_Sum = param_all,
              TimeSeries_Sum = ts_all))
} # end function
