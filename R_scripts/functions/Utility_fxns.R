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
                   iter.max = 1e5, eval.max = 1e5, 
                   silent = TRUE, getsdrep = TRUE) {
  
  # Make AD Function here
  model_fxn <- TMB::MakeADFun(data, parameters, map, random = random,
                              DLL="EM", silent = silent, 
                              checkParameterOrder = TRUE)
  
  # Optimize model here w/ nlminb
  mle_optim <- stats::nlminb(model_fxn$par, model_fxn$fn, model_fxn$gr, 
                     control = list(iter.max = iter.max, eval.max = eval.max), 
                     lower = -15, upper = 15)
  
  # Take additional newton steps
  add_newton(n.newton = n.newton, ad_model = model_fxn, mle_optim = mle_optim)
  
  if(getsdrep == TRUE) {
    # Get report with mle optimized parameters
    model_fxn$rep <- model_fxn$report(par = model_fxn$env$last.par.best) # Need to pass both fixed and random effects!!!
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
  
  if(mle_optim$convergence == 0 &
     sd_rep$pdHess == TRUE & 
     max_grad_val < min_grad &
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
#' @param log whether or not parameter is in log space - if TRUE, include backtransformed values
#' @return Dataframe with MLE values and error estimates, with 95% confidence intervals
#' @export
#'
#' @examples

extract_parameter_vals <- function(sd_rep, par, log) {
  
  # Get parameter MLE estimates
  par_vals <- sd_rep$par.fixed[names(sd_rep$par.fixed) == par]
  
  # Get parameter variance estimates
  par_var_vals <- diag(sd_rep$cov.fixed)[names(diag(sd_rep$cov.fixed)) == par]
    
  # Put these values into a dataframe
  if(log == TRUE) { # include exponentiated values and errors if log = TRUE
    mle_df <- data.frame(ln_mle_val = par_vals, mle_val = exp(par_vals), 
                         mle_var = par_var_vals, mle_sd = sqrt(par_var_vals), 
                         lwr_95 = exp(par_vals - (1.96 * sqrt(par_var_vals))),
                         upr_95 = exp(par_vals + (1.96 * sqrt(par_var_vals))))
  } else{ 
    mle_df <- data.frame(mle_val = par_vals, mle_var = par_var_vals, mle_sd = sqrt(par_var_vals), 
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
          true_ages[i, f, s] <- sum((prop_ages / sum(prop_ages)) * 1:30)
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
          true_ages[i, 1, s] <- sum((prop_ages / sum(prop_ages)) * 1:30)
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
          true_ages[i, 1, s] <- sum((prop_ages / sum(prop_ages)) * 1:30)
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
            true_ages[i, f, s] <- sum((prop_ages / sum(prop_ages)) * 1:30)
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

