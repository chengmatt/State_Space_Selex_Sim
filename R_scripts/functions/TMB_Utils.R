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

#' Title Check TMB Model Convergence (PD Hessian and parameter gradients)
#'
#' @param mle_optim MLE optimized object by nlminb or optim
#' @param sd_rep sd report object from TMB 
#' @param min_grad Minimum gradient we want to determine convergence
#'
#' @return List of objects (dim 1 = Convergence status; dim 2 = maximum gradient of the model,
#' dim 3 = parameter with the maximum gradient)
#' @export 
#'
#' @examples

check_model_convergence <- function(mle_optim, sd_rep, min_grad = 0.001) {
  
  # Maximum gradient of the model
  max_grad_val <- max(abs(sd_rep$gradient.fixed))
  
  # Parameter with maximum gradient
  max_grad_par <- names(sd_rep$par.fixed)[which.max(abs(sd_rep$gradient.fixed))]
  
  if(mle_optim$convergence == 0 &
     sd_rep$pdHess == TRUE & 
     max_grad_val < min_grad) convergence = "Converged"
  else convergence = "Not Converged"
  
  return(list(Convergence = convergence, 
              Max_Grad = max_grad_val, 
              Max_Grad_Par = max_grad_par))
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

