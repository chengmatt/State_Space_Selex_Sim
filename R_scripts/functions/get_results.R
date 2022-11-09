# Purpose: To collate results from simulation study - meant to be used
# after OM-EM models have been run! 
# Within this script, we have the following functions:
# 1) get_sdrep_values (to get sdrep values), 2) get_results (to get results from models - uses
# the get_sdrep_values function)
# Creator: Matthew LH. Cheng
# Date 11/4/22


# Get sdrep values --------------------------------------------------------

#' @param parameter Parameter name from sdreport
get_sdrep_values <- function(sd_rep = sd_rep, parameter) {
  
  idx <- which(rownames(sd_rep) == parameter) # find out where in the sdrep the 
  # parameter is
  
  if(parameter == "log_FAA_tot") { # Fishing mortality by age
    
    # The code here is adapted from WHAM - from plot.SSB.F.trend()
    
    # Stick these values into a matrix (nrow = years, length = ages)
    log_faa <- matrix(sd_rep[idx,1], length(Fish_Start_yr:(n_years-1)), length(ages))
    
    # Get CIs for the above
    faa_ci <- matrix(sd_rep[idx,2], length(Fish_Start_yr:(n_years-1)), length(ages))
    
    # Get max F - fully selected F 
    full_selex_F_index <- apply(log_faa,1, function(x) max(which(x == max(x)))) 
    
    # Now, get fully selected F values and associated SEs
    # Create empty vectors to store vlaues in
    value <- vector()
    se <- vector()
    
    # Loop through this because we might get years where the fully selected F differs by yrs
    for(y in 1:length(full_selex_F_index)) {
      
      # Get year index for fully selected F
      f_index <- full_selex_F_index[y]
      
      # search for the yth index and put into our vec
      value[y] <- log_faa[y,f_index]
      se[y] <- faa_ci[y,f_index]
      
    } # end y loop

  } else{ 
    value <- sd_rep[idx,1] # Get estimated values from sdrep
    se <- sd_rep[idx,2] # get standard errors
  } # else - not log_FAA_tot
  
  return(list(value, se)) # return values
  
} # function to get values from sdrep


# Get general results from sdrep -----------------------------------------------------

#' @param model Model object estimated from WHAM

get_results <- function(model) {
  
  require(tidyverse)
  
  # Get report file here from model
  report <- model$report()
  
  # Get sdreport file here from model
  sd_rep <- summary(model$sdrep)
  
# SSB, F, and q -----------------------------------------------------------

  # Get SSB values from sdrep
  ssb_values <- get_sdrep_values(sd_rep = sd_rep, parameter = "log_SSB")
  
  # Stick these into a dataframe and calculate CIs
  SSB_rep <- data.frame(Year = Fish_Start_yr:(n_years-1), log_SSB = ssb_values[[1]], SE = ssb_values[[2]]) %>% 
    mutate(SSB = exp(log_SSB),
           lwr = exp(log_SSB - (1.96*SE)),
           upr = exp(log_SSB + (1.96*SE)))
  
  # Get fishing mortality rates
  f_values <- get_sdrep_values(sd_rep = sd_rep, parameter = "log_FAA_tot")
  F_rep <- data.frame(Year = Fish_Start_yr:(n_years-1), log_F = f_values[[1]], SE = f_values[[2]]) %>% 
    mutate(F_val = exp(log_F),
           lwr = exp(log_F - (1.96*SE)),
           upr = exp(log_F + (1.96*SE)))

# Selectivity -------------------------------------------------------------
  
  # Get fishery selectivity
  fish_sel <- melt(report$selAA[[1]]) %>% 
    rename(Year = Var1, Age = Var2, Selex = value) %>% 
    mutate(Type = "Fishery")
  
  # Get survey selectivity
  surv_sel <- melt(report$selAA[[2]]) %>% 
    rename(Year = Var1, Age = Var2, Selex = value) %>% 
    mutate(Type = "Survey")
  
  selex_rep <- rbind(fish_sel, surv_sel) # bind these together!
  

# Get catchability --------------------------------------------------------
  if("q_re" %in% model$input$random) { # taken from WHAM plot_q function
    se <- as.list(model$sdrep, "Std. Error", report=TRUE)$logit_q_mat
  } else {
    se <- t(matrix(as.list(model$sdrep, "Std. Error")$logit_q, nrow = NCOL(model$rep$logit_q_mat), 
                     ncol = NROW(model$rep$logit_q_mat)))  
  } # if else statement as to whether this is a random effect or not
  
  # Get bounds of logit
  logit_q_lo <- model$rep$logit_q_mat - qnorm(0.975)*se
  logit_q_hi <- model$rep$logit_q_mat + qnorm(0.975)*se
  
  # Now, get inverse logit q
  q <- t(model$input$data$q_lower + (model$input$data$q_upper - model$input$data$q_lower)/
                                                      (1+exp(-t(model$rep$logit_q_mat))))
  q_lo <- t(model$input$data$q_lower + (model$input$data$q_upper - model$input$data$q_lower)/
                                                                    (1+exp(-t(logit_q_lo))))
  q_hi <- t(model$input$data$q_lower + (model$input$data$q_upper - model$input$data$q_lower)/
                                                                    (1+exp(-t(logit_q_hi))))
  
  # Now, put these into a dataframe
  # fishery
  q_rep_Fish <- data.frame(Year = Fish_Start_yr:(n_years-1), q_Fish = q[,1], 
                           q_Fish_lwr = q_lo[,1], q_Fish_upr = q_hi[,1])
  # survey
  q_rep_Surv <- data.frame(Fish_Start_yr:(n_years-1), q_Surv = q[,1], 
                           q_Surv_lwr = q_lo[,1], q_Surv_upr = q_hi[,1])
  
  
  return(list(SSB = SSB_rep, F = F_rep, Selex = selex_rep, q_Fish = q_rep_Fish,
              q_Surv = q_rep_Surv))
}



# OM-EM comparisons -------------------------------------------------------

#' @param n_sims Number of total simulations we ran
#' @param EM_variable Variable name in the get_results object for the EM
#' @param OM_df Variable name for the OM

# Returns a df of unaggregated resutls paired with the OM EM, and 95%, 50% CIs w/ relative error
# Right now, this function only does F and SSB and q

om_em_results <- function(n_sims, EM_variable, OM_df) {
  
  df <- data.frame() # create empty dataframe to hold values in
  
  # Get SSB from EMs and compare to OMs
  for(sim in 1:n_sims) {
    
    # Get ssb here from EM
    em_df <- results_list[[sim]][[EM_variable]]
    
    # Next, get ssb from the corresponding OM simulation
    om_df <- data.frame(OM_df[em_df$Year,sim]) %>% 
      rename(Truth = OM_df.em_df.Year..sim.) %>% 
      mutate(Sim = sim)
    
    # Now, cbind these together
    om_em_df <- cbind(em_df, om_df)
    
    # And now rbind all of these together
    df <- rbind(om_em_df, df) # Unaggregated data
    
  } # end sim loop
  
  # Return aggregated data w/ relative error, 95% and 50% CIs
  if(EM_variable == "SSB") {
    df_agg <- df %>% 
      drop_na() %>% 
      mutate(RE = (SSB - Truth)/Truth) %>% 
      group_by(Year) %>% 
      summarize(mean = mean(RE), sd = sd(RE),  n = n(),
                lwr_95 = mean - (1.96 * (sd / sqrt(n))),
                up_95 = mean + (1.96 * (sd / sqrt(n))),
                lwr_50 = mean - (0.674 * (sd / sqrt(n))),
                up_50 = mean + (0.674 * (sd / sqrt(n))))
  } # end if variable == SSB
  
  if(EM_variable == "F") {
    
    df_agg <- df %>% 
      drop_na() %>% 
      mutate(RE = (F_val - Truth)/Truth) %>% 
      group_by(Year) %>% 
      summarize(mean = mean(RE), sd = sd(RE),  n = n(),
                lwr_95 = mean - (1.96 * (sd / sqrt(n))),
                up_95 = mean + (1.96 * (sd / sqrt(n))),
                lwr_50 = mean - (0.674 * (sd / sqrt(n))),
                up_50 = mean + (0.674 * (sd / sqrt(n))))
    
  } # end if variable == fishing mortality
  
  if(str_detect(EM_variable, "q")) {
    
    df_agg <- df %>% 
      drop_na() # drops nas and nans
    
    # Next, create our relative error metric - we need to use indexing
    # because we the q's differ by survey and fishery
    df_agg$RE <- (df_agg[,2] - df_agg$Truth) / df_agg$Truth
    
    # Finally, summarize this!
    df_agg <- df_agg %>% 
      group_by(Year) %>% 
      summarize(mean = mean(RE), sd = sd(RE),  n = n(),
                lwr_95 = mean - (1.96 * (sd / sqrt(n))),
                up_95 = mean + (1.96 * (sd / sqrt(n))),
                lwr_50 = mean - (0.674 * (sd / sqrt(n))),
                up_50 = mean + (0.674 * (sd / sqrt(n))))
    
  } # if this is a catchability variable
  
  return(list(df, df_agg))
  
} # end function
