# Purpose: Run EMs that are 1 Fleet and are all random walks
# Creator: Matthew LH. Cheng
# Date 6/12/23

# Set up ------------------------------------------------------------------
rm(list=ls()) # remove objects prior to running

library(here)
library(tidyverse)
library(TMB)
library(doParallel)
library(parallel)

# Load in all functions into the environment
fxn_path <- here("R_scripts", "functions")
# Load in all functions from the functions folder
files <- list.files(fxn_path)
for(i in 1:length(files)) source(here(fxn_path, files[i]))

compile_tmb(wd = here("src"), cpp = "EM.cpp")

# Read in OM and EM Scenarios
om_scenarios <- readxl::read_excel(here('input', "OM_EM_Scenarios_v3.xlsx"), sheet = "OM") 
em_scenarios <- readxl::read_excel(here('input', "OM_EM_Scenarios_v3.xlsx"), sheet = "EM_1Fl_RW") %>% 
  filter(str_detect(EM_Scenario, "SP"))

# Read in spreadsheet for life history parameters
lh_path <- here("input", "Sablefish_Inputs.xlsx")

# Get number of OM and EM scenarios
n_OM_scen <- length(om_scenarios$OM_Scenarios)
n_EM_scen <- length(em_scenarios$EM_Scenario)

for(n_om in 1:n_OM_scen) {
  
# Load OM -----------------------------------------------------------------
  
  # Load in OM data
  om_path <- here("output", "OM_Scenarios", om_scenarios$OM_Scenarios[n_om])
  load(here(om_path, paste(om_scenarios$OM_Scenarios[n_om],".RData",sep = "")))
  list2env(oms,globalenv()) # output into global environment
  
  for(n_em in 1:n_EM_scen) {
    
    rm(sim_models)
    
    # register cores here
    ncores <- detectCores() 
    cl <- makeCluster(ncores - 2)
    registerDoParallel(cl)
    
# Pre-process EM ----------------------------------------------------------

    # Pre-processing to specify EM here
    n_fleets <- em_scenarios$n_fleets[n_em] # Number of fleets to model
    
    # Specify year options - differs if it is Fast vs Slow OM Scenario
    years_opt <- as.numeric(unlist(strsplit(em_scenarios$n_years[n_em], ",")) )
    if(str_detect(om_scenarios$OM_Scenarios[n_om], "Fast")) years <- 1:years_opt[1]
    if(str_detect(om_scenarios$OM_Scenarios[n_om], "Slow")) years <- 1:years_opt[2]
    
    # Fishery selectivity options
    fish_selex_opt <- unlist(strsplit(em_scenarios$Selex[n_em], ",")) # get fishery selectivity options
    # munge into matrix format here
    fish_selex_opt <- matrix(c(rep(fish_selex_opt, length(years))),  nrow = length(years), ncol = n_fleets)
    
    # Time-varying fishery selectivity options
    time_selex <- em_scenarios$Time_Selex[n_em] # Get time-varying selectivity options
    time_selex_npars <- em_scenarios$Time_Selex_Npars[n_em] # Number of time-varying selectivity parameters 
    random_fish_sel <- em_scenarios$Random_Effects[n_em] # wheter or not fishery selex is time-varying via random effects
    fixed_sigma_re_fish <- em_scenarios$Sigma_Fixed[n_em] # get sigma value to fix at
    share_ages_em <- em_scenarios$share_ages[n_em] # get ages to share for semi-parametric selex (if SP)

    # Run Simulations ---------------------------------------------------------
    
    sim_models <- foreach(sim = 1:n_sims, .packages = c("TMB", "here", "tidyverse")) %dopar% {
    
    # for(sim in 1:n_sims) {
    compile_tmb(wd = here("src"), cpp = "EM.cpp")
    dyn.unload(dynlib("EM"))
    dyn.load(dynlib("EM"))
    
    # Choose which parameters to fix
    if(fixed_sigma_re_fish == "NA") fix_pars = c( "ln_SigmaRec", "ln_q_fish", "ln_h", "ln_M")
    
    # If we want to estimate as penalized likelihood
    if(fixed_sigma_re_fish != "NA") {
      fix_pars = c( "ln_SigmaRec", "ln_q_fish", "ln_h", "ln_fixed_sel_re_fish", "ln_M")
      random_fish_sel = NULL # also change RE to null - estimating this as fixed effects
    }
    
    # Prepare inputs for EM
    input <- prepare_EM_input(years = years,
                              n_fleets = n_fleets, 
                              Fish_Start_yr = 1,
                              ages = ages,
                              catch_cv = rep(0.01, n_fleets),
                              
                              # Selectivity Blocks
                              F_Slx_Blocks_Input = matrix(c(rep(0)),
                                                          nrow = length(years), 
                                                          ncol = n_fleets), # fishery blocks
                              S_Slx_Blocks_Input = matrix(c(0), # selectivity blocks
                                                          nrow = length(years), 
                                                          ncol = 1),
                              use_fish_index = FALSE,
                              Fish_Comp_Like_Model = "multinomial",
                              Srv_Comp_Like_Model = "multinomial",
                              rec_model = "BH", 
                              F_Slx_Model_Input = fish_selex_opt,
                              S_Slx_Model_Input = c("logistic"), # logistic
                              time_selex = time_selex,  
                              share_ages = as.numeric(share_ages_em),
                              n_time_selex_pars = as.numeric(time_selex_npars),
                              fix_pars = fix_pars,
                              sim = sim)
    

    # Fix sigma if sigma_fixed is not NA 
    if(fixed_sigma_re_fish != "NA") input$parameters$ln_fixed_sel_re_fish[] <- log(as.numeric(fixed_sigma_re_fish))

    # Run EM model here and get sdrep
    model = list() # to ensure foreach loop doesnt break b/c of not detecting exported obj
    tryCatch(expr = model <- run_EM(data = input$data, parameters = input$parameters, 
                                    map = input$map, n.newton = 3,
                                    random = random_fish_sel, 
                                    silent = T, getsdrep = TRUE), error = function(e){e}) 
    
    # plot(model$model_fxn$rep$F_Slx[1,,,1])
    # model$model_fxn$rep$ln_fish_selpars_re[]
    # model$sd_rep
    # #
    # image(model$model_fxn$rep$F_Slx[,,1,1])
    # for(i in 1:49) {
    #   if(i==1) plot(model$model_fxn$rep$F_Slx[i,,1,1], type = "l", ylim = c(0,3))
    #   else lines(model$model_fxn$rep$F_Slx[i,,1,1], type = "l")
    # }

    # Check model convergence
    convergence_status = list() # to ensure foreach loop doesnt break b/c of not detecting exported obj
    tryCatch(expr = convergence_status <- check_model_convergence(mle_optim = model$mle_optim,
                                                  mod_rep = model$model_fxn,
                                                  sd_rep = model$sd_rep,
                                                  min_grad = 0.001), error = function(e){e}) 
    
    # Get quantities
    quants_df = list() # to ensure foreach loop doesnt break b/c of not detecting exported obj
    tryCatch(expr = quants_df <- get_quants(sd_rep = model$sd_rep,
                                            model_fxn = model$model_fxn, 
                                            sim = sim,  
                                            n_sex = n_sex, 
                                            sex_ratio = 0.5,
                                            conv = convergence_status$Convergence,
                                            n_years = length(years), 
                                            fish_selex_opt = fish_selex_opt,
                                            ntrue_fish_fleets = dim(Fish_Age_Comps)[3], 
                                            nmod_fish_fleets = input$data$n_fleets,
                                            ages = ages, 
                                            F_x = 0.4), error = function(e){e}) 
    
    # Combine objects to save
    list(model, quants_df$Par_df, quants_df$TS_df)
    
  } # end foreach loop
    
    # After we're done running EMs, output objects to save in folder
    
    # Pre-processing here
    model_list <- list()
    params <- data.frame()
    time_series <- data.frame()
    
    for(s in 1:n_sims) {
      model_list[[s]] <- sim_models[[s]][[1]] # save model objects
      params <- rbind(params, sim_models[[s]][[2]]) # save parameters
      time_series <- rbind(time_series, sim_models[[s]][[3]]) # save time series
    } # end s loop
    
    # Now, save our results - create directory to store results first
    em_path_res <- here(om_path, em_scenarios$EM_Scenario[n_em])
    dir.create(em_path_res)
    save(model_list, file = here(em_path_res, paste(em_scenarios$EM_Scenario[n_em], ".RData", sep = ""))) # save models
    
    # Differentiate OM and EMs
    params <- params %>% dplyr::mutate(OM_Scenario = om_scenarios$OM_Scenarios[n_om],
                                       EM_Scenario = em_scenarios$EM_Scenario[n_em])
    time_series <- time_series %>% dplyr::mutate(OM_Scenario = om_scenarios$OM_Scenarios[n_om],
                                                 EM_Scenario = em_scenarios$EM_Scenario[n_em])
    
    # Also, get AIC results
    AIC_df <- data.frame()
    # Get convergence summary here
    Convergence = params %>% group_by(sim) %>% summarize(conv = unique(conv))
    
    for(s in 1:n_sims) {
      # To deal with empty lists - foreach loop issue with objects not exporting
      if(length(model_list[[s]]) != 0) {
        TMB_AIC_val = TMB_AIC(model_list[[s]]$mle_optim)
        conv_val = Convergence$conv[Convergence$sim == s]
      }
      if(length(model_list[[s]]) == 0) {
        TMB_AIC_val = NA
        conv_val = "Not Converged"
      }
      # Get AIC here
      AIC_df_tmp <- data.frame(
        AIC = TMB_AIC_val,
        OM_Scenario = om_scenarios$OM_Scenarios[n_om],
        EM_Scenario = em_scenarios$EM_Scenario[n_em],
        sim = s, conv = conv_val
      )
      AIC_df <- rbind(AIC_df, AIC_df_tmp)
    } # end s sim AIC loop
    
    # Now, output these as csvs
    write.csv(AIC_df, here(em_path_res, "AIC_Results.csv"))   
    write.csv(params, here(em_path_res, "Param_Results.csv"))
    write.csv(time_series, here(em_path_res, "TimeSeries_Results.csv"))
    
    # Progress
    cat(crayon::green("OM", n_om, "out of", n_OM_scen))
    cat(crayon::yellow("EM", n_em, "out of", n_EM_scen))
    stopCluster(cl)
    
  } # end n_em loop
} # end n_om loop


