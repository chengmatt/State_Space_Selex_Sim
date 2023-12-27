  # Purpose: To generate OM datasets for use in EMs
  # Creator: Matthew LH. Cheng
  # Date: 3.24.23
  
  # Set Up ------------------------------------------------------------------
  
  library(here)
  library(tidyverse)
  
  # Load in all functions into the environment
  fxn_path <- here("R_scripts", "functions")
  source(here(fxn_path, "simulate_data.R"))
  source(here(fxn_path, "Utility_fxns.R"))
  source(here(fxn_path, "prepare_EM_input.R"))
  
  # Create folder for OMs
  dir.create(here("output", "OM_Scenarios"))
  
  # Read in OM scenarios
  om_scenarios <- readxl::read_excel(here('input', "OM_EM_Scenarios_v3.xlsx"), sheet = "OM") %>% 
    filter(str_detect(OM_Scenarios, "Ext"))
  
  # Read in spreadsheet for life history parameters
  lh_path <- here("input", "Sablefish_Inputs.xlsx")
  
  # Get number of OM and EM scenarios
  n_OM_scen <- length(om_scenarios$OM_Scenarios)
  
  for(n_om in 1:n_OM_scen) {
  
    # Operating Model Loop + Set Up -------------------------------------------
    
    # Create file directory to save model outputs
    om_path <- here("output", "OM_Scenarios", om_scenarios$OM_Scenarios[n_om])
    dir.create(om_path)
    
    # Pre-Processing
    n_years <- om_scenarios$n_years[n_om] # number of years to run om for
    n_fleets <- om_scenarios$n_fleets[n_om] # number of fleets (fishery)
    yr_chng_end <- om_scenarios$yr_chng_end[n_om] # the year we want the fleet transition to finish
    F_type <- c(om_scenarios$Fl_1_Ftype[n_om], om_scenarios$Fl_2_Ftype[n_om]) # Get Fishing Mortality Pattern
    Start_F <- c(om_scenarios$Fl_1_StartF[n_om], om_scenarios$Fl_2_StartF[n_om])  # Get Starting Fs
    Max_Rel_F_M <- c(om_scenarios$Fl_1_MaxRelM[n_om], om_scenarios$Fl_2_MaxRelM[n_om]) # Get max F relative to M
    Input_Fish_N_Max <- c(om_scenarios$Fl_1_InputN_Max[n_om], om_scenarios$Fl_2_InputN_Max[n_om]) # Get max Input N for comps
    Input_N_Fish_Fixed <- c(om_scenarios$Fl_1_InputN_Fixed[n_om], om_scenarios$Fl_2_InputN_Fixed[n_om]) # Get min Input N for comps
    Input_N_Srv <- om_scenarios$Srv_InputN_Max[n_om] # Survey Input N
    Srv_CV <- om_scenarios$Srv_CV[n_om] # Survey Index CV
    Fish_Selex_Opt <- c(om_scenarios$Fl_1_Selex[n_om], om_scenarios$Fl_2_Selex[n_om]) # Fishery Selectivity Options
    Input_N_Fish_Time = om_scenarios$Input_N_Fish_Time[n_om] # Whether comps are varying or constant
    
    # Selectivity parameters for Fleet 1 Females, followed by Males
    Fish_Fleet1_SelPars <- c(om_scenarios$Fl_1_Slx_Par1_F[n_om], 
                             om_scenarios$Fl_1_Slx_Par2_F[n_om],
                             om_scenarios$Fl_1_Slx_Par1_M[n_om], 
                             om_scenarios$Fl_1_Slx_Par2_M[n_om])
    
    # Selectivity parameters for Fleet 2 Females, followed by Males
    Fish_Fleet2_SelPars <- c(om_scenarios$Fl_2_Slx_Par1_F[n_om], 
                             om_scenarios$Fl_2_Slx_Par2_F[n_om],
                             om_scenarios$Fl_2_Slx_Par1_M[n_om], 
                             om_scenarios$Fl_2_Slx_Par2_M[n_om])  
  
    # Simulate data here
    oms <- simulate_data(fxn_path = fxn_path, 
                  spreadsheet_path = lh_path,
                  n_years = n_years,
                  
                  # Biological controls
                  rec_type = "BH",
                  f_ratio = 0.5, # Assuming 50_50 Sex Ratio
                  m_ratio = 0.5,
                  
                  # Fishery Controls
                  Fish_Start_yr = c(1, 1), 
                  F_type = F_type,
                  Start_F = Start_F, 
                  max_rel_F_M = Max_Rel_F_M, 
                  desc_rel_F_M = c(NULL, NULL), 
                  
                  yr_chng = rep(25, n_fleets), # when the fleet transition begins
                  yr_chng_end = rep(yr_chng_end, n_fleets), # when we want the flee transition to stop
                  fish_likelihood = "multinomial",
                  Input_Fish_N_Max = Input_Fish_N_Max, 
                  fish_CV = c(0.1, 0.1), # not really used
                  Input_N_Fish_Time = Input_N_Fish_Time, 
                  Input_N_Fish_Fixed = Input_N_Fish_Fixed, # min of Neff comps when allowed to vary
                  catch_CV = c(0.0, 0.0), 
                  q_Mean_Fish = c(0.01, 0.01),
                  
                  # Survey Controls
                  Surv_Start_yr = c(1),  
                  srv_likelihood = "multinomial",
                  Input_Srv_N_Max = Input_N_Srv,
                  srv_CV = Srv_CV, 
                  q_Mean_Surv = 0.01, # doesn't change
                  
                  # Selectivity Controls
                  fish_selex = Fish_Selex_Opt, 
                  srv_selex = c("logistic"), 
                  fish_pars = list(Fish_Fleet1 = matrix(data = c(Fish_Fleet1_SelPars), 
                                                        nrow = 2, byrow = TRUE), # Fishery Fleet 1
                                   Fish_Fleet2 = matrix(data = c(Fish_Fleet2_SelPars), 
                                                        nrow = 2, byrow = TRUE)), 
                  srv_pars = list(Srv_Fleet1 = matrix(data = c(2.5, 0.65, 4.5, 0.85), 
                                                      nrow = 2, byrow = TRUE))) # survey fleet 1
  
    # Save as RData file - ifelse for sensitivity tests
      save(oms, file = here(om_path, paste(om_scenarios$OM_Scenarios[n_om], ".RData", sep = "")))
      
      # Plot OM here
      plot_OM(path = here("output", "OM_Scenarios", om_scenarios$OM_Scenarios[n_om]), 
              file_name = "OM_Plots.pdf")
  } # end n_om loop
