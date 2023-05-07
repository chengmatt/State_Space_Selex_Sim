# Purpose: Generate OMs as test cases (quick shift w/ shortened period in which the first fleet is dominant)
# Creator: Matthew LH. Cheng
# Date 4/14/23


# Set up ------------------------------------------------------------------

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
om_scenarios <- readxl::read_excel(here('input', "OM_EM_Scenarios.xlsx"), sheet = "OM_TestCases")

# Read in spreadsheet for life history parameters
lh_path <- here("input", "Sablefish_Inputs.xlsx")

# Get number of OM and EM scenarios
n_OM_scen <- length(om_scenarios$OM_Scenarios)


# Operating Model ---------------------------------------------------------

for(n_om in 1:n_OM_scen) {
  
  ### Operating Model Loop + Set Up -------------------------------------------
  
  # Create file directory to save model outputs
  om_path <- here("output", "OM_Scenarios", om_scenarios$OM_Scenarios[n_om])
  dir.create(om_path)
  
  # Pre-Processing
  F_type <- c(om_scenarios$Fl_1_Ftype[n_om], om_scenarios$Fl_2_Ftype[n_om]) # Get Fishing Mortality Pattern
  Start_F <- c(om_scenarios$Fl_1_StartF[n_om], om_scenarios$Fl_2_StartF[n_om])  # Get Starting Fs
  Max_Rel_F_M <- c(om_scenarios$Fl_1_MaxRelM[n_om], om_scenarios$Fl_2_MaxRelM[n_om]) # Get max F relative to M
  Input_Fish_N_Max <- c(om_scenarios$Fl_1_InputN_Max[n_om], om_scenarios$Fl_2_InputN_Max[n_om]) # Get max Input N for comps
  Input_N_Fish_Fixed <- c(om_scenarios$Fl_1_InputN_Fixed[n_om], om_scenarios$Fl_2_InputN_Fixed[n_om]) # Get min Input N for comps
  Input_N_Srv <- om_scenarios$Srv_InputN_Max[n_om] # Survey Input N
  Srv_CV <- om_scenarios$Srv_CV[n_om] # Survey Index CV
  Fish_Selex_Opt <- c(om_scenarios$Fl_1_Selex[n_om], om_scenarios$Fl_2_Selex[n_om]) # Fishery Selectivity Options
  
  # Selectivity parameters for Fleet 1 Males, followed by Females
  Fish_Fleet1_SelPars <- c(om_scenarios$Fl_1_Slx_Par1_M[n_om], 
                           om_scenarios$Fl_1_Slx_Par2_M[n_om],
                           om_scenarios$Fl_1_Slx_Par1_F[n_om], 
                           om_scenarios$Fl_1_Slx_Par2_F[n_om])
  
  # Selectivity parameters for Fleet 2 Males, followed by Females
  Fish_Fleet2_SelPars <- c(om_scenarios$Fl_2_Slx_Par1_M[n_om], 
                           om_scenarios$Fl_2_Slx_Par2_M[n_om],
                           om_scenarios$Fl_2_Slx_Par1_F[n_om], 
                           om_scenarios$Fl_2_Slx_Par2_F[n_om])  
  # Simulate data here
  oms <- simulate_data(fxn_path = fxn_path, 
                       spreadsheet_path = lh_path, 
                       
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
                       yr_chng = c(5, 5), 
                       yr_chng_end = c(10, 10),
                       fish_likelihood = "multinomial",
                       Input_Fish_N_Max = Input_Fish_N_Max, 
                       fish_CV = c(0.1, 0.1), # not really used
                       Input_N_Fish_Time = "F_Vary", 
                       Input_N_Fish_Fixed = Input_N_Fish_Fixed, # min of Neff comps when allowed to vary
                       catch_CV = c(0.0, 0.0), 
                       q_Mean_Fish = c(0.05, 0.05), # not really used
                       
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
                       srv_pars = list(Srv_Fleet1 = matrix(data = c(2, 0.85, 5, 0.3), 
                                                           nrow = 2, byrow = TRUE))) # survey fleet 1
  
  # Save as RData file
  save(oms, file = here(om_path, paste(om_scenarios$OM_Scenarios[n_om], ".RData", sep = "")))
  
  # Plot OM here
  plot_OM(path = here("output", "OM_Scenarios", om_scenarios$OM_Scenarios[n_om]), 
          file_name = "OM_Plots.pdf")
  
} # end n_om loop



# Estimation Model --------------------------------------------------------


