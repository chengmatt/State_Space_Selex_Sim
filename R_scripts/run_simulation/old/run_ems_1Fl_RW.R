# Purpose: Run EMs that are 1 Fleet and are all random walks
# Creator: Matthew LH. Cheng
# Date 6/12/23


# Set up ------------------------------------------------------------------
rm(list=ls()) # remove objects prior to running

library(here)
library(tidyverse)
library(TMB)
library(doSNOW)
library(parallel)

ncores <- detectCores() 
cl <- makeCluster(ncores - 2)
registerDoSNOW(cl)

# Load in all functions into the environment
fxn_path <- here("R_scripts", "functions")
# Load in all functions from the functions folder
files <- list.files(fxn_path)
for(i in 1:length(files)) source(here(fxn_path, files[i]))

compile_tmb(wd = here("src"), cpp = "EM.cpp")

# Read in OM and EM Scenarios
om_scenarios <- readxl::read_excel(here('input', "OM_EM_Scenarios_v2.xlsx"), sheet = "OM")
em_scenarios <- readxl::read_excel(here('input', "OM_EM_Scenarios_v2.xlsx"), sheet = "EM_1Fl_TI_Blk")

# Read in spreadsheet for life history parameters
lh_path <- here("input", "Sablefish_Inputs.xlsx")

# Get number of OM and EM scenarios
n_OM_scen <- length(om_scenarios$OM_Scenarios)
n_EM_scen <- length(em_scenarios$EM_Scenario)


