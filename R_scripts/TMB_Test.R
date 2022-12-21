# Purpose: To test TMB Stock assessment Model
# Creator: Matthew LH. Cheng


# Set up ------------------------------------------------------------------

library(here)
library(tidyverse)
library(TMB)

# Need to set wd in source folder
setwd(here("src"))

# Compile and run TMB model -----------------------------------------------

compile("EM.cpp") # Compile .cpp file
dyn.load(dynlib("EM")) # Load in .cpp



# Load in data ------------------------------------------------------------
sim <- 1

ages <- 1:30
years <- Fish_Start_yr:(n_years - 1)
n_sexes <- n_sex
n_fleets <- n_fish_fleets
n_fish_indices <- 1
n_srv_indices <- 1

# Calculate catches 
obs_catches <- melt(Catch_at_age) %>%
  drop_na() %>% # drop nas in the last year
  rename(Year = Var1, Age = Var2, # Rename varialbes
         Fleet = Var3, Sex = Var4, Sim = Var5, Catch = value) %>%
  mutate(Year = parse_number(as.character(Year)),
         Sim = parse_number(as.character(Sim))) %>%  # Parse number for year and simulation
  filter(Year >= Fish_Start_yr[1],
         Sim == 1)
obs_catches <- matrix(with(obs_catches, tapply(Catch, list(Year), FUN = sum)))
# obs_catches <-obs_catches + rnorm(length(obs_catches), 0, 1)

# Observed fishery age comps
obs_fish_age_comps <- t(apply(Fish_Age_Comps[Fish_Start_yr:(n_years - 1),,,,sim], MARGIN = 1, 
                              FUN=function(x) { x/sum(x) }))
obs_fish_age_Neff <- fish_Neff

# Observed survey age comps
obs_srv_age_comps <- t(apply(Survey_Age_Comps[Fish_Start_yr:(n_years - 1),,,,sim], MARGIN = 1, 
                              FUN=function(x) { x/sum(x) }))
obs_srv_age_Neff <- srv_Neff

# Observed fishery indices
obs_fish_indices <-as.matrix( Fishery_Index_Agg[Fish_Start_yr:(n_years - 1),,sim])
obs_srv_indices <- as.matrix(Survey_Index_Agg[Fish_Start_yr:(n_years - 1),,sim])

WAA <- wt_at_age[Fish_Start_yr:(n_years),,1,sim]
MatAA <- mat_at_age[Fish_Start_yr:(n_years),,1,sim]
Sex_Ratio <- c(1,0)

catch_cv <- 0.3
fish_cv <- 0.1
srv_cv <- 0.1

# Testing
F_Slx <- Fish_selex_at_age[Fish_Start_yr:(n_years-1),,,,]
S_Slx <- Surv_selex_at_age[Fish_Start_yr:(n_years-1),,,,]
F_Mort <- fish_mort[Fish_Start_yr:(n_years-1),,]
init_model <- 0
rec_model <- 0

Init_N_at_age <- N_at_age[(Fish_Start_yr-1),,,]

# TMB Section -------------------------------------------------------------


# Fill in list for data
data <- list( ages = ages, years = years,
             n_sexes = n_sexes, n_fleets = n_fleets,
             n_fish_indices = n_fish_indices, n_srv_indices = n_srv_indices,
             obs_catches = obs_catches, obs_fish_age_comps = obs_fish_age_comps,
             obs_fish_age_Neff = obs_fish_age_Neff, obs_srv_age_comps = obs_srv_age_comps,
             obs_srv_age_Neff = obs_srv_age_Neff, obs_fish_indices =  obs_fish_indices,
             obs_srv_indices = obs_srv_indices, WAA = array(WAA, dim = c(32, 30, 1)), 
             MatAA = array(MatAA, dim = c(32, 30, 1)),
             Sex_Ratio = Sex_Ratio, init_model = 0, rec_model = 0, catch_cv = 0.05,
             fish_cv = fish_cv, srv_cv = srv_cv, F_Slx = array(F_Slx, dim = c(31, 30, 1, 1)),
             S_Slx = array(S_Slx, dim = c(31, 30, 1, 1)), F_Mort = as.matrix(F_Mort),
             Init_N_at_age = as.matrix(Init_N_at_age)
             )

# Define parameter inits here
parameters <- list(ln_R0 = 20, ln_SigmaRec = 1.2,ln_MeanRec = log(2),
                   ln_M = log(0.1), 
                   ln_F_y = as.matrix(c(rnorm(length(Fish_Start_yr:(n_years - 1)), 1, 0.05))),
                   log_q_fish = as.matrix(log(rep(0.03, n_fish_fleets))), 
                   log_q_srv = as.matrix(log(rep(0.05, n_srv_fleets))))

compile("EM.cpp") # Compile .cpp file
dyn.load(dynlib("EM")) # Load in .cpp

# Make ADFun
mod <- MakeADFun(data, parameters, DLL="EM")
mod$opt <- stats::nlminb(mod$par, mod$fn, mod$gr, 
                     control = list(iter.max = 1500, eval.max = 1000))
sd_rep <- TMB::sdreport(mod)

plot(exp(sd_rep$par.fixed[3:33]))
lines(fish_mort[(Fish_Start_yr:(n_years - 1)),,], col = "red")
sd_rep$gradient.fixed


