# Purpose: To specify data quality and quantity scenarios for the fishery
# and the survey, start year for the fishery, and start year for the survey
# Creator: Matthew LH. Cheng
# Date: 10/30/22

#' @param Fish_Start_yr Fishery start year
#' @param Surv_Start_yr Survey start year
#' @param n_years Number of years within the simulation
#' @param Fish_freq How frequently the fishery samples comps and index
#' @param Surv_freq How frequently the survey samples comps and index
#' @param Scenario Data Quality and Quantity Scenarios (Low, Medium, High)

fish_surv_data_scenarios <- function(Fish_Start_yr, Surv_Start_yr, n_years,
                                     Fish_freq, Surv_freq, Scenario) {
  
  # Create objects for Fishery Start Year
  Fish_Start_yr <<- Fish_Start_yr
  
  # Survey Start Year
  Surv_Start_yr <<- Surv_Start_yr
  
  # Frequency of fishing years that data are collected
  Fish_yrs <<- seq(Fish_Start_yr, n_years, Fish_freq)
  
  # Frequency of survey years that data are collected
  Surv_yrs <<- seq(Surv_Start_yr, n_years, Surv_freq)
  
  # Specify Data quality and quantity scenarios
  if(Scenario == "Low") {
    Fishery_CV <<- 0.3
    Survey_CV <<- 0.25
    N_eff_fish <<- 30
    N_eff_surv <<- 20
  } # Low data quality and quantity scenario
  
  if(Scenario == "Medium") {
    Fishery_CV <<- 0.2
    Survey_CV <<- 0.2
    N_eff_fish <<- 60
    N_eff_surv <<- 40
  } # Medium data quality and quantity scenario
  
  if(Scenario == "High") {
    Fishery_CV <<- 0.1
    Survey_CV <<- 0.1
    N_eff_fish <<- 200
    N_eff_surv <<- 150
  } # Low data quality and quantity scenario
  
} # end function
