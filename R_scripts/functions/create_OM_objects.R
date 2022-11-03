# Purpose: To create objects to hold for the opearting model and observation models
# Creator: Matthew LH. Cheng
# Date: 10/29/22

#' @param n_years Number of years for the simulation
#' @param ages Vector of ages
#' @param n_sims Number of simulations we want to run

create_OM_objects <<- function(n_years, ages, n_sims) {
  
  # Numbers at age
  N_at_age <<- array(dim = c(n_years, length(ages), n_sims), 
                    dimnames = list( # Set up dimensions names
                      c(paste("Year", 1:n_years, sep = "_")), # Years 
                      c(paste("Age", ages, sep = "_")), # Years 
                      c(paste("Sim",1:n_sims)) # Simulation 
                    ))
  
  # Biomass at age
  Biom_at_age <<- array(dim = c(n_years, length(ages), n_sims), 
                       dimnames = list( # Set up dimensions names
                         c(paste("Year", 1:n_years, sep = "_")), # Years 
                         c(paste("Age", ages, sep = "_")), # Years 
                         c(paste("Sim",1:n_sims)) # Simulation 
                       ))
  
  # SSB across years
  SSB <<- array(dim = c(n_years, n_sims), 
               dimnames = list( # Set up dimensions names
                 c(paste("Year", 1:n_years, sep = "_")), # Years 
                 c(paste("Sim",1:n_sims)) # Simulation 
               ))
  
  # Total Recruits across years
  rec_total <<- array(dim = c(n_years, n_sims), 
                     dimnames = list( # Set up dimensions names
                       c(paste("Year", 1:n_years, sep = "_")), # Years 
                       c(paste("Sim",1:n_sims)) # Simulation 
                     ))
  
  #Recruitment deviates
  rec_devs <<- array(dim = c(n_years, n_sims),
                    dimnames = list(
                      c(paste("Year", 1:n_years, sep = "_")),
                      c(paste("Sim",1:n_sims))
                    ))
  
  # Catch at age across years
  Catch_at_age <<- array(dim = c(n_years, length(ages), n_sims), 
                        dimnames = list( # Set up dimensions names
                          c(paste("Year", 1:n_years, sep = "_")), # Years 
                          c(paste("Age", ages, sep = "_")), # Years 
                          c(paste("Sim",1:n_sims)) # Simulation 
                        ))
  

# Observation model -------------------------------------------------------

  # Fishery Age Comps across ages and years
  Fish_Age_Comps <<- array(dim = c(n_years, length(ages), n_sims), 
                          dimnames = list( # Set up dimensions names
                            c(paste("Year", 1:n_years, sep = "_")), # Years 
                            c(paste("Age", ages, sep = "_")), # Years 
                            c(paste("Sim",1:n_sims)) # Simulation 
                          ))
  
  # Fishery Index across years
  Fishery_Index <<- array(dim = c(n_years, n_sims), 
                         dimnames = list( # Set up dimensions names
                           c(paste("Year", 1:n_years, sep = "_")), # Years 
                           c(paste("Sim",1:n_sims)) # Simulation 
                         ))
  
  Survey_Age_Comps <<- array(dim = c(n_years, length(ages), n_sims), 
                          dimnames = list( # Set up dimensions names
                            c(paste("Year", 1:n_years, sep = "_")), # Years 
                            c(paste("Age", ages, sep = "_")), # Years 
                            c(paste("Sim",1:n_sims)) # Simulation 
                          ))
  
  # Fishery Index across years
  Survey_Index <<- array(dim = c(n_years, n_sims), 
                         dimnames = list( # Set up dimensions names
                           c(paste("Year", 1:n_years, sep = "_")), # Years 
                           c(paste("Sim",1:n_sims)) # Simulation 
                         ))
  
  
  # Input parameters here ---------------------------------------------------

  # Natural mortality at age and across years
  Mort_at_age <<- array(dim = c(n_years, length(ages), n_sims), 
                       dimnames = list( # Set up dimensions names
                         c(paste("Year", 1:n_years, sep = "_")), # Years 
                         c(paste("Age", ages, sep = "_")), # Years 
                         c(paste("Sim",1:n_sims)) # Simulation 
                       ))
  
  # Maturity at age across years
  mat_at_age <<- array(dim = c(n_years, length(ages), n_sims), 
                      dimnames = list( # Set up dimensions names
                        c(paste("Year", 1:n_years, sep = "_")), # Years 
                        c(paste("Age", ages, sep = "_")), # Years 
                        c(paste("Sim",1:n_sims)) # Simulation 
                      ))
  
  # Weight at age across years
  wt_at_age <<- array(dim = c(n_years, length(ages), n_sims), 
                     dimnames = list( # Set up dimensions names
                       c(paste("Year", 1:n_years, sep = "_")), # Years 
                       c(paste("Age", ages, sep = "_")), # Years 
                       c(paste("Sim",1:n_sims)) # Simulation 
                     ))
  
  
  # Fishing mortality across years
  fish_mort <<- array(dim = c(n_years, n_sims),
                     dimnames = list(
                       c(paste("Year", 1:n_years, sep = "_")),
                       c(paste("Sim",1:n_sims))
                     ))
  
  # Fishery selectivity across ages and years
  Fish_selex_at_age <<- array(dim = c(n_years, length(ages), n_sims), 
                             dimnames = list( # Set up dimensions names
                               c(paste("Year", 1:n_years, sep = "_")), # Years 
                               c(paste("Age", ages, sep = "_")), # Years 
                               c(paste("Sim",1:n_sims)) # Simulation 
                             ))
  
  # Survey selectivity across ages and years
  Surv_selex_at_age <<- array(dim = c(n_years, length(ages), n_sims), 
                              dimnames = list( # Set up dimensions names
                                c(paste("Year", 1:n_years, sep = "_")), # Years 
                                c(paste("Age", ages, sep = "_")), # Years 
                                c(paste("Sim",1:n_sims)) # Simulation 
                              ))
  
  # Catchability across time for the fishery
  q_Fish <<- array(dim = c(n_years, n_sims), 
             dimnames = list( # Set up dimensions names
               c(paste("Year", 1:n_years, sep = "_")), # Years 
               c(paste("Sim",1:n_sims)) # Simulation 
             ))
  
  # Catchability across time for the survey
  q_Surv <<- array(dim = c(n_years, n_sims), 
                  dimnames = list( # Set up dimensions names
                    c(paste("Year", 1:n_years, sep = "_")), # Years 
                    c(paste("Sim",1:n_sims)) # Simulation 
                  ))
  
} # end function
