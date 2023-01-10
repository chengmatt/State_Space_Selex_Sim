# Creator: Matthew LH. Cheng
# Purpose: To preapre inputs for the estimation model in TMB
# Date 1/10/23

#' Title Prepares data inputs for the EM 
#'
#' @param ages vector of ages
#' @param years vector of years
#' @param n_sexes integer of sexes
#' @param n_fleets integer of fleets
#' @param n_fish_comps integer for fishery comps (needs to match number of fleets)
#' @param n_srv_comps integer for survey comps (needs to match number of surveys)
#' @param n_fish_indices integer for number of fishery indices fit
#' @param n_srv_indices integer for number of survey indices fit
#' @param Fish_Start_yr vector of fishery start years
#' @param catch_cv vector of catch CVs
#' @param F_Slx_Blocks matrix of selectivity blocks, where it is encoded as:
#' matrix(c(0), nrow = length(years), ncol = n_fleets), with unique numbers representing unique time-blocks
#' @param S_Slx_Blocks same as above
#' @param use_catch matrix of whether or not to use catch data 0 = don't use, 1 = use, where it is encoded as:
#' use_catch <- matrix(1, nrow = n_years, ncol = n_fish_fleets)
#' @param use_fish_index same as above
#' @param use_srv_index same as above
#' @param use_fish_comps array of whether or not to use comp data 0 = don't use, 1 = use, where it is encoded as:
#' use_fish_comps <- array(1, dim = c(n_years, n_fish_comps, n_sex))
#' @param use_srv_comps same as abovce
#' @param rec_model recruitment model == 0 (mean recruitment)
#' @param F_Slx_model Fishery selex model ==0, logistic, == 1 gamma, ==2 double logistic
#' @param S_Slx_model same as above
#' @param sim simulation indexing
#'
#' @return
#' @export
#'
#' @examples
prepare_EM_input <- function(ages, 
                          years,
                          n_sexes, 
                          n_fleets, 
                          n_fish_comps,
                          n_srv_comps,
                          n_fish_indices, 
                          n_srv_indices, 
                          Fish_Start_yr,
                          catch_cv,
                          F_Slx_Blocks,
                          S_Slx_Blocks,
                          use_catch,
                          use_fish_index,
                          use_srv_index, 
                          use_fish_comps,
                          use_srv_comps,
                          F_Slx_model,
                          S_Slx_model = 0,
                          rec_model = 0,
                          sim) {
  
  # Make input into a list 
  input <- list()
  

# Catch -------------------------------------------------------------------

  if(n_fleets == 1) { # if single fleet - sum across or leave as is
    obs_catches <- apply( as.matrix(Catch_agg[(Fish_Start_yr[1]:(n_years-1)),,sim]), 1, FUN = sum)
  } else{ # multi fleet - leave as is
    obs_catches <- matrix(Catch_agg[(Fish_Start_yr[1]:(n_years-1)),,sim], nrow = length(years), ncol = n_fleets)
  } # end else
  

# Fishery Age Comps -------------------------------------------------------

  obs_fish_age_comps <- array(data = 0, dim = c(length(years), length(ages), n_fleets, n_sexes))
  
  if(n_fleets == 1) { # one fleet = sum across
    
    # Loop through to sum - index the 3rd dimension to get fleets
    for(f in 1:dim(Fish_Age_Comps)[3]) {
      for(s in 1:n_sexes) {
        # Filter to save as an object
        fish_age_comps <-  Fish_Age_Comps[Fish_Start_yr[1]:(n_years - 1),,f,s,sim]
        # Increment comps - fixing fleet index to 1 here
        obs_fish_age_comps[,,1,s] <- obs_fish_age_comps[,,1,s] + fish_age_comps
      } # end s loop
    } # end f loop
    
    # Now, apply the proportion function over a single fleet
    obs_fish_age_comps <- array(t(apply(obs_fish_age_comps, MARGIN = 1, FUN=function(x) { x/sum(x) })),
                                dim = c(n_years, n_fleets, n_sexes))
    
    # Effective Sample Sizes
    obs_fish_age_Neff <- as.matrix(apply(fish_Neff[Fish_Start_yr[1]:(n_years - 1),], MARGIN = 1, FUN = sum), 
                                   nrow = length(years), ncol = n_fleets)
    
  } else{ # more than one fleet here
    
    for(f in 1:n_fleets) { # needs to loop through transpose and apply for easy retnetion of array dimensions
      for(s in 1:n_sexes) {
        obs_fish_age_comps[,,f,s] <- t(
          apply(Fish_Age_Comps[Fish_Start_yr[1]:(n_years - 1),,f,s,sim], MARGIN = 1, 
                FUN=function(x) { x/sum(x) })
        )
      } # end s loop
    } # end f loop
    
    # Effective Sample Sizes
    obs_fish_age_Neff <- matrix(fish_Neff[Fish_Start_yr[1]:(n_years - 1),], 
                                    nrow = length(years), ncol = n_fleets)
    
  } # end else
  

# Survey Age Comps --------------------------------------------------------
# Right now, this is only able to accomdate one single survey fleet
  
  # Apply function to get proportions and munge into matrix
  obs_srv_age_comps <- array(
    t(apply(Survey_Age_Comps[Fish_Start_yr[1]:(n_years - 1),,,,sim], MARGIN = 1, 
            FUN=function(x) { x/sum(x) })),
    dim = c(length(years), length(ages), n_srv_comps, n_sexes)
  )
  
  # Get survey age neffs
  obs_srv_age_Neff <- matrix(srv_Neff[Fish_Start_yr[1]:(n_years - 1),], 
                             nrow = length(years), ncol = n_srv_comps)
  

# Abundance Indices -------------------------------------------------------

  # Fishery index
  obs_fish_indices <-as.matrix( Fishery_Index_Agg[Fish_Start_yr[1]:(n_years - 1),,sim], 
                                nrow = length(years), ncol = n_fish_indices)
  
  # Survey index
  obs_srv_indices <- as.matrix(Survey_Index_Agg[Fish_Start_yr[1]:(n_years - 1),,sim], 
                               nrow = length(years), ncol = n_srv_indices)
  
  
# Biological inputs -------------------------------------------------------

  # Weight at age
  WAA <- array(wt_at_age[Fish_Start_yr[1]:(n_years),,1,sim], 
               dim = c(length(years), length(ages), n_sexes))
  
  # Maturity at age
  MatAA <- array(mat_at_age[Fish_Start_yr[1]:(n_years),,1,sim]  ,
                 dim = c(length(years), length(ages), n_sexes))
  
  # Sex Ratios
  Sex_Ratio <- as.vector(Sex_Ratio[1:n_sexes])
  

# CV inputs ---------------------------------------------------------------
  
  # Fishery, survey, and catch CV
  fish_cv <- as.vector(fish_CV)
  srv_cv <- as.vector(srv_CV)
  catch_cv <- as.vector(catch_cv)

# Selectivity blocks ------------------------------------------------------

  F_Slx_Blocks <- F_Slx_Blocks # fishery blocks
  S_Slx_Blocks <- S_Slx_Blocks # survey blocks

# Data Indicators ---------------------------------------------------------

  use_fish_index <- use_fish_index # fishery index
  use_srv_index <- use_srv_index # survey index
  use_fish_comps <- use_fish_comps # fishery comps
  use_srv_comps <- use_srv_comps  # survey comps

  # Input these data into a list object
  input$ages <- ages
  input$years <- years
  input$n_sexes <- n_sexes
  input$n_fleets = n_fleets
  input$n_fish_comps = n_fish_comps
  input$n_srv_comps = n_srv_comps
  input$n_fish_indices = n_fish_indices
  input$n_srv_indices = n_srv_indices
  input$obs_catches <- obs_catches
  input$obs_fish_age_comps <- obs_fish_age_comps
  input$obs_fish_age_Neff <- obs_fish_age_Neff
  input$obs_srv_age_comps <- obs_srv_age_comps
  input$obs_srv_age_Neff <- obs_srv_age_Neff
  input$obs_fish_indices <- obs_fish_indices
  input$obs_srv_indices <- obs_srv_indices
  input$WAA <- WAA
  input$MatAA <- MatAA
  input$Sex_Ratio <- Sex_Ratio
  input$fish_cv <- fish_cv
  input$srv_cv <- srv_cv
  input$catch_cv <- catch_cv
  input$F_Slx_Blocks <- F_Slx_Blocks
  input$S_Slx_Blocks <- S_Slx_Blocks
  input$use_catch <- use_catch
  input$use_fish_index <- use_fish_index
  input$use_srv_index  <- use_srv_index 
  input$use_fish_comps <- use_fish_comps
  input$use_srv_comps  <- use_srv_comps 
  input$rec_model <- rec_model
  input$S_Slx_model <- S_Slx_model
  input$F_Slx_model <- F_Slx_model
  
  return(input)
  
}
