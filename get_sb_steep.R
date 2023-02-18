# Purpose: To get steepness values for sablefish using SR model outputs from Goethel et al. 2021
# Creator: Matthew LH. Cheng
# Date: 2/17/22


# Set Up ------------------------------------------------------------------

library(here)
library(bbmle)
library(tidyverse)

# Read in SR data
sr_dat <- read.csv(here("input", "Sablefish_SR_Outputs.csv")) 


# Function to minimize over for BH SR (Francis 1992 Parameterization)
nll_bh_sr <- function(SSB_y, Rec_y, ln_R0, logit_h, ln_sigma, SSB0) {
  
  # Exponenitate parameters
  R0 <- exp(ln_R0)
  h <- 0 + (1 - 0)/(1 + exp(-logit_h))
  sigma <- exp(ln_sigma)
  
  # Get Observed SSB and Rec
  obs_SSB <- SSB_y
  obs_Rec <- Rec_y
  
  # Get predicted recruitment
  BH_first_part <- 4 * h * R0 * obs_SSB
  BH_sec_part <- (SSB0 * (1 - h)) + obs_SSB * ((5*h) - 1)
  
  # Now calculate BH parameterization
  pred_rec <- BH_first_part / BH_sec_part
  
  plot(SSB_y, pred_rec)

  # Calculate log likelihood
  logLike <- dnorm(x=log(obs_Rec), mean=log(pred_rec), sd=sigma, log=TRUE)

  # Calculate total negative log-likelihood
  NLL <- -1*sum(logLike)
  
} # end function


# Fit model now
mod_fit <- mle2(nll_bh_sr, 
                start=list(logit_h = 0.3, 
                           ln_R0 = 2.8587),
                data=list(SSB_y = sr_dat$SSB, Rec_y = sr_dat$Rec,
                          SSB0 = 409.7684, ln_sigma = 1.2),
                optimizer = "nlminb",
                method="Nelder-Mead",
                control=list(maxit=1e8))

# Look at coefficients
summary(mod_fit)
