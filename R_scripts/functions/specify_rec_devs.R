# Purpose: To specify type of recruitment deviations
# Creator: Matthew LH. Cheng
# Date: 10/31/22


#' @param Rec_Dev_Type Recruitment deviation type - options are iid or Auto_Cor
#' @param rho_rec If rec dev type is autocorrelated, then specify a correlation parameter
#' that is between 0 and 1

specify_rec_devs <- function(Rec_Dev_Type, rho_rec = 0) {
  
  # If correaltion parameter is not between 0 and 1
  if(rho_rec < 0 & rho_rec > 1 & Rec_Dev_Type == "Auto_Cor") stop("Correlation parameter is not between 0 and 1, please respecify!")
  
  if(Rec_Dev_Type == "iid") {
    
    for(sim in 1:n_sims) {
      
      rec_devs[,sim] <- rnorm(n_years, 0, sigma_rec)
      
    } # end simulation loop
     
  } # if we want random iid recrutiment deviates
  
  if(Rec_Dev_Type == "Auto_Cor") {
    
    for(sim in 1:n_sims) {
      
      # Create deviations here first
      devs <- rnorm(n_years, 0, sigma_rec)
      
      # Put these generated deviations in our rec devs dataframe
      rec_devs[1,sim] <- devs[1]
      
      for(y in 2:n_years) {

        # Next, calculate the true rec devs given autocorrelation parameter rho_rec
        # adding some random noise to it multiplied by some degree of disassociation
        # Degree of diassociation - basically if rho_rec = 1, then it ensures that there aren't 
        # any deviations (because it would be devs * 0, conversely, if rho_rec = 0, then the noise
        # term gets applied to its full extent, because there isn't any correlation)
        
        rec_devs[y,sim] <- ((rho_rec * rec_devs[y-1,sim]) + (devs[y] * sqrt(1 - rho_rec^2)))
        
      } # end years loop

    } # end simulation loop
    
  } # if recruitment deviates are autocorrelated
  
  # Output rec devs into environment here
  rec_devs <<- rec_devs
  
  print(paste("Recruitment deviations are specified as:", Rec_Dev_Type))
  
} # end function

