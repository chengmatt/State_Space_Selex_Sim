# Purpose: To specify sex ratios for females and males
# Creator: Matthew LH. Cheng
# Date: 11/10/22


#' @param f_ratio Female sex ratios
#' @param m_ratio Male sex ratios
#' If this is a single sex, it doesn't matter what the sex ratios are set to, it'll automatically be 1

specify_sex <- function(f_ratio, m_ratio) {
  
  if(sum(f_ratio, m_ratio) != 1 & n_sex != 1) stop("Sex ratios do not sum to 1")
  
  # If this is single sex
  if(n_sex == 1) {
    sex_ratio <- matrix(1, nrow = n_years, ncol = n_sex, byrow = TRUE)
  } else {
    # Create matrix and ouput to environment - first col = female, second col = male
    sex_ratio <- matrix(c(f_ratio, m_ratio), nrow = n_years, ncol = n_sex, byrow = TRUE)
  }

  sex_ratio <<- sex_ratio
  
} 
