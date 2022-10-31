# Purpose: To feed the outputs of our OM into WHAM format!
# Creator: Matthew LH. Cheng
# Date: 10/31/22

feed_WHAM_input <- function(n_fleets, n_indices, Catch_CV_Val, Selex_Opt = "Asymp",
                            n_indices) {
  
  require(tidyverse)
  
  # The structure of this is taken from https://timjmiller.github.io/wham/reference/prepare_wham_input.html (basic_info section)
  
  # Note: A lot of these objects are getting called in from the environment
  # when we first specify the OM!
  
  # Create a list to store objects in
  input <- list()
  
  # Vector of ages - the last being the plus group
  input$ages <- ages 
  
  # Specify years within our assessment model - Our model start year
  # will be when fishing first starts - the terminal year of the OM
  # N years - 1 because we have no catch or sampling events in the last year
  input$years <- Fish_Start_yr:(n_years-1)
  
  # Number of fishing fleets
  input$n_fleets <- n_fleets
  
# Catch at age munging --------------------------------------------------
  
  # Specify the aggregated catch for each fleet rows = length(years), col = number of fleets
  agg_catch <- matrix(nrow = length(input$years), ncol = input$n_fleets)
  
  # Now munge the catch at age output into an aggregate catch
  catch_df <- melt(Catch_at_age) %>% 
    drop_na() %>% # drop nas in the last year
    rename(Year = Var1, Age = Var2, # Rename varialbes
           Sim = Var3, Catch = value) %>% 
    mutate(Year = parse_number(as.character(Year)),
           Sim = parse_number(as.character(Sim))) %>% # Parse number for year and simulation
    group_by(Year, Sim) %>% # summarize and aggregate catch
    summarize(Catch = sum(Catch, na.rm = TRUE)) %>% 
    filter(Year >= Fish_Start_yr) # Filter out to the start of the fishing year
  
    # Now stick that into the matrix and put it into the input list
    agg_catch[, 1] <- catch_df$Catch
    input$agg_catch <- agg_catch
    
    
# Specify Composition Data in Numbers for Fishery (only 1 fleet for now) -------------------------------------
    
    # Create array to store values in Rows = Fleets, Columns  = Years, 3D = Ages
    # comp_paa[fleet, year, age]
    catch_paa <- array(dim = c(input$n_fleets, length(input$years), length(input$ages)))
    
    # Munge the comp dataframe into the correct format
    fish_comps <- melt(Fish_Age_Comps) %>% 
      drop_na() %>% # drop na values
      rename(Year = Var1, Age = Var2,
             Sim = Var3, Count = value) %>% # Parse number for year, age, and simulation
      mutate(Year = parse_number(as.character(Year)),
             Sim = parse_number(as.character(Sim)),
             Age = parse_number(as.character(Age)))
    
    # Create a vector of years and ages we can loop through
    year_vec <- unique(fish_comps$Year)
    age_vec <- unique(fish_comps$Age)
    
    # Now, feed it into the catch_paa array
    for(y in 1:ncol(catch_paa)) {
      
      for(a in 1:length(input$ages)) { 
        
        # Filter the dataset to a given year and age and stick it into the dataframe
        catch_paa[1, y,a] <- fish_comps[fish_comps$Year == year_vec[y] &  fish_comps$Age == age_vec[a], ]$Count
        
        # print(a)
      } # end age loop
      
    } # end year loop
    
    # Put this into the list
    input$catch_paa <- catch_paa
    
# Catch CV ----------------------------------------------------------------
    
    # Create matrix of catch CVs Rows = years, Columns = Number of fleets
    catch_cv <- matrix(nrow = length(input$years), ncol = input$n_fleets)
    
    # Input CV values in the matrix
    catch_cv[,] <- Catch_CV_Val

    # Now, put this into a list
    input$catch_cv <- catch_cv
    

# Effective Sample Size for Proportion at Age -----------------------------

    # Create matrix for effective sample sizes for each fleets age comps
    # rows = years, columns = fleet
    catch_Neff <- matrix(nrow = length(input$years), ncol = input$n_fleets)
    
    # Stick the effective sample size into the matrix and put matrix into a list
    catch_Neff[,] <- N_eff_fish
    input$catch_Neff <- catch_Neff
    
    # Next, create a matrix to specify whehter or not to use the comps for a given year
    # For now, let's just specify 1s - which is to use comps
    # Rows = years, Columns = Fleets
    use_catch_paa <- matrix(nrow = length(input$years), ncol = input$n_fleets)
    
    # Specify 1s and stick into list
    use_catch_paa[,] <- 1
    input$use_catch_paa <- use_catch_paa
    

# Selectivity options for each year and fleet -----------------------------
    
    # Age specific
    if(Selex_Opt == "AgeSpec") pointer <- 1
    # Asymptotic
    if(Selex_Opt == "Asymp") pointer <- 2
    # Double logistic
    if(Selex_Opt == "DubLogist") pointer <- 3
    # Decreasing logistic
    if(Selex_Opt == "DecLogist") pointer <- 4
    # Dome shaped
    if(Selex_Opt == "Dome") pointer <- 5 

    # Create matrix to store selex pointers in
    selblock_pointer_fleets <- matrix(nrow = length(input$years), ncol = input$n_fleets)
    
    # Specify the pointer/selex option here and stick this matrix into our input list
    selblock_pointer_fleets[,] <- pointer
    input$selblock_pointer_fleets <- selblock_pointer_fleets
    

# Specify Fs to initialize model ------------------------------------------

    # Create matrix to store values in 
    F_mat <- matrix(nrow = length(input$years), ncol = input$n_fleets)
    
    # Fill this in with random normal draws and stick into our list
    F_mat[,] <- abs(rnorm(n = length(input$years), mean = 0, sd = 0.3))
    input$`F` <- F_mat # "Random" starting values
  

# Indices of abundance ----------------------------------------------------

  # Specify number of indices we want to model in the assessment
  input$n_indices <- n_indices
  
  
  
  
}