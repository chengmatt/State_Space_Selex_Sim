# Purpose: Plotting outputs of operating model
# Date: 10/29/22
# Creator: Matthew LH. Cheng

#' @param path Path to output folder

plot_OM <- function(path) {
  
  pdf(here(path, "OM_output_figs.pdf"), height = 15, width = 10)
  
  require(tidyverse)

# Recruitment -------------------------------------------------------------

  # Plot Recruitment scenarios
  rec_df <- melt(rec_total)
  names(rec_df) <- c("Year", "Sim", "Rec")
  
  # Total recruitment across time and simulations
print(  ggplot(rec_df, aes(x = as.numeric(Year), y = Rec, color = Sim, group = 1)) +
          geom_line() +
          facet_wrap(~Sim, ncol = 1) +
          theme_bw() +
          theme(legend.position = "none") +
          labs(x = "Year",y = "Recruitment Total") )
  
# SSB ---------------------------------------------------------------------

  # Plot SSB
  ssb_df <- melt(SSB)
  names(ssb_df) <- c("Year", "Sim", "SSB")
  
print(
  ggplot(ssb_df, aes(x = as.numeric(Year), y = SSB, color = Sim, group = 1)) +
    geom_line(size = 1.5) +
    facet_wrap(~Sim, ncol = 1) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "Year",y = "SSB") 
)
  

# Recruitment ~ SSB -------------------------------------------------------

  spr <- data.frame(SSB = ssb_df$SSB, Rec = rec_df$Rec, Sim = rec_df$Sim)

print(
  
  ggplot(spr, aes(x = SSB, y = Rec)) +
    geom_point() +
    geom_smooth() +
    facet_wrap(~Sim, ncol = 1, scales = "free") +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "SSB",y = "Recruitment Total") 
)
  

# Biomass and Numbers at Age ----------------------------------------------
  
  # biomass at age
  biom_df <- melt(Biom_at_age)
  names(biom_df) <- c("Year", "Age", "Sim", "Biomass")
  
print(
  # Biomass at age
  ggplot(biom_df, aes(x = as.numeric(Year), y = Biomass, color = Age)) +
    geom_line() +
    facet_wrap(~Sim, ncol = 1, scales = "free") +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "Year",y = "Biomass at age")
)
  
print(
  # Biomass aggregated
  biom_df %>% 
    group_by(Year, Sim) %>% 
    summarize(Biomass = sum(Biomass, na.rm = TRUE)) %>% 
    ggplot(aes(x = as.numeric(Year), y = Biomass, color = Sim)) +
    geom_line(size = 1.5) +
    facet_wrap(~Sim, ncol = 1, scales = "free") +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "Year",y = "Biomass")
)
  
  # Plot Numbers at age
  natage <- melt(N_at_age)
  names(natage) <- c("Year", "Age", "Sim", "Numbers")
  
print(
  # Numbers across ages
  ggplot(natage, aes(x = as.numeric(Year), y = Numbers, color = Age)) +
    geom_line() +
    facet_wrap(~Sim, ncol = 1, scales = "free") +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "Year",y = "Numbers at age")
)
  
print(
  # Numbers aggregated by ages
  natage %>% 
    group_by(Year, Sim) %>% 
    summarize(Numbers = sum(Numbers, na.rm = TRUE)) %>% 
    ggplot(aes(x = as.numeric(Year), y = Numbers, color = Sim)) +
    geom_line(size = 1.5) +
    facet_wrap(~Sim, ncol = 1, scales = "free") +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "Year",y = "Numbers")
)
  

# Proportion at Age -------------------------------------------------------
  
print(
  biom_df %>% 
    group_by(Sim, Year) %>% 
    mutate(total = sum(Biomass)) %>% 
    group_by(Age, Sim) %>% 
    summarize(prop_age = mean(Biomass/total, na.rm = TRUE)) %>% 
    ggplot(aes(x = as.numeric(Age), y = prop_age, group = Sim, color = Sim )) +
    geom_line(size = 1.5) +
    facet_wrap(~Sim, ncol = 1, scales = "free") +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "Age",y = "Average proportion at age")
)
  

# Catch at age ------------------------------------------------------------

  catch_df <- melt(Catch_at_age)
  names(catch_df) <- c("Year", "Age", "Sim", "Catch")
  
print(
  catch_df %>% 
    drop_na() %>% 
    group_by(Year, Sim) %>% # Aggregated catch
    summarize(Catch  = sum(Catch, na.rm = T)) %>% 
    ggplot(aes(x = as.numeric(Year), y = Catch, color = Sim)) +
    geom_line(size = 1.5) +
    facet_wrap(~Sim, ncol = 1, scales = "free") +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "Age",y = "Catch")
)
  

# Composition data --------------------------------------------------------

tryCatch({ 
  
  # Fishery Comps
  comps_fish <- melt(Fish_Age_Comps)
  names(comps_fish) <- c("Year", "Age", "Sim", "Count")
  
print(
  ggplot(comps_fish, aes(x = as.numeric(Year), y = Count, fill = Age)) +
    geom_col(position = "stack") +
    facet_wrap(~Sim, ncol = 1, scales = "free") +
    theme_bw() +
    theme(legend.position = "top") +
    labs(x = "Year",y = "Count", title = "Fishery Comps")
)
  
  comps_surv <- melt(Survey_Age_Comps)
  names(comps_surv) <- c("Year", "Age", "Sim", "Count")
  
print(
  ggplot(comps_surv, aes(x = as.numeric(Year), y = Count, fill = Age)) +
    geom_col(position = "stack") +
    facet_wrap(~Sim, ncol = 1, scales = "free") +
    theme_bw() +
    theme(legend.position = "top") +
    labs(x = "Year",y = "Count", title = "Survey Comps") 
)
  

# Indices of abundance ----------------------------------------------------

  # Fishery index
  idx_fish <- melt(Fishery_Index_at_age) 
  names(idx_fish) <- c("Year", "Age", "Sim", "Index")
  
print(
  # Across ages
  ggplot(idx_fish %>% drop_na(), aes(as.numeric(Year), y = Index, color = Age))+
    geom_line(size = 1.1) +
    facet_wrap(~Sim, ncol = 1, scales = "free") +
    theme_bw() +
    theme(legend.position = "top") +
    labs(x = "Year",y = "Biomass index at age", title = "Fishery Index") 
)
  
print(
  # Aggregated
  idx_fish %>% 
    drop_na() %>% 
    group_by(Year, Sim) %>% 
    summarize(Index = sum(Index, na.rm = TRUE)) %>% 
    ggplot(aes(as.numeric(Year), y = Index, color = Sim))+
    geom_line(size = 1.1, alpha = 0.75) +
    geom_point(size = 2) +
    facet_wrap(~Sim, ncol = 1, scales = "free") +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "Year",y = "Biomass index", title = "Fishery Index") 
)
  
  # Survey index
  idx_surv <- melt(Survey_Index_at_age)
  names(idx_surv) <- c("Year", "Age", "Sim", "Index")
  
print(
  # Across ages
  ggplot(idx_surv %>% drop_na(), aes(as.numeric(Year), y = Index, color = Age))+
    geom_line(size = 1.1) +
    facet_wrap(~Sim, ncol = 1, scales = "free") +
    theme_bw() +
    theme(legend.position = "top") +
    labs(x = "Year",y = "Biomass index at age", title = "Survey Index") 
)
  
print(
  # Aggregated
  idx_surv %>% 
    drop_na() %>% 
    group_by(Year, Sim) %>% 
    summarize(Index = sum(Index, na.rm = TRUE)) %>% 
    ggplot(aes(as.numeric(Year), y = Index, color = Sim))+
    geom_line(size = 1.1, alpha = 0.75) +
    geom_point(size = 2) +
    facet_wrap(~Sim, ncol = 1, scales = "free") +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "Year",y = "Biomass index", title = "Survey Index") 
)

  }, error = function(error) {cat("ERROR :",conditionMessage(error), "\n")}) # end try catch statement
  


# Parameterizations -------------------------------------------------------


  ### Selectivity -------------------------------------------------------------

  # Fishery selectivity
  fish_sel_df <- melt(Fish_selex_at_age)
  names(fish_sel_df) <- c("Year", "Age", "Sim", "Selex")
  
  fish_sel_df <- fish_sel_df %>% 
    mutate(Age = parse_number(substr(Age, 5, 10)))
  
print(
  ggplot(fish_sel_df, aes(x = as.numeric(Year), y = Age,
                          fill = Selex)) +
    geom_tile() +
    facet_wrap(~Sim, ncol = 1, scales = "free") +
    scale_fill_viridis_c() +
    theme_bw() +
    theme(legend.position = "top") +
    labs(x = "Year",y = "Age", title = "Fishery Selectivity") 
)

  # Survey selectivity
  surv_sel_df <- melt(Surv_selex_at_age)
  names(surv_sel_df) <- c("Year", "Age", "Sim", "Selex")
  
  surv_sel_df <- surv_sel_df %>% 
    mutate(Age = parse_number(substr(Age, 5, 10)))

print(
  ggplot(surv_sel_df, aes(x = as.numeric(Year), y = Age,
                          fill = Selex)) +
    geom_tile() +
    facet_wrap(~Sim, ncol = 1, scales = "free") +
    scale_fill_viridis_c() +
    theme_bw() +
    theme(legend.position = "top") +
    labs(x = "Year",y = "Age", title = "Survey Selectivity") 
)


# Natural Mortality ---------------------------------------------------------------

  nat_mort <- melt(Mort_at_age)
  names(nat_mort) <- c("Year", "Age", "Sim", "Mortality")
  
  nat_mort <- nat_mort %>% 
    mutate(Age = parse_number(substr(Age, 5, 10)))

print( ggplot(nat_mort, aes(x = as.numeric(Year), y = Age,
                          fill = Mortality)) +
    geom_tile() +
    facet_wrap(~Sim, ncol = 1, scales = "free") +
    scale_fill_viridis_c() +
    theme_bw() +
    theme(legend.position = "top") +
    labs(x = "Year",y = "Age", title = "Natural Mortality") 
  )


# Catchability + Fishing Mortality------------------------------------------------------------


  ### Fishery -----------------------------------------------------------------

  q_df_fish <- melt(q_Fish)
  names(q_df_fish) <- c("Year", "Sim", "q")
  
print(
  ggplot(q_df_fish, aes(x = as.numeric(Year), y = q, color = Sim)) +
    geom_line(size = 1.3) +
    facet_wrap(~Sim, ncol = 1, scales = "free") +
    theme_bw() +
    theme(legend.position = "top") +
    labs(x = "Year",y = "q", title = "Fishery Catchability") 
)


  ### Survey ------------------------------------------------------------------
  
  q_df_surv <- melt(q_Surv)
  names(q_df_surv) <- c("Year", "Sim", "q")
  
  print(
    ggplot(q_df_surv, aes(x = as.numeric(Year), y = q, color = Sim)) +
      geom_line(size = 1.3) +
      facet_wrap(~Sim, ncol = 1, scales = "free") +
      theme_bw() +
      theme(legend.position = "top") +
      labs(x = "Year",y = "q", title = "Survey Catchability") 
  )
  

  ### Fishing Mortality -------------------------------------------------------
  
  fmort_df <- melt(fish_mort)
  names(fmort_df) <- c("Year", "Sim", "F_Mort")
    
  print(
    ggplot(fmort_df, aes(x = as.numeric(Year), y = F_Mort, color = Sim)) +
      geom_line(size = 1.3) +
      facet_wrap(~Sim, ncol = 1, scales = "free") +
      theme_bw() +
      theme(legend.position = "top") +
      labs(x = "Year",y = "Fishing Mortality", title = "Fishing Mortality") 
  )

  dev.off()
  
} # end function