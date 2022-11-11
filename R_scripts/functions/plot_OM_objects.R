# Purpose: Plotting outputs of operating model
# Date: 10/29/22
# Creator: Matthew LH. Cheng

#' @param path Path to output folder
#' @param file_name Name of pdf output

plot_OM <- function(path, file_name) {
  
  pdf(here(path, file_name), height = 7.5, width = 13)
  
  require(tidyverse)

# Recruitment -------------------------------------------------------------

  # Plot Recruitment scenarios
  rec_df <- melt(rec_total)
  names(rec_df) <- c("Year", "Sim", "Rec")
  
  # Total recruitment across time and simulations
print(  ggplot(rec_df, aes(x = as.numeric(Year), y = Rec, group = Sim)) +
          geom_line(size = 1.5, color = "grey", alpha = 0.75) +
          theme_bw() +
          theme(legend.position = "none") +
          labs(x = "Year",y = "Recruitment Total") )
  
# SSB ---------------------------------------------------------------------

  # Plot SSB
  ssb_df <- melt(SSB)
  names(ssb_df) <- c("Year", "Sim", "SSB")
  
print(
  ggplot(ssb_df, aes(x = as.numeric(Year), y = SSB)) +
    geom_line(size = 1.5, color = "grey", aes(group = Sim), alpha = 0.75) +
    geom_smooth(se = FALSE) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "Year",y = "SSB") 
)
  

# Recruitment ~ SSB -------------------------------------------------------

  spr <- data.frame(SSB = ssb_df$SSB, Rec = rec_df$Rec, Sim = rec_df$Sim)

print(
  
  ggplot(spr, aes(x = SSB, y = Rec, group = Sim)) +
    geom_point(alpha = 0.3) +
    geom_smooth(aes(color = Sim), se = FALSE) +
    # facet_wrap(~Sim, ncol = 1, scales = "free") +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "SSB",y = "Recruitment Total") 
)
  

# Biomass and Numbers at Age ----------------------------------------------
  
  # biomass at age
  biom_df <- melt(Biom_at_age)
  names(biom_df) <- c("Year", "Age", "Sex", "Sim", "Biomass")
  
print(
  # Biomass at age across sexes
  ggplot(biom_df %>% 
           filter(Sim %in% c("Sim 1", "Sim 2", "Sim 3")), 
         aes(x = as.numeric(Year), y = Biomass, fill = Age, group = Sim)) +
    geom_col(position = "stack") +
    facet_wrap(~Sim, ncol = 3, scales = "free") +
    theme_bw() +
    theme(legend.position = "top") +
    labs(x = "Year",y = "Biomass at age")
)


print(
  # Biomass at age differntiated by sex
  ggplot(biom_df %>% 
           filter(Sim %in% c("Sim 1", "Sim 2", "Sim 3")), 
         aes(x = as.numeric(Year), y = Biomass, fill = Age, group = Sim)) +
    geom_col(position = "stack") +
    facet_grid(Sex~Sim, scales = "free") +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "Year",y = "Biomass at age")
)
  
print(
  # Biomass aggregated
  biom_df %>% 
    group_by(Year, Sim) %>% 
    summarize(Biomass = sum(Biomass, na.rm = TRUE)) %>% 
    ggplot(aes(x = as.numeric(Year), y = Biomass,  group = Sim)) +
    geom_line(size = 1.5, color = "grey", alpha = 0.75) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "Year",y = "Biomass")
)
  
  # Plot Numbers at age
  natage <- melt(N_at_age)
  names(natage) <- c("Year", "Age", "Sex", "Sim", "Numbers")
  
print(
  # Numbers across ages - aggregated sex
  ggplot(natage %>% 
           filter(Sim %in% c("Sim 1",  "Sim 2", "Sim 3")), 
         aes(x = as.numeric(Year), y = Numbers, fill = Age, group = Sim)) +
    geom_col(position = "stack") +
    facet_grid(Sex~Sim,  scales = "free") +
    theme_bw() +
    theme(legend.position = "top") +
    labs(x = "Year",y = "Numbers at age")
)
  
print(
  # Numbers aggregated by ages
  natage %>% 
    group_by(Year, Sim) %>% 
    summarize(Numbers = sum(Numbers, na.rm = TRUE)) %>% 
    ggplot(aes(x = as.numeric(Year), y = Numbers, group = Sim)) +
    geom_line(size = 1.5, color = "grey", alpha = 0.75) +
    geom_smooth(group = 1, se = F, alpha = 0.5)+
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "Year",y = "Numbers")
)

print(
  # Numbers aggregated by ages differenitaed by sex
  natage %>% 
    group_by(Year, Sim, Sex) %>% 
    summarize(Numbers = sum(Numbers, na.rm = TRUE)) %>% 
    ggplot(aes(x = as.numeric(Year), y = Numbers, group = Sim)) +
    geom_line(size = 1.5, color = "grey", alpha = 0.75) +
    geom_smooth(group = 1, se = F, alpha = 0.5)+
    facet_wrap(~Sex) +
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
    ggplot(aes(x = as.numeric(Age), y = prop_age, group = Sim, group = Sim )) +
    geom_line(size = 1.5, color = "grey", alpha = 0.75) +
    # facet_wrap(~Sim, ncol = 1, scales = "free") +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "Age",y = "Average proportion at age")
)

  

# Catch at age ------------------------------------------------------------

  catch_df <- melt(Catch_at_age)
  names(catch_df) <- c("Year", "Age", "Fleet", "Sex", "Sim", "Catch")
  
print(
  catch_df %>% 
    group_by(Year, Fleet, Sex, Sim) %>% # Aggregated catch
    summarize(Catch  = sum(Catch, na.rm = T)) %>% 
    ggplot(aes(x = as.numeric(Year), y = Catch, group = Sim)) +
    geom_line(size = 1.5, color = "grey", alpha = 0.75) +
    facet_grid(Sex~Fleet, scales = "free") +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "Age",y = "Catch")
)
  

# Composition data --------------------------------------------------------

tryCatch({ 
  
  # Fishery Comps
  comps_fish <- melt(Fish_Age_Comps)
  names(comps_fish) <- c("Year", "Age", "Fleet", "Sex","Sim", "Count")
  
print(
  ggplot(comps_fish %>% 
           filter(Sim %in% c("Sim 1")),
         aes(x = as.numeric(Year), y = Count, fill = Age, group = Sim)) +
    geom_col(position = "stack") +
    facet_grid(Sex~Fleet, scales = "free") +
    theme_bw() +
    theme(legend.position = "top") +
    labs(x = "Year",y = "Count", title = "Fishery Comps")
)
  
  comps_surv <- melt(Survey_Age_Comps)
  names(comps_surv) <- c("Year", "Age", "Fleet", "Sex","Sim", "Count")
  
print(
  ggplot(comps_surv %>% 
           filter(Sim %in% c("Sim 1")), 
         aes(x = as.numeric(Year), y = Count, fill = Age, group = Sim)) +
    geom_col(position = "stack") +
    facet_grid(Sex~Fleet, scales = "free") +
    theme_bw() +
    theme(legend.position = "top") +
    labs(x = "Year",y = "Count", title = "Survey Comps") 
)
  

# Indices of abundance ----------------------------------------------------

  # Fishery index
  idx_fish <- melt(Fishery_Index_Agg) 
  names(idx_fish) <- c("Year", "Fleet","Sim", "Index")
  
print(
  # Across ages
  ggplot(idx_fish, aes(as.numeric(Year), y = Index, group = Sim))+
    geom_line(size = 1.1, color = "grey", alpha = 0.75) +
    facet_wrap(~Fleet, scales = "free_x") +
    theme_bw() +
    theme(legend.position = "top") +
    labs(x = "Year",y = "Biomass index", title = "Fishery Index") 
)
  
  # Survey index
  idx_surv <- melt(Survey_Index_Agg)
  names(idx_surv) <- c("Year", "Fleet","Sim", "Index")
  
print(
  # Across ages
  ggplot(idx_surv, aes(as.numeric(Year), y = Index, group = Sim))+
    geom_line(size = 1.1, color = "grey", alpha = 0.75) +
    facet_wrap(~Fleet, scales = "free_x") +
    theme_bw() +
    theme(legend.position = "top") +
    labs(x = "Year",y = "Biomass index at age", title = "Survey Index") 
)

  }, error = function(error) {cat("ERROR :",conditionMessage(error), "\n")}) # end try catch statement
  


# Parameterizations -------------------------------------------------------


  ### Selectivity -------------------------------------------------------------

  # Fishery selectivity
  fish_sel_df <- melt(Fish_selex_at_age)
  names(fish_sel_df) <- c("Year", "Age", "Fleet", "Sex", "Sim", "Selex")
  
  fish_sel_df <- fish_sel_df %>% 
    mutate(Age = parse_number(substr(Age, 5, 10)))
  
print( # print tiles here
  ggplot(fish_sel_df %>% filter(Sim %in% c("Sim 1")),
         aes(x = as.numeric(Year), y = Age, fill = Selex)) +
    geom_tile() +
    facet_grid(Sex~Fleet, scales = "free") +
    scale_fill_viridis_c() +
    theme_bw() +
    theme(legend.position = "top") +
    labs(x = "Year",y = "Age", title = "Fishery Selectivity") 
)

print( # Print lines here!
  ggplot(fish_sel_df, aes(x = Age, y = Selex, group = as.numeric(Year))) +
    geom_line(color = "grey", size = 1.5) +
    facet_grid(Sex~Fleet, scales = "free") +
    scale_fill_viridis_c() +
    theme_bw() +
    theme(legend.position = "top") +
    labs(x = "Age",y = "Prop of Selex", title = "Fishery Selectivity") 
)


  # Survey selectivity
  surv_sel_df <- melt(Surv_selex_at_age)
  names(surv_sel_df) <- c("Year", "Age", "Fleet", "Sex", "Sim", "Selex")
  
  surv_sel_df <- surv_sel_df %>% 
    mutate(Age = parse_number(substr(Age, 5, 10)))

print( # print tiles here
  ggplot(surv_sel_df %>% 
           filter(Sim %in% c("Sim 1")),
         aes(x = as.numeric(Year), y = Age, fill = Selex)) +
    geom_tile() +
    facet_grid(Sex~Fleet, scales = "free") +
    scale_fill_viridis_c() +
    theme_bw() +
    theme(legend.position = "top") +
    labs(x = "Year",y = "Age", title = "Survey Selectivity") 
)

print( # Print lines here!
  ggplot(surv_sel_df, aes(x = Age, y = Selex, group = as.numeric(Year))) +
    geom_line(color = "grey", size = 1.5) +
    facet_grid(Sex~Fleet, scales = "free") +
    scale_fill_viridis_c() +
    theme_bw() +
    theme(legend.position = "top") +
    labs(x = "Age",y = "Prop of Selex", title = "Survey Selectivity") 
)


# Natural Mortality ---------------------------------------------------------------

  nat_mort <- melt(Mort_at_age)
  names(nat_mort) <- c("Year", "Age", "Sim", "Mortality")
  
  nat_mort <- nat_mort %>% 
    mutate(Age = parse_number(substr(Age, 5, 10)))

print( ggplot(nat_mort %>% 
                filter(Sim %in% c("Sim 1")),
              aes(x = as.numeric(Year), y = Age, fill = Mortality, group = Sim)) +
    geom_tile() +
    facet_wrap(~Sim, ncol = 3, scales = "free") +
    scale_fill_viridis_c() +
    theme_bw() +
    theme(legend.position = "top") +
    labs(x = "Year",y = "Age", title = "Natural Mortality") 
  )


# Catchability + Fishing Mortality------------------------------------------------------------


  ### Fishery -----------------------------------------------------------------

  q_df_fish <- melt(q_Fish)
  names(q_df_fish) <- c("Year", "Fleet", "Sim", "q")
  
print(
  ggplot(q_df_fish %>% 
           filter(Sim %in% c("Sim 1")), 
         aes(x = as.numeric(Year), y = q)) +
    geom_line(size = 1.3) +
    facet_wrap(~Fleet, ncol = 3, scales = "free_y") +
    theme_bw() +
    theme(legend.position = "top") +
    labs(x = "Year",y = "q", title = "Fishery Catchability") 
)


  ### Survey ------------------------------------------------------------------
  
  q_df_surv <- melt(q_Surv)
  names(q_df_surv) <- c("Year", "Fleet", "Sim", "q")
  
  print(
    ggplot(q_df_surv %>% 
             filter(Sim %in% c("Sim 1")),
           aes(x = as.numeric(Year), y = q)) +
      geom_line(size = 1.3) +
      facet_wrap(~Fleet, ncol = 3, scales = "free_y") +
      theme_bw() +
      theme(legend.position = "top") +
      labs(x = "Year",y = "q", title = "Survey Catchability") 
  )
  

  ### Fishing Mortality -------------------------------------------------------
  
  fmort_df <- melt(fish_mort)
  names(fmort_df) <- c("Year", "Fleet", "Sim", "F_Mort") 

  print(
    ggplot(fmort_df %>% 
             filter(Sim %in% c("Sim 1")), 
           aes(x = as.numeric(Year), y = F_Mort)) +
      geom_line(size = 1.3) +
      facet_wrap(~Fleet, ncol = 3, scales = "free") +
      theme_bw() +
      theme(legend.position = "top") +
      labs(x = "Year",y = "Fishing Mortality", title = "Fishing Mortality") 
  ) +
    xlim(150, 200)
  
  dev.off()
  
} # end function
