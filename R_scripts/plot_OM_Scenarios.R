# Purpose: To plot and illustrate OM scenarios
# Creator: Matthew LH. Cheng (UAF-CFOS)
# Date 3/13/23

# Set Up ------------------------------------------------------------------

library(here)
library(tidyverse)
library(ggpubr)
library(geomtextpath)
library(patchwork)

# Get path to OMs
om_scenario_path <- here("output", "OM_Scenarios")

# Selectivity Differences -------------------------------------------------
bins <- 1:30
### Logistic-Logistic ----------------------------------------------------------------

# Get a50 value males
a50m1 <- 5.5
a50m2 <- 10.5
deltam1 <- 0.75
deltam2 <- 1.25

# Get a50 value females
a50f1 <- 3.5
a50f2 <- 8.5
deltaf1 <- 0.65
deltaf2 <- 1.15

# Compute selex 
selexm1 <- cbind(1 / (1 + exp(-1 * ((bins - a50m1)/deltam1) )), "Males", "Fleet 1", age = bins)
selexm2 <- cbind(1 / (1 + exp(-1 * ((bins - a50m2)/deltam2) )) , "Males", "Fleet 2", age = bins)
selexf1 <- cbind(1 / (1 + exp(-1 * ((bins - a50f1)/deltaf1) )) , "Females", "Fleet 1", age = bins)
selexf2 <- cbind(1 / (1 + exp(-1 * ((bins - a50f2)/deltaf2) )) , "Females", "Fleet 2", age = bins)

# Bind these together
selex_logist <- data.frame(rbind(selexm1, selexm2, selexf1, selexf2), Type = "Logist_Logist")
colnames(selex_logist) <- c("Selex", "Sex", "Fleet", "Age", "Type")

### Logistic-Gamma (Old) ----------------------------------------------------------

# Gamma parameters
amaxm2 <- 19
deltam2 <- 8
amaxf2 <- 15.5
deltaf2 <- 8

# Get Selex here
pm2 <- (0.5 * (sqrt(amaxm2^2 + 4*deltam2^2) - amaxm2))
selex_m2 <- cbind((bins/amaxm2) ^ (amaxm2/pm2) * exp((amaxm2 - bins) / pm2), "Males", "Fleet 2", age = bins)
pf2 <- (0.5 * (sqrt(amaxf2^2 + 4*deltaf2^2) - amaxf2))
selex_f2 <- cbind((bins/amaxf2) ^ (amaxf2/pf2) * exp((amaxf2 - bins) / pf2), "Females", "Fleet 2", age = bins)

# Bind these together
selex_gamma_old <- data.frame(rbind(selexm1, selex_m2, selexf1, selex_f2), Type = "Logist-Gamma-Old")
colnames(selex_gamma_old) <- c("Selex", "Sex", "Fleet", "Age", "Type")


# Logistic-Gamma (Young) --------------------------------------------------

amaxm2 <- 7
deltam2 <- 6.5
amaxf2 <- 5
deltaf2 <- 5

# Get Selex here
pm2 <- (0.5 * (sqrt(amaxm2^2 + 4*deltam2^2) - amaxm2))
selex_m2 <- cbind((bins/amaxm2) ^ (amaxm2/pm2) * exp((amaxm2 - bins) / pm2), "Males", "Fleet 2", age = bins)
pf2 <- (0.5 * (sqrt(amaxf2^2 + 4*deltaf2^2) - amaxf2))
selex_f2 <- cbind((bins/amaxf2) ^ (amaxf2/pf2) * exp((amaxf2 - bins) / pf2), "Females", "Fleet 2", age = bins)

# Bind these together
selex_gamma_young <- data.frame(rbind(selexm1, selex_m2, selexf1, selex_f2), Type = "Logist-Gamma-Young")
colnames(selex_gamma_young) <- c("Selex", "Sex", "Fleet", "Age", "Type")


# Plot Selectivity Scenarios -------------------------------------------------------------------

# Coerce into dataframe
selex_all <- rbind(selex_logist, selex_gamma_old, selex_gamma_young)

pdf(here("figs", "OM_Scenarios", "selex_scenario.pdf"), width = 15, height = 10)
(selex_plot <- ggplot(selex_all %>% 
                        mutate(Type = factor(Type,
                                             levels = c("Logist_Logist",
                                                        "Logist-Gamma-Old",
                                                        "Logist-Gamma-Young"))), aes(x = as.numeric(Age), 
                                     y = as.numeric(paste(Selex)), 
                                     lty = factor(Fleet))) +
  geom_line(size = 1.5, alpha = 1) + 
  facet_grid(Type~Sex) +
  theme_bw() +
  labs(x = "Ages", y = "Selectivity", lty = "Fleet") +
  theme(legend.position = "none",
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 15)))
dev.off()

# Presentation Figure
pdf(here("figs", "Presentation_Figures", "selex_scenario.pdf"), width = 15, height = 5)
(selex_plot_1 <- ggplot(selex_all %>% filter(Sex == "Females") %>% 
                        mutate(Type = factor(Type,
                                             levels = c("Logist_Logist",
                                                        "Logist-Gamma-Old",
                                                        "Logist-Gamma-Young"))), 
                      aes(x = as.numeric(Age),  y = as.numeric(paste(Selex)), lty = factor(Fleet))) +
    geom_line(size = 1.5, alpha = 1) + 
    facet_grid(Sex~Type) +
    theme_bw() +
    labs(x = "Ages", y = "Selectivity", lty = "Fleet") +
    theme(legend.position = "none",
          axis.title = element_text(size = 17),
          axis.text = element_text(size = 15, color = "black"),
          legend.title = element_text(size = 17),
          legend.text = element_text(size = 15),
          strip.text = element_text(size = 15)))
dev.off()

# Get Quantities from OMs -------------------------------------------------------
# Get OM files
OMs <- list.files(here("output", "OM_Scenarios"))

# empty containers
harv_rates_df <- data.frame() # Harvest Rate dataframe
neff_df <- data.frame() # Effective sample size dataframe
catch_df <- data.frame() # Catch dataframe
fmort_df <- data.frame() # Fishing mortality dataframe
ssb_df <- data.frame() # spawning stock biomass
naa_df <- data.frame() # spawning stock biomass

for(i in 1:length(OMs)) {
  
  # Load in OMs
  load(here(om_scenario_path, OMs[i], paste(OMs[i], ".RData", sep = "")))
  
  # Get effective sample sizes
  nf_df <- reshape2::melt(oms$Input_N_Fish)
  colnames(nf_df) <- c("Year", "Fleet", "Value")
  
  # Melt to df and rename - get harvest rate
  hr_df <- reshape2::melt(oms$Harvest_Rate)
  colnames(hr_df) <- c("Year", "Fleet", "Sim", "Value")
  
  # Get fishing mortality rate
  fmort_om <- reshape2::melt(oms$fish_mort)
  colnames(fmort_om) <- c("Year", "Fleet", "Sim", "Value")
  # Parse out numbers
  fmort_om <- fmort_om %>% 
    mutate(Year = parse_number(as.character(Year)),
           Fleet = parse_number(as.character(Fleet)),
           Sim = parse_number(as.character(Sim)))
  
  # Get Catch
  cat_df <- reshape2::melt(oms$Catch_agg)
  colnames(cat_df) <- c("Year", "Fleet", "Sim", "Value")
  
  # Get SSB
  ssb_om <- reshape2::melt(oms$SSB)
  colnames(ssb_om) <- c("Year", "Sim", "Value")
  
  # Parse out numbers
  ssb_om <- ssb_om %>% 
    mutate(Year = parse_number(as.character(Year)),
           Sim = parse_number(as.character(Sim)))
  
  # Get Numbers at age
  naa_om <- reshape2::melt(oms$N_at_age)
  colnames(naa_om) <- c("Year", "Age", "Sex", "Sim", "Value")
  
  # parse out numbers
  naa_om <- naa_om %>% 
    mutate(Year = parse_number(as.character(Year)),
           Age = parse_number(as.character(Age)),
           Sex = parse_number(as.character(Sex)),
           Sim = parse_number(as.character(Sim)))
  
  hr_df$OM_Scenario <- OMs[i]
  nf_df$OM_Scenario <- OMs[i]
  cat_df$OM_Scenario <- OMs[i]
  fmort_om$OM_Scenario <- OMs[i]
  ssb_om$OM_Scenario <- OMs[i]
  naa_om$OM_Scenario <- OMs[i]
  
  # Bind together
  harv_rates_df <- rbind(harv_rates_df, hr_df)
  neff_df <- rbind(nf_df, neff_df)
  catch_df <- rbind(cat_df, catch_df)
  fmort_df <- rbind(fmort_om, fmort_df)
  ssb_df <- rbind(ssb_df, ssb_om)
  naa_df <- rbind(naa_df, naa_om)
}


# Fishing Mortality Plot --------------------------------------------------

fmort_plot_df <- fmort_df %>% 
  filter(Sim == 1,
         str_detect(OM_Scenario, "High"),
         !str_detect(OM_Scenario, "_O|_LL")) %>% 
mutate(
  OM_Scenario = str_remove(OM_Scenario, "_High"),
  OM_Scenario = str_remove(OM_Scenario, "_Y"),
  OM_Scenario = str_remove(OM_Scenario, "_LL|_LG"),
  BreakPoint = case_when(
    str_detect(OM_Scenario, "Fast") & Year == 28 ~ 28,
    str_detect(OM_Scenario, "Fast") & Year == 30 ~ 30,
    str_detect(OM_Scenario, "Fast") & Year == 50 ~ 50,
    str_detect(OM_Scenario, "Slow") & Year == 40 ~ 40,
    str_detect(OM_Scenario, "Slow") & Year == 50 ~ 50,
    str_detect(OM_Scenario, "Slow") & Year == 70 ~ 70),
  # Names for breakpoints
  BreakPoint_Name = case_when(
    str_detect(OM_Scenario, "Fast") & Year == 28 ~ "Fleet Intersect",
    str_detect(OM_Scenario, "Fast") & Year == 30 ~ "Fleet Trans End",
    str_detect(OM_Scenario, "Fast") & Year == 50 ~ "Terminal",
    str_detect(OM_Scenario, "Slow") & Year == 40 ~ "Fleet Intersect",
    str_detect(OM_Scenario, "Slow") & Year == 50 ~ "Fleet Trans End",
    str_detect(OM_Scenario, "Slow") & Year == 70 ~ "Terminal"))

pdf(here("figs", "OM_Scenarios", "F_Scenario.pdf"), width = 20, height = 10)

(fmort_plot <- ggplot(fmort_plot_df, aes(x = Year, y = Value,  lty = factor(Fleet))) +
  geom_line(size = 0.9) +
  coord_cartesian(ylim = c(0, 0.13)) + 
  geom_segment(aes(x = BreakPoint, xend = BreakPoint,
                   y=-Inf, yend=0.115, color = BreakPoint_Name), lty = 1, lwd = 1.3) +
  scale_color_manual(values = viridis::viridis(n = 50)[c(1, 20, 43)],
                     na.translate = F) +
  facet_wrap(~OM_Scenario, ncol = 2, scales = "free_x") +
  theme_bw() +
  labs(x = "Year", y =  "Fishing Mortality Multiplier", linetype = "Fleet",
       color = "Assessment Period") +
  theme(legend.position = "top",
        axis.title = element_text(size = 17),
        legend.key.size = unit(1.5, "cm"),
        axis.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 17),
        strip.text = element_text(size = 17)))

dev.off()

# Presentation Figure
pdf(here("figs", "Presentation_Figures", "F_Scenario.pdf"), width = 8, height = 10)
(fmort_plot_1 <- ggplot(fmort_plot_df, aes(x = Year, y = Value,  lty = factor(Fleet))) +
    geom_line(size = 0.9) +
    coord_cartesian(ylim = c(0, 0.13)) + 
    facet_wrap(~OM_Scenario, nrow = 2, scales = "free_x") +
    theme_bw() +
    labs(x = "Year", y =  "Fishing Mortality Multiplier", linetype = "Fleet") +
    theme(legend.position = "left",
          axis.title = element_text(size = 17),
          legend.key.size = unit(1.5, "cm"),
          axis.text = element_text(size = 15, color = "black"),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 17),
          strip.text = element_text(size = 17)))
dev.off()

# Harvest Rate Plot -------------------------------------------------------

# Residual munging for numbers, etc
harv_rates_df_sum <- harv_rates_df %>% 
  drop_na() %>% 
  filter(str_detect(OM_Scenario, "High")) %>% 
  mutate(Year = parse_number(as.character(Year)),
         Fleet = parse_number(as.character(Fleet)),
         OM_Scenario = str_remove(OM_Scenario, "_High")) %>% 
  group_by(Year, OM_Scenario, Fleet) %>% 
  summarize(Median_HR = median(Value),
            Lwr_95_HR = quantile(Value, 0.025),
            Upr_95_HR = quantile(Value, 0.975)) %>% 
  # Renaming for ordering purposes
  mutate(# When breakpoints occur in doing assessment
         BreakPoint = case_when(
           str_detect(OM_Scenario, "Fast") & Year == 28 ~ 28,
           str_detect(OM_Scenario, "Fast") & Year == 30 ~ 30,
           str_detect(OM_Scenario, "Fast") & Year == 50 ~ 50,
           str_detect(OM_Scenario, "Slow") & Year == 40 ~ 40,
           str_detect(OM_Scenario, "Slow") & Year == 50 ~ 50,
           str_detect(OM_Scenario, "Slow") & Year == 70 ~ 70),
         
         # Names for breakpoints
         BreakPoint_Name = case_when(
           str_detect(OM_Scenario, "Fast") & Year == 28 ~ "Int",
           str_detect(OM_Scenario, "Fast") & Year == 30 ~ "TrxE",
           str_detect(OM_Scenario, "Fast") & Year == 50 ~ "Term",
           str_detect(OM_Scenario, "Slow") & Year == 40 ~ "Int",
           str_detect(OM_Scenario, "Slow") & Year == 50 ~ "TrxE",
           str_detect(OM_Scenario, "Slow") & Year == 70 ~ "Term"))

pdf(here("figs", "OM_Scenarios", "HarvestRate_Scenarios.pdf"), width = 28, height = 13)

(harv_rate_plot <- ggplot(harv_rates_df_sum, aes(x = Year, y = Median_HR,
                          color = factor(Fleet), fill = factor(Fleet), 
                          ymin = Lwr_95_HR, ymax = Upr_95_HR)) +
  geom_ribbon(alpha = 0.35) +
  geom_line(size = 1.5) +
  geom_text(aes(x = BreakPoint, y = 0.126, label = BreakPoint_Name), size = 7,
            color = "black", angle = 90) +
  geom_segment(aes(x = BreakPoint, xend= BreakPoint,y=-Inf, yend=0.115), color="black",
               lty = 2) + 
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  coord_cartesian(ylim = c(0, 0.15)) +
  facet_wrap(~OM_Scenario, scales = "free_x") +
  labs(x = "Year", y = "Harvest Rate", fill = "Fleet",
       color = "Fleet", linetype = "Fleet") +
  theme_bw() +
  theme(legend.position = "top",
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 11)) )
dev.off()
  
# Changes in catch --------------------------------------------------------
# Residual munging for numbers, etc
catch_df_sum <- catch_df %>% 
  drop_na() %>% 
  filter(str_detect(OM_Scenario, "_High")) %>% 
  mutate(Year = parse_number(as.character(Year)),
         Fleet = parse_number(as.character(Fleet)),
         OM_Scenario = str_remove(OM_Scenario, "_High"),
         Speed = ifelse(str_detect(OM_Scenario, "Fast"), "Fast", "Slow"),
         Selex = ifelse(str_detect(OM_Scenario, "LG_O"), "Logist-Gamma-Old",
                 ifelse(str_detect(OM_Scenario, "LG_Y"), 'Logist-Gamma-Young', "Logist-Logist")),
         Selex = factor(Selex, levels = c("Logist-Logist",
                                          "Logist-Gamma-Old",
                                          "Logist-Gamma-Young"))
         ) %>% 
  group_by(Year, OM_Scenario, Fleet, Speed) %>% 
  summarize(Median_HR = median(Value),
            Lwr_95_HR = quantile(Value, 0.025),
            Upr_95_HR = quantile(Value, 0.975)) 

pdf(here("figs", "OM_Scenarios", "Catch_Scenarios.pdf"), width = 28, height = 10)
(catch_plot <- ggplot(catch_df_sum, aes(x = Year, y = Median_HR,  lty = factor(Fleet),
                      ymin = Lwr_95_HR, ymax = Upr_95_HR)) +
    geom_ribbon(alpha = 0.35, color = NA) +
    geom_line(size = 1.5) +
    facet_grid(OM_Scenario~Speed, scales = "free_x") +
    labs(x = "Year", y = "Catch", fill = "Fleet",
         color = "Fleet", linetype = "Fleet") +
    theme_bw() +
    theme(legend.position = "none",
          axis.title = element_text(size = 17),
          axis.text = element_text(size = 15, color = "black"),
          strip.text = element_text(size = 17)) )
dev.off()

# Aggregated catch
agg_catch = catch_df %>% 
  drop_na() %>% 
  filter(str_detect(OM_Scenario, "_High")) %>% 
  mutate(Year = parse_number(as.character(Year)),
         Fleet = parse_number(as.character(Fleet)),
         OM_Scenario = str_remove(OM_Scenario, "_High"),
         Speed = ifelse(str_detect(OM_Scenario, "Fast"), "Fast", "Slow"),
         Selex = ifelse(str_detect(OM_Scenario, "LG_O"), "Logist-Gamma-Old", 
                        ifelse(str_detect(OM_Scenario, "LG_Y"), 'Logist-Gamma-Young', "Logist-Logist")),
         Selex = factor(Selex, levels = c("Logist-Logist", 
                                          "Logist-Gamma-Old",
                                          "Logist-Gamma-Young"))) %>% 
  group_by(Year, OM_Scenario, Speed, Selex, Sim) %>% 
  mutate(Sum_Catch = sum(Value)) %>% 
  ungroup() %>% 
  group_by(Year, OM_Scenario, Speed, Selex) %>% 
  summarize(
    Median_HR = median(Sum_Catch),
    Lwr_95_HR = quantile(Sum_Catch, 0.025),
    Upr_95_HR = quantile(Sum_Catch, 0.975))

(ggplot(agg_catch, aes(x = Year, y = Median_HR,  
                                        ymin = Lwr_95_HR, ymax = Upr_95_HR,
                                        fill = Speed, color = Speed)) +
    geom_ribbon(alpha = 0.35, color = NA) +
    geom_line(size = 1.5) +
    facet_grid(Selex~., scales = "free_x") +
    labs(x = "Year", y = "Catch", fill = "Fleet",
         color = "Fleet", linetype = "Fleet") +
    theme_bw() +
    theme(legend.position = "top",
          axis.title = element_text(size = 17),
          axis.text = element_text(size = 15, color = "black"),
          legend.title = element_text(size = 17),
          legend.text = element_text(size = 15),
          strip.text = element_text(size = 11)) )

# Catch Ratio
ratio_df <- catch_df %>% 
  drop_na() %>% 
  filter(str_detect(OM_Scenario, "High")) %>% 
  pivot_wider(names_from = "Fleet", values_from = Value) %>% 
  rowwise() %>% 
  mutate(`Fleet 1` = Fish_Fleet_1 / (Fish_Fleet_1 + Fish_Fleet_2),
         `Fleet 2` = Fish_Fleet_2 / (Fish_Fleet_1 + Fish_Fleet_2)) %>% 
  dplyr::select(-Fish_Fleet_1, -Fish_Fleet_2) %>% 
  pivot_longer(names_to = "Fleet", values_to = "Value", cols = c("Fleet 1",
                                                                 "Fleet 2")) %>% 
  group_by(Year, OM_Scenario, Fleet) %>% 
  summarize(Median= median(Value),
            Lwr_95= quantile(Value, 0.025),
            Upr_95= quantile(Value, 0.975)) %>% 
  mutate(Year = parse_number(as.character(Year)),
         OM_Scenario = str_remove(OM_Scenario, "_High")) 

pdf(here("figs", "OM_Scenarios", "Catch_Ratio.pdf"), width = 28, height = 10)
# Plot out the catch ratio now!
print(
  ggplot(ratio_df, aes(x = Year, y = Median, 
                       ymin = Lwr_95, ymax = Upr_95,
                       fill = factor(Fleet), color = factor(Fleet))) +
    geom_ribbon(alpha = 0.35) +
    geom_line(size = 1.5) +
    scale_color_brewer(palette = "Set2") + 
    scale_fill_brewer(palette = "Set2") + 
    labs(x = "Year", y = "Catch Ratio", color = "Fleet",
         lty = "Fleet", fill = "Fleet") +
    facet_wrap(~OM_Scenario, scales = "free_x") +
    theme_bw() +
    theme(legend.position = "top",
          axis.title = element_text(size = 17),
          axis.text = element_text(size = 15, color = "black"),
          legend.title = element_text(size = 17),
          legend.text = element_text(size = 15),
          strip.text = element_text(size = 11)) 
)
dev.off()

# Get SSB Trajectories ----------------------------------------------------

# Munging for plot dataframe
ssb_plot_df <- ssb_df %>% 
  drop_na() %>% 
  filter(str_detect(OM_Scenario, "High")) %>% 
  mutate(OM_Scenario = str_remove(OM_Scenario, "_High"),
         Speed = ifelse(str_detect(OM_Scenario, "Fast"), "Fast", "Slow"),
         Sel_Type = ifelse(str_detect(OM_Scenario, "LG_O"), "Logist-Gamma-Old", 
                           ifelse(str_detect(OM_Scenario, "LG_Y"),
                                  'Logist-Gamma-Young', "Logist-Logist")),
         Sel_Type = factor(Sel_Type, levels = c("Logist-Logist", 
                                                "Logist-Gamma-Old",
                                                "Logist-Gamma-Young"))) %>% 
  group_by(OM_Scenario, Year, Sel_Type, Speed) %>% 
  summarize(median = median(Value),
         Lwr_95 = quantile(Value, 0.025),
         Upr_95 = quantile(Value, 0.975))

pdf(here("figs", "OM_Scenarios", "SBB_Trajectory.pdf"), width = 18, height = 10)
(ssb_plot <- ssb_plot_df %>% 
  ggplot(aes(x = Year, y = median, ymin = Lwr_95, ymax = Upr_95)) +
  geom_line(lwd = 1.5) +
  geom_ribbon(alpha = 0.35) +
  theme_bw() +
  facet_grid(Sel_Type~Speed, scales = "free_x") +
  labs(x = 'Year', y = "Spawning Stock Biomass",
     color = "OM Scenario") +
  theme(legend.position = "top",
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 15)) )
dev.off()

# Numbers at Age ----------------------------------------------------------

med_naa_df <- naa_df %>% 
  filter(str_detect(OM_Scenario, "High")) %>% 
  mutate(Sex = ifelse(Sex == 1, "Female", "Male"),
         OM_Scenario = str_remove(OM_Scenario, "_High")
         # OM_Scenario = factor(OM_Scenario,
         #                      levels = c("Fast_LG_O", "Fast_LG_Y",
         #                                 "Slow_LG_O", "Slow_LG_Y",
         #                                 "Fast_LL", "Slow_LL"))
         ) %>% 
  group_by(Age, Year, Sex, OM_Scenario) %>% 
  summarize(median = median(Value))

pdf(here("figs", "OM_Scenarios", "NAA_plot.pdf"), width = 15, height = 10)

(naa_plot = ggplot(med_naa_df %>% filter(Age %in% c(seq(2,30, 3), 30),
                                         # str_detect(OM_Scenario, "LL"),
                                         Sex == "Female") %>% 
                     mutate(OM_Scenario = factor(OM_Scenario,
                                                 levels = c(
                                                   "Fast_LL",
                                                   "Fast_LL_Rev",
                                                   "Fast_LG_O",
                                                   "Fast_LG_O_Rev",
                                                   "Fast_LG_Y",
                                                   "Fast_LG_Y_Rev",
                                                   "Slow_LL",
                                                   "Slow_LL_Rev",
                                                   "Slow_LG_O",
                                                   "Slow_LG_O_Rev",
                                                   "Slow_LG_Y",
                                                   "Slow_LG_Y_Rev"
                                                 ))), 
       aes(x = Year, y = median, fill = factor(Age))) +
  geom_col(color = "black", alpha = 0.85) +
  facet_wrap(~OM_Scenario, scales = "free") +
  labs(x = "Year", y = "Median Numbers-at-age", fill = "Age") +
  scale_fill_viridis_d() +
  theme_bw() +
  theme(legend.position = "top",
        axis.title = element_text(size = 25),
        axis.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 25),
        legend.text = element_text(size = 23),
        strip.text = element_text(size = 18)) )

dev.off()



# Effective Sample Size ---------------------------------------------------

pdf(here("figs", "OM_Scenarios", "Neff_plot.pdf"), width = 15, height = 10)

print(
  neff_df %>% 
    ggplot(aes(x = Year, y = Value, color = factor(Fleet))) +
    geom_line(size = 1.5) +
    facet_wrap(~OM_Scenario) +
    scale_color_brewer(palette = "Set2") + 
    labs(x = "Year", y = "Effective Sample Size", color = "Fleet") +
    theme_bw() +
    theme(legend.position = "top",
          axis.title = element_text(size = 17),
          axis.text = element_text(size = 15, color = "black"),
          legend.title = element_text(size = 17),
          legend.text = element_text(size = 15),
          strip.text = element_text(size = 11)) 
)
dev.off()


# Combine Plots -----------------------------------------------------------


pdf(here("figs", "OM_Scenarios", "Fig1_OM_Scenarios.pdf"), width = 17, height = 15)
# Combine plots
comb_plots = ggpubr::ggarrange(selex_plot, catch_plot,
                  ssb_plot, ncol = 3, 
                  labels = c("B", "C", "D"),
                  common.legend = FALSE,  label.x = -0.005,
                  font.label=list(color="black",size = 25))
# Now get all plots
ggpubr::ggarrange(fmort_plot, comb_plots, 
                  labels = "A",
                  label.x = 0.0, 
                  label.y = 0.95,
                  font.label=list(color="black",size = 25),
                  ncol = 1, heights = c(0.6, 1))
dev.off()


