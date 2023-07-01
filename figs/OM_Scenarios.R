# Purpose: To plot and illustrate OM scenarios
# Creator: Matthew LH. Cheng (UAF-CFOS)
# Date 3/13/23


# Set Up ------------------------------------------------------------------

library(here)
library(tidyverse)
library(ggpubr)
library(geomtextpath)

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
selex_logist <- data.frame(rbind(selexm1, selexm2, selexf1, selexf2), Type = "Logistic-Logistic (LL)")
colnames(selex_logist) <- c("Selex", "Sex", "Fleet", "Age", "Type")


### Logistic-Gamma ----------------------------------------------------------

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
selex_gamma <- data.frame(rbind(selexm1, selex_m2, selexf1, selex_f2), Type = "Logistic-Gamma (LG)")
colnames(selex_gamma) <- c("Selex", "Sex", "Fleet", "Age", "Type")

# Plot Selectivity Scenarios -------------------------------------------------------------------

# Coerce into dataframe
selex_all <- rbind(selex_logist, selex_gamma)

(selex_plot <- ggplot(selex_all, aes(x = as.numeric(Age), 
                                     y = as.numeric(paste(Selex)), lty = Fleet)) +
  geom_line(size = 1, alpha = 0.85) + 
  facet_grid(Sex~Type) +
  scale_linetype_manual(values = c(10,2)) +
  theme_bw() +
  labs(x = "Ages", y = "Selectivity") +
  theme(legend.position = "none",
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 17)))

# Get Quantities from OMs -------------------------------------------------------

# Load unique OMs
OMs <- c("Fast_LL_95_High", "Fast_LG_95_High", "Slow_LL_95_High", "Slow_LG_95_High",
         "Fast_LL_75_High", "Fast_LG_75_High", "Slow_LL_75_High", "Slow_LG_75_High")

# Harvest Rate dataframe
harv_rates_df <- data.frame()
# Effective sample size dataframe
neff_df <- data.frame()
# Catch dataframe
catch_df <- data.frame()
# Fishing mortality dataframe
fmort_df <- data.frame()
ssb_df <- data.frame() # spawning stock biomass

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
  
  
  hr_df$OM_Scenario <- OMs[i]
  nf_df$OM_Scenario <- OMs[i]
  cat_df$OM_Scenario <- OMs[i]
  fmort_om$OM_Scenario <- OMs[i]
  ssb_om$OM_Scenario <- OMs[i]
  
  # Bind together
  harv_rates_df <- rbind(harv_rates_df, hr_df)
  neff_df <- rbind(nf_df, neff_df)
  catch_df <- rbind(cat_df, catch_df)
  fmort_df <- rbind(fmort_om, fmort_df)
  ssb_df <- rbind(ssb_df, ssb_om)
}


# Fishing Mortality Plot --------------------------------------------------

fmort_plot_df <- fmort_df %>% 
  filter(Sim == 1, OM_Scenario %in% c("Fast_LG_75_High", "Slow_LG_75_High",
                            "Fast_LG_95_High", "Slow_LG_95_High")) %>% 
mutate(OM_Scenario = case_when(
  OM_Scenario == "Fast_LG_75_High" ~ "Fast (75% switch)",
  OM_Scenario == "Slow_LG_75_High" ~ "Slow (75% switch)",
  OM_Scenario == "Fast_LG_95_High" ~ "Fast (95% switch)",
  OM_Scenario == "Slow_LG_95_High" ~ "Slow (95% switch)"),
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

(fmort_plot <- ggplot(fmort_plot_df, aes(x = Year, y = Value, lty = factor(Fleet))) +
  geom_line(size = 0.9) +
  coord_cartesian(ylim = c(0, 0.13)) + 
  geom_segment(aes(x = BreakPoint, xend = BreakPoint,
                   y=-Inf, yend=0.115, color = BreakPoint_Name), lty = 1, lwd = 1.3) +
  scale_color_manual(values = viridis::viridis(n = 50)[c(1, 20, 43)],
                     na.translate = F) +
  scale_linetype_manual(values = c(10,2)) +
  facet_wrap(~OM_Scenario, scales = "free_x") +
  theme_bw() +
  labs(x = "Year", y = "Relative Fishing Mortality", linetype = "Fleet",
       color = "Time Component") +
  theme(legend.position = "top",
        axis.title = element_text(size = 17),
        legend.key.size = unit(1.5, "cm"),
        axis.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 17)))

# Harvest Rate Plot -------------------------------------------------------

# Residual munging for numbers, etc
harv_rates_df_sum <- harv_rates_df %>% 
  drop_na() %>% 
  mutate(Year = parse_number(as.character(Year)),
         Fleet = parse_number(as.character(Fleet))) %>% 
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

pdf(here("figs", "Hypothetical_OM", "HarvestRate_Scenarios.pdf"), width = 28, height = 8)

(harv_rate_plot <- ggplot(harv_rates_df_sum, aes(x = Year, y = Median_HR,
                          lty = factor(Fleet), ymin = Lwr_95_HR, ymax = Upr_95_HR)) +
  geom_ribbon(alpha = 0.35) +
  geom_line(size = 1.5) +
  geom_text(aes(x = BreakPoint, y = 0.126, label = BreakPoint_Name), size = 7,
            color = "black", angle = 90) +
  geom_segment(aes(x = BreakPoint, xend= BreakPoint,y=-Inf, yend=0.115), color="black",
               lty = 2) + 
  scale_linetype_manual(values = c(10, 2)) +
  coord_cartesian(ylim = c(0, 0.15)) +
  facet_wrap(~OM_Scenario, scales = "free_x", ncol = 4) +
  labs(x = "Year", y = "Harvest Rate", fill = "Fleet",
       color = "Fleet", linetype = "Fleet") +
  theme_bw() +
  theme(legend.position = "top",
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 17)) )

dev.off()
  
# Changes in catch --------------------------------------------------------
# Residual munging for numbers, etc
catch_df_sum <- catch_df %>% 
  drop_na() %>% 
  mutate(Year = parse_number(as.character(Year)),
         Fleet = parse_number(as.character(Fleet))) %>% 
  group_by(Year, OM_Scenario, Fleet) %>% 
  summarize(Median_HR = median(Value),
            Lwr_95_HR = quantile(Value, 0.025),
            Upr_95_HR = quantile(Value, 0.975)) 

pdf(here("figs", "Hypothetical_OM", "Catch_Scenarios.pdf"), width = 28, height = 8)

(catch_plot <- ggplot(catch_df_sum, aes(x = Year, y = Median_HR, 
                                        linetype = factor(Fleet),
                                        ymin = Lwr_95_HR, ymax = Upr_95_HR)) +
    geom_ribbon(alpha = 0.35) +
    geom_line(size = 1.5) +
    scale_linetype_manual(values = c(10, 2)) +
    facet_wrap(~OM_Scenario, scales = "free_x",  ncol = 4) +
    labs(x = "Year", y = "Catch", fill = "Fleet",
         color = "Fleet", linetype = "Fleet") +
    theme_bw() +
    theme(legend.position = "top",
          axis.title = element_text(size = 17),
          axis.text = element_text(size = 15, color = "black"),
          legend.title = element_text(size = 17),
          legend.text = element_text(size = 15),
          strip.text = element_text(size = 17)) )

dev.off()

# Catch Ratio
ratio_df <- catch_df %>% 
  drop_na() %>% 
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
  mutate(Year = parse_number(as.character(Year))) %>% 
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

pdf(here("figs", "Hypothetical_OM", "Catch_Ratio.pdf"), width = 28, height = 8)
# Plot out the catch ratio now!
print(
  ggplot(ratio_df, aes(x = Year, y = Median, 
                       ymin = Lwr_95, ymax = Upr_95,
                       lty = factor(Fleet), group = factor(Fleet))) +
    geom_ribbon(alpha = 0.35) +
    geom_line(size = 1.5) +
    geom_segment(aes(x = BreakPoint, xend = BreakPoint,
                     y=-Inf, yend= 0.9, color = BreakPoint_Name), lty = 1, lwd = 1.3) +
    scale_color_manual(values = viridis::viridis(n = 50)[c(1, 43, 20)],
                       na.translate = F) +
    scale_linetype_manual(values = c(10, 2)) + 
    labs(x = "Year", y = "Catch Ratio", color = "Time Component",
         lty = "Fleet") +
    facet_wrap(~OM_Scenario, scales = "free_x", ncol = 4) +
    theme_bw() +
    theme(legend.position = "top",
          axis.title = element_text(size = 17),
          axis.text = element_text(size = 15, color = "black"),
          legend.title = element_text(size = 17),
          legend.text = element_text(size = 15),
          strip.text = element_text(size = 17)) 
)
dev.off()

# Get SSB Trajectories ----------------------------------------------------

# Munging for plot dataframe
ssb_plot_df <- ssb_df %>% 
  filter(OM_Scenario %in% c("Fast_LG_75_High", "Slow_LG_75_High",
                            "Fast_LG_95_High", "Slow_LG_95_High")) %>% 
  mutate(OM_Scenario = case_when(
    OM_Scenario == "Fast_LG_75_High" ~ "Fast (75% switch)",
    OM_Scenario == "Slow_LG_75_High" ~ "Slow (75% switch)",
    OM_Scenario == "Fast_LG_95_High" ~ "Fast (95% switch)",
    OM_Scenario == "Slow_LG_95_High" ~ "Slow (95% switch)")) %>% 
  group_by(OM_Scenario, Year) %>% 
  mutate(median = median(Value))

(ssb_plot <- ssb_plot_df %>% 
  ggplot(aes(x = Year, y = Value, group = Sim)) +
  geom_line(alpha = 0.25) +
  geom_line(aes(y = median), color = "white", 
            lty = 2, lwd = 1.5) +
  coord_cartesian(ylim = c(0, 450)) +
  facet_wrap(~OM_Scenario, scales = "free_x") +
  theme_bw() +
  labs(x = 'Year', y = "Relative Spawning Stock Biomass") +
  theme(legend.position = "none",
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 17)) )

# Combine plots -----------------------------------------------------------

pdf(here("figs", "Hypothetical_OM", "OM_Scenarios.pdf"), width = 28, height = 8)
ggarrange(fmort_plot, ssb_plot, selex_plot, ncol = 3,
          common.legend = TRUE, legend="top",
          labels = c("A", "B", "C"),
          font.label = list(size = 25),
          label.x = 0.02, widths = c(0.5, 0.5, 0.6))
dev.off()             
   


