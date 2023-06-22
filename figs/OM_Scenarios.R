# Purpose: To plot and illustrate OM scenarios
# Creator: Matthew LH. Cheng (UAF-CFOS)
# Date 3/13/23


# Set Up ------------------------------------------------------------------

library(here)
library(tidyverse)
library(ggpubr)

# Get path to OMs
om_scenario_path <- here("output", "OM_Scenarios")

# Fleet Structure Shift ---------------------------------------------------
### Fast --------------------------------------------------------------------

chngpoint <- 24     # Change point for fmort
chngpoint_end <- 29 # when does change point end
mean_nat_mort <- 0.108 # Mean M
min_rel_F1_M <- 0.01 # Min F relative to M
max_rel_F2_M <- 0.85  # Max F relative to M
Start_F1 <- 1 * mean_nat_mort # Start F1
Start_F2 <- 0.01 * mean_nat_mort # Start F2
Fish_Start_yr <- 0
n_years <- 50

# Decrease Change
# F stays constant initially
constant_F1 <- seq(Start_F1, Start_F1, length.out = length(Fish_Start_yr:chngpoint))
# Change point for F (Ramp)
change_F1 <- seq(Start_F1, (min_rel_F1_M * mean_nat_mort), length.out = length(chngpoint:chngpoint_end))
# Change point end (plateaus/stays constant now)
change_F1_const <- seq((min_rel_F1_M * mean_nat_mort), (min_rel_F1_M * mean_nat_mort),
                      length.out = length(chngpoint_end:(n_years-1)))

# Increase Change
constant_F2 <- seq(Start_F2, Start_F2, length.out = length(Fish_Start_yr:chngpoint))
# Change point for F (Ramp)
change_F2 <- seq(Start_F2, (max_rel_F2_M * mean_nat_mort), length.out = length(chngpoint:chngpoint_end))
# Change point end (plateaus/stays constant now)
change_F2_const <- seq((max_rel_F2_M * mean_nat_mort), (max_rel_F2_M * mean_nat_mort),
                       length.out = length(chngpoint_end:(n_years-1)))

# Make F vector here
F_vec1 <- as.numeric(c(constant_F1, change_F1[-1], change_F1_const[-1]))
F_vec2 <- as.numeric(c(constant_F2, change_F2[-1], change_F2_const[-1]))

# Coerce into dataframe
fleet_change_fast <- data.frame(Fleet1 = F_vec1, Fleet2 = F_vec2, Year = 1:n_years,
                           Type = "Fast") %>% 
  pivot_longer(cols = c(Fleet1, Fleet2), names_to = "Fleet", values_to = "Fs")


### Slow --------------------------------------------------------------------

chngpoint <- 24     # Change point for fmort
chngpoint_end <- 49 # when does change point end
mean_nat_mort <- 0.108 # Mean M
min_rel_F1_M <- 0.01 # Min F relative to M
max_rel_F2_M <- 0.85  # Max F relative to M
Start_F1 <- 1 * mean_nat_mort # Start F1
Start_F2 <- 0.01 * mean_nat_mort # Start F2
Fish_Start_yr <- 0
n_years <- 50

# Decrease Change
# F stays constant initially
constant_F1 <- seq(Start_F1, Start_F1, length.out = length(Fish_Start_yr:chngpoint))
# Change point for F (Ramp)
change_F1 <- seq(Start_F1, (min_rel_F1_M * mean_nat_mort), length.out = length(chngpoint:chngpoint_end))
# Change point end (plateaus/stays constant now)
change_F1_const <- seq((min_rel_F1_M * mean_nat_mort), (min_rel_F1_M * mean_nat_mort),
                       length.out = length(chngpoint_end:(n_years-1)))

# Increase Change
constant_F2 <- seq(Start_F2, Start_F2, length.out = length(Fish_Start_yr:chngpoint))
# Change point for F (Ramp)
change_F2 <- seq(Start_F2, (max_rel_F2_M * mean_nat_mort), length.out = length(chngpoint:chngpoint_end))
# Change point end (plateaus/stays constant now)
change_F2_const <- seq((max_rel_F2_M * mean_nat_mort), (max_rel_F2_M * mean_nat_mort),
                       length.out = length(chngpoint_end:(n_years-1)))

# Make F vector here
F_vec1 <- as.numeric(c(constant_F1, change_F1[-1], change_F1_const[-1]))
F_vec2 <- as.numeric(c(constant_F2, change_F2[-1], change_F2_const[-1]))
F_vec1 <- c(F_vec1, rep(0.0010800, 20))
F_vec2 <- c(F_vec2, rep(0.0918000, 20))
# Coerce into dataframe
fleet_change_slow <- data.frame(Fleet1 = F_vec1, Fleet2 = F_vec2, Year = 1:70,
                                Type = "Slow") %>% 
  pivot_longer(cols = c(Fleet1, Fleet2), names_to = "Fleet", values_to = "Fs")


# Plot Fleet Structure Change Scenarios -------------------------------------------------------------------

# Coerce into dataframe
fleet_change_all <- rbind(fleet_change_fast, fleet_change_slow)

# Get aggregate F
agg_F <- fleet_change_all %>% 
  group_by(Type, Year) %>% 
  summarize(Fs = sum(Fs)) %>% 
  mutate(Fleet = "Aggregate")

# Bind together
fleet_change_all <- rbind(fleet_change_all, agg_F)

(fleet_plot <- ggplot() +
  geom_line(fleet_change_all, mapping = aes(x = Year, y = Fs, 
                                            color = Fleet, lty = Fleet),  size = 1.3) +
  scale_color_manual(values = c("black", "red", "blue")) +
  scale_linetype_manual(values = c(1, 2, 1)) +
  facet_wrap(~Type, ncol = 1) +
  theme_bw() +
  labs(x = "Year", y = "Fishing Mortality Rate") +
  ylim(0, 0.15) +
  theme(legend.position = "top", 
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 17)))

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
selexm1 <- cbind(1 / (1 + exp(-1 * ((bins - a50m1)/deltam1) )), "Males", "Fleet 2", age = bins)
selexm2 <- cbind(1 / (1 + exp(-1 * ((bins - a50m2)/deltam2) )) , "Males", "Fleet 1", age = bins)
selexf1 <- cbind(1 / (1 + exp(-1 * ((bins - a50f1)/deltaf1) )) , "Females", "Fleet 2", age = bins)
selexf2 <- cbind(1 / (1 + exp(-1 * ((bins - a50f2)/deltaf2) )) , "Females", "Fleet 1", age = bins)

# Bind these together
selex_logist <- data.frame(rbind(selexm1, selexm2, selexf1, selexf2), Type = "Logistic-Logistic (LL)")
colnames(selex_logist) <- c("Selex", "Sex", "Fleet", "Age", "Type")


### Logistic-Gamma ----------------------------------------------------------

# Gamma parameters
amaxm2 <- 18
deltam2 <- 7.5
amaxf2 <- 16
deltaf2 <- 8

# Get Selex here
pm2 <- (0.5 * (sqrt(amaxm2^2 + 4*deltam2^2) - amaxm2))
selex_m2 <- cbind((bins/amaxm2) ^ (amaxm2/pm2) * exp((amaxm2 - bins) / pm2), "Males", "Fleet 1", age = bins)
pf2 <- (0.5 * (sqrt(amaxf2^2 + 4*deltaf2^2) - amaxf2))
selex_f2 <- cbind((bins/amaxf2) ^ (amaxf2/pf2) * exp((amaxf2 - bins) / pf2), "Females", "Fleet 1", age = bins)


# Bind these together
selex_gamma <- data.frame(rbind(selexm1, selex_m2, selexf1, selex_f2), Type = "Logistic-Gamma (LG)")
colnames(selex_gamma) <- c("Selex", "Sex", "Fleet", "Age", "Type")

# Plot Selectivity Scenarios -------------------------------------------------------------------

# Coerce into dataframe
selex_all <- rbind(selex_logist, selex_gamma)

(selex_plot <- ggplot(selex_all, aes(x = as.numeric(Age), y = as.numeric(paste(Selex)), color = Fleet, lty = Fleet)) +
  geom_line(size = 1.5, alpha = 0.75) + 
  scale_color_manual(values = c("blue", "red")) +
  scale_linetype_manual(values = c(2,1)) +
  facet_grid(Sex~Type) +
  theme_bw() +
  labs(x = "Ages", y = "Selectivity") +
  theme(legend.position = "none",
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 17)))

# Get Harvest Rates -------------------------------------------------------

# Load unique OMs
OMs <- c("Fast_LL", "Fast_LG", "Slow_LL", "Slow_LG")

# Harvest Rate dataframe
harv_rates_df <- data.frame()
# Effective sample size dataframe
neff_df <- data.frame()

for(i in 1:length(OMs)) {
  # Load in OMs
  load(here(om_scenario_path, OMs[i], paste(OMs[i], ".RData", sep = "")))
  
  # Get effective sample sizes
  nf_df <- reshape2::melt(oms$Input_N_Fish)
  colnames(nf_df) <- c("Year", "Fleet", "Value")
  
  # Melt to df and rename - get harvest rate
  hr_df <- reshape2::melt(oms$Harvest_Rate)
  colnames(hr_df) <- c("Year", "Fleet", "Sim", "Value")

  hr_df$OM_Scenario <- OMs[i]
  nf_df$OM_Scenario <- OMs[i]
  
  # Bind together
  harv_rates_df <- rbind(harv_rates_df, hr_df)
  neff_df <- rbind(nf_df, neff_df)
}

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
  mutate(OM_Scenario = factor(OM_Scenario, levels = c("Fast_LG", "Fast_LL",
                                                      "Slow_LG", "Slow_LL"),
                              labels = c("Fast: Logistic-Gamma (Fast_LG)",
                                         "Fast: Logistic-Logistic (Fast_LL)",
                                          "Slow: Logistic-Gamma (Slow_LG)",
                                          "Slow: Logistic-Logistic (Fast_LL)")))

(harv_rate_plot <- ggplot(harv_rates_df_sum, aes(x = Year, y = Median_HR, color = factor(Fleet),
                          fill = factor(Fleet), ymin = Lwr_95_HR, ymax = Upr_95_HR,
                          linetype = factor(Fleet))) +
  geom_ribbon(alpha = 0.35) +
  geom_line(size = 1.5) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  scale_linetype_manual(values = c(1,2)) +
  facet_wrap(~OM_Scenario) +
  labs(x = "Year", y = "Harvest Rate", fill = "Fleet",
       color = "Fleet", linetype = "Fleet") +
  theme_bw() +
  theme(legend.position = "top",
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 17)) )
  


# Changes in Effective Sample Size Plot -----------------------------------

neff_sum_df <- neff_df %>% 
  filter(OM_Scenario %in% c("Fast_LL", "Slow_LL")) %>% 
  mutate(OM_Scenario = str_remove(OM_Scenario, "_LL"))

(neff_plot <- ggplot(neff_sum_df, aes(x = Year, y = Value, color = factor(Fleet),
                        linetype = factor(Fleet))) +
  geom_line(size = 1.5) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  scale_linetype_manual(values = c(1,2)) +
  facet_wrap(~OM_Scenario, ncol = 1) +
  labs(x = "Year", y = "Multinomial Sample Size", fill = "Fleet",
       color = "Fleet", linetype = "Fleet") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 17)) )

# Combine plots -----------------------------------------------------------

pdf(here("figs", "Hypothetical_OM", "OM_Scenarios.pdf"), width = 25, height = 8)
ggarrange(neff_plot, harv_rate_plot, selex_plot, ncol = 3,
          common.legend = TRUE, legend="top",
          labels = c("A", "B", "C"),
          font.label = list(size = 25), label.x = 0.02, widths = c(0.3, 0.6, 0.6))
dev.off()             
  
