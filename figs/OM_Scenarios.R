# Purpose: To plot and illustrate OM scenarios
# Creator: Matthew LH. Cheng (UAF-CFOS)
# Date 3/13/23


# Set Up ------------------------------------------------------------------

library(here)
library(tidyverse)
library(ggpubr)


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

png(here("figs", "Hypothetical_OM", "Fleet_Str_Change.png"), width = 650, height = 350)
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
dev.off()


# Selectivity Differences -------------------------------------------------

bins <- 1:30

### Logistic-Logistic ----------------------------------------------------------------

# Get a50 value males
a50m1 <- 8.5
a50m2 <- 13.5
deltam1 <- 2.25
deltam2 <- 1.25

# Get a50 value females
a50f1 <- 5.5
a50f2 <- 10.5
deltaf1 <- 1.75
deltaf2 <- 1.15

# Compute selex 
selexm1 <- cbind(1 / (1 + exp(-1 * ((bins - a50m1)/deltam1) )), "Males", "Fleet 2", age = bins)
selexm2 <- cbind(1 / (1 + exp(-1 * ((bins - a50m2)/deltam2) )) , "Males", "Fleet 1", age = bins)
selexf1 <- cbind(1 / (1 + exp(-1 * ((bins - a50f1)/deltaf1) )) , "Females", "Fleet 2", age = bins)
selexf2 <- cbind(1 / (1 + exp(-1 * ((bins - a50f2)/deltaf2) )) , "Females", "Fleet 1", age = bins)

# Bind these together
selex_logist <- data.frame(rbind(selexm1, selexm2, selexf1, selexf2), Type = "Logistic-Logistic")
colnames(selex_logist) <- c("Selex", "Sex", "Fleet", "Age", "Type")


### Logistic-Gamma ----------------------------------------------------------

# Gamma parameters
amaxm2 <- 19.5
deltam2 <- 6
amaxf2 <- 16
deltaf2 <- 7

# Get Selex here
pm2 <- (0.5 * (sqrt(amaxm2^2 + 4*deltam2^2) - amaxm2))
selex_m2 <- cbind((bins/amaxm2) ^ (amaxm2/pm2) * exp((amaxm2 - bins) / pm2), "Males", "Fleet 1", age = bins)
pf2 <- (0.5 * (sqrt(amaxf2^2 + 4*deltaf2^2) - amaxf2))
selex_f2 <- cbind((bins/amaxf2) ^ (amaxf2/pf2) * exp((amaxf2 - bins) / pf2), "Females", "Fleet 1", age = bins)


# Bind these together
selex_gamma <- data.frame(rbind(selexm1, selex_m2, selexf1, selex_f2), Type = "Logistic-Gamma")
colnames(selex_gamma) <- c("Selex", "Sex", "Fleet", "Age", "Type")

# Plot Selectivity Scenarios -------------------------------------------------------------------

# Coerce into dataframe
selex_all <- rbind(selex_logist, selex_gamma)

png(here("figs", "Hypothetical_OM", "Selex_Scenario.png"), width = 650, height = 500)
(selex_plot <- ggplot(selex_all, aes(x = as.numeric(Age), y = as.numeric(paste(Selex)), color = Fleet, lty = Fleet)) +
  geom_line(size = 1.3, alpha = 0.75) + 
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
dev.off()


# Combine plots -----------------------------------------------------------

png(here("figs", "Hypothetical_OM", "OM_Scenarios.png"), width = 1000, height = 750)
ggarrange(fleet_plot, selex_plot, ncol = 2, 
          common.legend = TRUE, legend="bottom",
          labels = c("A", "B"), 
          font.label = list(size = 25), label.x = 0.02, widths = c(0.8, 1.5))
dev.off()             
  
