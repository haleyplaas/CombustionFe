setwd("C:\\Users\\haley\\OneDrive - North Carolina State University\\Code Repositories\\R\\coalflyash")
rm(list = ls())
#.rs.restartR() # turn on when need to clear up RAM space

library(dplyr);library(ggplot2);library(readxl);library(patchwork);library(ggpmisc);library(ggpattern);library(forcats) 

coal.ash.sims.3 <- read_excel("D:\\NCSU\\Collaborations\\Coal Fly Ash\\data_for_pub\\excel_files\\FeDepositionBudgets_Ocean_Final.xlsx", sheet = "ocean_dep_soluble_only", na = c("NA",""))

coal.ash.sims.3.filtered <- coal.ash.sims.3 %>% 
 # filter(Variable != "FESOLDEP_mean" , Variable !="FETOTDEP_mean") %>% 
  mutate(settings_combo = paste(Scenario, Simulation, sep = "_")) 

simulation.color.palette <- c("PD_V1"   = "#a8ffa8",
                              "PD_V3" = "#04ff00",
                              "PD_V2" = "#04ff90",
                              "PD_V4" = "#029000",
                             # "PD_V4" = "#0c3b0c"
                             # "PD_V4" = "#0c6b0d", 
                              "PI_V1" = "#1A1A1A80", 
                              "PI_V3" = "#2e2e2e", 
                              "SSP370mid_V1" = "#cc5a0080", 
                              "SSP370mid_V3" = "#cc5a00",
                              "SSP370end_V1" = "#a93b0080",
                              "SSP370end_V3" = "#a93b00")

period.color.palette <- c("PI"= "#000000",
                          "PD"= "#0a4c00",
                          "SSP370mid"= "#ab4c00",
                          "SSP370end"= "#832e00")

period.color.palette2 <- c("PI"= "white",
                          "PD"= "white",
                          "SSP370mid"= "white",
                          "SSP370end"= "white")

aerosol.source.patterns <- c("FEDUSOLDEP_mean" = "stripe",
                             "FEDUTOTDEP_mean" = "stripe",
                             "FEANSOLDEP_mean"   = "circle",
                             "FEANTOTDEP_mean"   = "circle",
                             "FEBBSOLDEP_mean"   = NA,
                             "FEBBTOTDEP_mean"   = NA)


coal.ash.sims.3.soluble <- coal.ash.sims.3.filtered %>% mutate(settings_combo = fct_relevel(settings_combo, "PI_V1", "PI_V3", "PD_V1", "PD_V2", "PD_V3", "PD_V4", "SSP370mid_V1", "SSP370mid_V3", "SSP370end_V1", "SSP370end_V3"))

#coal.ash.sims.3.total <- coal.ash.sims.3.filtered %>% 
#  filter(!grepl("SOL", Variable)) %>%
 # mutate(settings_combo = fct_relevel(settings_combo, "PI_V1", "PI_V3", "PD_V1", "PD_V2", "PD_V3", "PD_V4", "SSP370mid_V1", #"SSP370mid_V3", "SSP370end_V1", "SSP370end_V3"))

#ggplot(coal.ash.sims.3.total, aes(x = settings_combo, y = Ocean_budget, 
                                #  fill = settings_combo, color = Scenario, 
                                #  pattern = Variable, pattern_fill = Scenario)) +  # Map pattern fill
  #geom_bar_pattern(stat = "identity", position = "stack", linewidth = 1.2, 
                  # pattern_density = 0.3) +  # Outline color of the pattern  
 # scale_fill_manual(values = simulation.color.palette) +
 # scale_color_manual(values = period.color.palette) +  
 # scale_pattern_manual(values = aerosol.source.patterns) +  
 # scale_pattern_fill_manual(values = period.color.palette) +  # Map pattern fill
 # theme_bw() # + 
 # ylim(c(0,38))

#GLOBAL
ggplot(coal.ash.sims.3.soluble, aes(x = settings_combo, y = Ocean_budget,
                                                  fill = settings_combo, color = Scenario, 
                                                  pattern = Variable, pattern_fill = Scenario)) +  # Map pattern fill
  geom_bar_pattern(stat = "identity", position = "stack", linewidth = 1.2, 
                   pattern_density = 0.3) +  # Outline color of the pattern  
  scale_fill_manual(values = simulation.color.palette) +
  scale_color_manual(values = period.color.palette) +  
  scale_pattern_manual(values = aerosol.source.patterns) +  
  scale_pattern_fill_manual(values = period.color.palette) +  # Map pattern fill
  theme_bw()  + 
  ylim(c(0,0.5)) +
  ggtitle("Global") +
  theme(legend.position = "none")

#SO
ggplot(coal.ash.sims.3.soluble, aes(x = settings_combo, y = SO_Budget, 
                                          fill = settings_combo, color = Scenario, 
                                          pattern = Variable, pattern_fill = Scenario)) +  # Map pattern fill
  geom_bar_pattern(stat = "identity", position = "stack", linewidth = 1.2, 
                   pattern_density = 0.3) +  # Outline color of the pattern  
  scale_fill_manual(values = simulation.color.palette) +
  scale_color_manual(values = period.color.palette) +  
  scale_pattern_manual(values = aerosol.source.patterns) +  
  scale_pattern_fill_manual(values = period.color.palette) +  # Map pattern fill
  theme_bw() + 
  ylim(c(0,0.04)) +
  ggtitle("SO") +
  theme(legend.position = "none")

#SEAS
ggplot(coal.ash.sims.3.soluble, aes(x = settings_combo, y = SEAS_Budget, 
                                            fill = settings_combo, color = Scenario, 
                                            pattern = Variable, pattern_fill = Scenario)) +  # Map pattern fill
  geom_bar_pattern(stat = "identity", position = "stack", linewidth = 1.2, 
                   pattern_density = 0.3) +  # Outline color of the pattern  
  scale_fill_manual(values = simulation.color.palette) +
  scale_color_manual(values = period.color.palette) +  
  scale_pattern_manual(values = aerosol.source.patterns) +  
  scale_pattern_fill_manual(values = period.color.palette) +  # Map pattern fill
  theme_bw() + 
  ylim(c(0,0.04)) +
  ggtitle("SEAS") +
  theme(legend.position = "none")

#BB
ggplot(coal.ash.sims.3.soluble, aes(x = settings_combo, y = BB_Budget, 
                                          fill = settings_combo, color = Scenario, 
                                          pattern = Variable, pattern_fill = Scenario)) +  # Map pattern fill
  geom_bar_pattern(stat = "identity", position = "stack", linewidth = 1.2, 
                   pattern_density = 0.3) +  # Outline color of the pattern  
  scale_fill_manual(values = simulation.color.palette) +
  scale_color_manual(values = period.color.palette) +  
  scale_pattern_manual(values = aerosol.source.patterns) +  
  scale_pattern_fill_manual(values = period.color.palette) +  # Map pattern fill
  theme_bw() + 
  ylim(c(0,0.04)) +
  ggtitle("BB") +
  theme(legend.position = "none")

#AUSP
ggplot(coal.ash.sims.3.soluble, aes(x = settings_combo, y = AUSP_Budget, 
                                            fill = settings_combo, color = Scenario, 
                                            pattern = Variable, pattern_fill = Scenario)) +  # Map pattern fill
  geom_bar_pattern(stat = "identity", position = "stack", linewidth = 1.2, 
                   pattern_density = 0.3) +  # Outline color of the pattern  
  scale_fill_manual(values = simulation.color.palette) +
  scale_color_manual(values = period.color.palette) +  
  scale_pattern_manual(values = aerosol.source.patterns) +  
  scale_pattern_fill_manual(values = period.color.palette) +  # Map pattern fill
  theme_bw() + 
  ylim(c(0,0.04)) +
  ggtitle("AUSP") +
  theme(legend.position = "none")

#NATL
ggplot(coal.ash.sims.3.soluble, aes(x = settings_combo, y = NATL_Budget, 
                                            fill = settings_combo, color = Scenario, 
                                            pattern = Variable, pattern_fill = Scenario)) +  # Map pattern fill
  geom_bar_pattern(stat = "identity", position = "stack", linewidth = 1.2, 
                   pattern_density = 0.3) +  # Outline color of the pattern  
  scale_fill_manual(values = simulation.color.palette) +
  scale_color_manual(values = period.color.palette) +  
  scale_pattern_manual(values = aerosol.source.patterns) +  
  scale_pattern_fill_manual(values = period.color.palette) +  # Map pattern fill
  theme_bw() + 
  ylim(c(0,0.04)) +
  ggtitle("NATL") +
  theme(legend.position = "none")

#SATL
ggplot(coal.ash.sims.3.soluble, aes(x = settings_combo, y = SATL_Budget, 
                                    fill = settings_combo, color = Scenario, 
                                    pattern = Variable, pattern_fill = Scenario)) +  # Map pattern fill
  geom_bar_pattern(stat = "identity", position = "stack", linewidth = 1.2, 
                   pattern_density = 0.3) +  # Outline color of the pattern  
  scale_fill_manual(values = simulation.color.palette) +
  scale_color_manual(values = period.color.palette) +  
  scale_pattern_manual(values = aerosol.source.patterns) +  
  scale_pattern_fill_manual(values = period.color.palette) +  # Map pattern fill
  theme_bw() + 
  ylim(c(0,0.04)) +
  ggtitle("SATL") +
  theme(legend.position = "none")

#NPAC
ggplot(coal.ash.sims.3.soluble, aes(x = settings_combo, y = NPAC_Budget, 
                                            fill = settings_combo, color = Scenario, 
                                            pattern = Variable, pattern_fill = Scenario)) +  # Map pattern fill
  geom_bar_pattern(stat = "identity", position = "stack", linewidth = 1.2, 
                   pattern_density = 0.3) +  # Outline color of the pattern  
  scale_fill_manual(values = simulation.color.palette) +
  scale_color_manual(values = period.color.palette) +  
  scale_pattern_manual(values = aerosol.source.patterns) +  
  scale_pattern_fill_manual(values = period.color.palette) +  # Map pattern fill
  theme_bw() + 
  ylim(c(0,0.04)) +
  ggtitle("NPAC") +
  theme(legend.position = "none")

#AS
ggplot(coal.ash.sims.3.soluble, aes(x = settings_combo, y = AS_Budget, 
                                    fill = settings_combo, color = Scenario, 
                                    pattern = Variable, pattern_fill = Scenario)) +  # Map pattern fill
  geom_bar_pattern(stat = "identity", position = "stack", linewidth = 1.2, 
                   pattern_density = 0.3) +  # Outline color of the pattern  
  scale_fill_manual(values = simulation.color.palette) +
  scale_color_manual(values = period.color.palette) +  
  scale_pattern_manual(values = aerosol.source.patterns) +  
  scale_pattern_fill_manual(values = period.color.palette) +  # Map pattern fill
  theme_bw() + 
  ylim(c(0,0.2)) +
  ggtitle("AS") +
  theme(legend.position = "none")

#ENPAC
ggplot(coal.ash.sims.3.soluble, aes(x = settings_combo, y = ENPAC_Budget, 
                                    fill = settings_combo, color = Scenario, 
                                    pattern = Variable, pattern_fill = Scenario)) +  # Map pattern fill
  geom_bar_pattern(stat = "identity", position = "stack", linewidth = 1.2, 
                   pattern_density = 0.3) +  # Outline color of the pattern  
  scale_fill_manual(values = simulation.color.palette) +
  scale_color_manual(values = period.color.palette) +  
  scale_pattern_manual(values = aerosol.source.patterns) +  
  scale_pattern_fill_manual(values = period.color.palette) +  # Map pattern fill
  theme_bw() + 
  ylim(c(0,0.04)) +
  ggtitle("ENPAC") +
  theme(legend.position = "none")

#WNPAC
ggplot(coal.ash.sims.3.soluble, aes(x = settings_combo, y = WNPAC_Budget, 
                                    fill = settings_combo, color = Scenario, 
                                    pattern = Variable, pattern_fill = Scenario)) +  # Map pattern fill
  geom_bar_pattern(stat = "identity", position = "stack", linewidth = 1.2, 
                   pattern_density = 0.3) +  # Outline color of the pattern  
  scale_fill_manual(values = simulation.color.palette) +
  scale_color_manual(values = period.color.palette) +  
  scale_pattern_manual(values = aerosol.source.patterns) +  
  scale_pattern_fill_manual(values = period.color.palette) +  # Map pattern fill
  theme_bw() + 
  ylim(c(0,0.04)) +
  ggtitle("WNPAC") +
  theme(legend.position = "none")

#CPAO
ggplot(coal.ash.sims.3.soluble, aes(x = settings_combo, y = CPAO_Budget, 
                                    fill = settings_combo, color = Scenario, 
                                    pattern = Variable, pattern_fill = Scenario)) +  # Map pattern fill
  geom_bar_pattern(stat = "identity", position = "stack", linewidth = 1.2, 
                   pattern_density = 0.3) +  # Outline color of the pattern  
  scale_fill_manual(values = simulation.color.palette) +
  scale_color_manual(values = period.color.palette) +  
  scale_pattern_manual(values = aerosol.source.patterns) +  
  scale_pattern_fill_manual(values = period.color.palette) +  # Map pattern fill
  theme_bw() + 
  ylim(c(0,0.2)) +
  ggtitle("CPAO") +
  theme(legend.position = "none")





