################################################################
# Code for analyzing and plotting emergence flux data,
# crayfish potential food source taxa, light intensity,
# and environmental parameters from an ecological experiment.
################################################################

######Load required packages##########
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(car)

########Read in the emergence data##########

emerg_r = readRDS("data/other data_output/emergence.rds")

########Labels for facets#########

em_lab = c("Week 1" = "Week 1 (1 week after treatment start)", "Week 4" = "Week 4 (4 weeks after treatment start)")

########Plot emergence flux Figure 5##########

emerg_r |>
  dplyr::group_by(treat, week) |>
  rstatix::get_summary_stats(abund_rate, type = "mean_se") |>
  dplyr::mutate(letter = case_when(
    treat == "control" & week == "Week 1" ~ "a",  # Assign 'a' to Control group in Family 2
    treat == "control" & week == "Week 4" ~ "a",  # Assign 'b' to Treated group in Family 2
    treat == "alan+cray" & week == "Week 1" ~ "b",
    treat == "alan+cray" & week == "Week 4" ~ "b",
    treat == "cray" & week == "Week 1" ~ "b",
    treat == "cray" & week == "Week 4" ~ "b",
    TRUE ~ ""
  )) |>
  ggplot(aes(x = treat, y = mean)) +
  geom_point(size = 1.5) +  # Add points
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.1) +  # Error bars
  geom_text(aes(label = letter, y = mean + se+3.3), size = 5, color = "black")+
  scale_x_discrete(labels = c("Control", "ALAN", "Crayfish", "ALAN + Crayfish"))+
  labs(x = "Treatment",
       y = expression(bold("Emergence flux (" * ind ~ m^{-2} ~ day^{-1} * ")"))) +
  ggh4x::facet_wrap2(vars(week), labeller = as_labeller(em_lab))+
  theme_rsm() +
  theme(
    strip.text = element_text(face = "bold", size = 11, family = "Helvetica"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 11),
    axis.title.x = element_text(face = "bold", size = 11, color = "black"),
    axis.line = element_line(color = "black"),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text.y.left = element_text(angle = 0, face = "bold"),
    strip.text.x = element_text(face = "bold", size = 11)
  )

#######Save the Figure 5######
ggsave("data/other data_output/output/emergence_flux.png",
       height = 5,
       width = 7,
       units = "in",
       dpi = 600)

# Statistical analysis
# Fit a basic GLMM without autocorrelation structure
mod1 = glmmTMB(abund_rate ~ treat + week +
                 (1|flume),
               family = Gamma(link = "log"),
               data = emerg_r)
summary(mod1)
res <- residuals(mod1)
auto_cor <- acf(res, lag.max=40, plot=FALSE)
plot(auto_cor, main="ACF of Residuals", xlab="Lag", ylab="ACF")

Box.test(res, lag = 20, type = "Ljung-Box")

# Fit the model with AR(1) correlation structure

model <- glmmTMB(abund_rate ~ treat + week +
                   (1 | flume) + ar1(week + 0 |flume),
                 #ziformula = ~1,
                 dispformula = ~week,
                 data = emerg_r,
                 family = Gamma(link = "log"))  # Negative binomial family for count data
summary(model)
testDispersion(model)
# Check model diagnostics
simulationOutput <- simulateResiduals(fittedModel = model, n = 1000)
plot(simulationOutput)
# Test for overdispersion
testDispersion(simulationOutput)
# Test for zero-inflation
#testZeroInflation(simulationOutput)

####test for interaction####
Anova(model, type = 2)
# Post-hoc tests
e_results <- emmeans(model, ~ treat + week,
                     adjust = "bonferroni")
summary(e_results)

result = contrast(e_results, method = "pairwise", by = "week", adjust = "bonferroni")
result
# Plotting the results
openxlsx::write.xlsx(result, "data/other data_output/output/emerg_emm.xlsx")

########Macroinverbrate as crayfish potential food source##########

mzb_taxa = readRDS("data/other data_output/cray_source_taxa.rds")

t_clrs = c("#E69F00", "#0072B2", "#009E73", "#D55E00")  

labels = c("week0" = "Week 0 (before treatment)", "week6" = "Week 6 (end of treatment)")

####Plot the crayfish potential food source taxa abundance####

mzb_taxa |>
  dplyr::filter(!taxa %in% c("Other", "Trichoptera")) |> 
  dplyr::mutate(taxa = factor(taxa, 
                              levels = c("Chironomidae", "Gammaridae", 
                                         "Ephemeroptera", "Simuliidae"))) |> 
  ggplot(aes(Treatment, abundance, fill = taxa)) + 
  geom_boxplot(width = 0.6, position = position_dodge(width = 0.8))+
  scale_fill_manual(values = t_clrs)+
  ggh4x::facet_wrap2(vars(time), labeller = as_labeller(labels))+
  labs(y = expression(bold("Total abundance (ind. m"^-2*")")),
       x = "Treatment",
       fill = "Crayfish potential food source")+
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        strip.background = element_rect(fill = "grey80", color = NA),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
        strip.text = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.line = element_line(color = "black"),
        strip.placement = "outside",  
        strip.text.y.left = element_text(angle = 0, face = "bold"), 
        strip.text.x = element_text(face = "bold") 
  )

####### Save the Figure S3 ######
ggsave("data/other data_output/output/mzb_source_taxa.png",
       height = 6, width = 8,
       dpi = 600, 
       units = "in")


######### Light intensity plot Figure S1#########
light_int = readRDS("data/other data_output/light_intensity.rds")

######### Plot light intensity timeseries #########

light_int |>
  dplyr::group_by(date, day_night, Treatment) %>%
  dplyr::summarise(
    total_light = sum(intensity * 1800, na.rm = TRUE), .groups = "drop"  # Lux * seconds
  ) %>%
  dplyr::mutate(total_light_hours = total_light / 3600) |>
  dplyr::ungroup() |>
  ggplot(aes(x = date, y = total_light_hours)) +
  geom_line() +
  facet_grid(day_night~Treatment, scale = "free_y") +
  labs(x = "Date",
       y = expression(bold("Total light exposure (Lux·hours)"))) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.line = element_line(color = "black"),
    strip.background = element_rect(fill = "grey90", color = "black"),  
    strip.placement = "outside",  
    strip.text.y.left = element_text(angle = 0, face = "bold"), 
    strip.text.x = element_text(face = "bold"))

####### Save the Figure S1######

ggsave(
  "data/other data_output/output/light_intensity.png",
  height = 12,
  width = 16.5,
  units = "cm"
)


######### Environmental data summary Figure S2#########

#############Read in the data##########

env_para = readRDS("data/other data_output/env_para.rds")

######Colors for each treatment#########
my_pal <- c(
  "#0072B2",  # blue
  "#D55E00",  # vermillion (red-orange)
  "#009E73",  # bluish green
  "#000000"   # orange (safe orange, not yellow)
)

###########Labels#########
f_lab = c("DOM" = "Dissolved Oxgen (% active saturation)", "pH" = "pH", "Temp" =  "Temperature (\u00B0C)",
          "Day" = "Daytime", "Night" = "Nighttime")

######### Plot environmental parameters timeseries #########
env_para |> 
  dplyr::select(Flume,Treatment,date, time, hour,DOM, Temp, pH) |> 
  pivot_longer(6:8,
               names_to = "phch",
               values_to = "value") |>
  dplyr::select(Treatment, time, hour, date, phch, value) |> 
  dplyr::mutate(day_night = if_else(
    time >= hms::as_hms("05:00:00") & time < hms::as_hms("21:30:00"),
    "Day",
    "Night"
  )) |> 
  mutate(across(c("phch", "day_night"), as.factor)) |>
  dplyr::group_by(Treatment,date,day_night, phch) |> 
  dplyr::summarise(mean = mean(value), .groups = "drop") |>
  ggplot(aes(x = date, y = mean, color = Treatment, group = Treatment)) +
  geom_line(linewidth = 0.8) + 
  scale_color_manual(values = my_pal)+
  ggh4x::facet_grid2(day_night ~ phch, scales = "free_y", independent = "y",
                     labeller = as_labeller(f_lab))+
  labs(
    x = "Date",
    y = "Mean environmental parameters")+
  theme_minimal(base_size = 14)+
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color = "black"), 
    axis.title.y = element_text(),
    axis.title = element_text(size = 14, color = "black", face = "bold"),
    axis.text = element_text(size = 14, color = "black", face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(size = 14, color = "black", face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    plot.tag = element_text(size = 14, face = "bold"),
    strip.text = element_text(face = "bold", size = 14),
    plot.margin = margin(t = 10, r = 20, b = 10, l = 10), 
    strip.background = element_rect(fill = "grey90", color = "black")
  )

####### Save the Figure S2######
ggsave("data/other data_output/output/env_para_plot.png", width = 13, height = 9, units = "in", dpi = 600)

