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
library(openxlsx)

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
       dpi = 300)

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

####Test for interaction (Table S1)####
Anova(model, type = 2)

# Post-hoc tests for treatment effects within each week (Table S2)
e_results <- emmeans(model, ~ treat + week,
                     adjust = "bonferroni")
summary(e_results)

result = contrast(e_results, method = "pairwise", by = "week", adjust = "bonferroni")
result

######### Save the emmeans results for emergence flux analysis (Table S2) ##########

#openxlsx::write.xlsx(result, "data/other data_output/output/data/emerg_emm.xlsx")

########Macroinverbrate as crayfish potential food source##########

mzb_taxa = readRDS("data/other data_output/cray_source_taxa.rds")

t_clrs = c("#E69F00", "#0072B2", "#009E73", "#D55E00")  

labels = c("week0" = "Week 0 (before treatment)", "week6" = "Week 6 (end of treatment)")

####Plot the crayfish potential food source taxa abundance (Figure S3)####

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
       dpi = 300, 
       units = "in")


######### Light intensity plot Figure S1#########
light_int = readRDS("data/other data_output/light_intensity.rds")

######### Plot light intensity timeseries (Figure S1) #########

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
  dpi = 300,
  units = "cm"
)


######### Environmental data summary Figure S2#########

#############Read in the data##########

env_para = readRDS("data/other data_output/env_para.rds")

######Colors for each treatment#########
my_pal <- c(
  "#0072B2",  # blue
  "#D55E00",  # red-orange
  "#009E73",  # bluish green
  "#000000"   # black
)

###########Labels#########
f_lab = c("DOM" = "Dissolved Oxgen (% active saturation)", "pH" = "pH", "Temp" =  "Temperature (\u00B0C)",
          "Day" = "Daytime", "Night" = "Nighttime")

######### Plot environmental parameters timeseries (Figure S2) #########
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
ggsave("data/other data_output/output/env_para_plot.png", 
       width = 13, 
       height = 9, 
       units = "in", 
       dpi = 300)

#############Flume flow velocity data analysis##############
flow_vel = readRDS("data/other data_output/flow_velocity.rds") |>
  dplyr::rename(distance_m = `Distance(m)` , 
                velocity_cm_s = `Velocity(cm/s)`) |>
  dplyr::summarise(mean_vel = mean(velocity_cm_s),
                   sd_vel = sd(velocity_cm_s))

print(flow_vel)

##############Crayfish mortality analysis (Table S4)##############
cray_mort = readRDS("data/other data_output/cray_mortality.rds") |> 
  dplyr::mutate(mortality_per = (died/initial)*100) |>
  dplyr::mutate(treatment = factor(treatment, 
                                   levels = c("cray", "alan+cray"))) |> 
  dplyr::mutate(sex = factor(sex, levels = c("f", "m")))


##########Summary statistics for crayfish mortality (Table S4)##########
cray_mort |>
  dplyr::group_by(treatment) |>
  dplyr::summarise(n = n(),
    mean_mort = mean(mortality_per),
    sd_mort = sd(mortality_per),
    se_mort = sd_mort/sqrt(n))

########Create a summary table for crayfish survival and mortality data (Table S4)##########
mort_tab = cray_mort |>
  dplyr::select(treatment, flume, sex, initial, survived, died, mortality_per) |> 
  dplyr::arrange(treatment, flume)

####### Save the crayfish survival and mortality data (Table S4) ######
#openxlsx::write.xlsx(mort_tab, "data/other data_output/output/data/cray_surv_mort_data.xlsx")

# Statistical analysis of crayfish mortality

######Fitted GLMM model for crayfish mortality######

mod_surv = glmmTMB(cbind(died, survived) ~ treatment + sex + (1|flume), 
                   family = binomial,
                   data = cray_surv)
summary(mod_surv)
# Check model diagnostics
simulationOutput <- simulateResiduals(fittedModel = mod_surv, n = 1000)
plot(simulationOutput)                                      
# Test for overdispersion
testDispersion(simulationOutput)
# Test for zero-inflation
testZeroInflation(simulationOutput)
# ANOVA results for crayfish mortality (Table S5)
Anova(mod_surv, type = 2)

################Other emerging aquatic insects data analysis (Table S3)##############

emerg_other = readRDS("data/other data_output/other_emergence.rds")

emerg_dat = emerg_other |>
  dplyr::group_by(treat, week, taxa) |>
  rstatix::get_summary_stats(abund_rate, type = "mean_sd") |> 
  dplyr::arrange(treat, week, taxa)


######## Save the other emerging aquatic insects data (Table S3) ##########
#openxlsx::write.xlsx(emerg_dat, "data/other data_output/data/other_emergence_data.xlsx")


#########Crayfish activity data analysis##########
cray_activity = readRDS("data/other data_output/crayfish_activity.rds")

cray_act = cray_activity |>
  #dplyr::select(-Total) |> 
  pivot_longer(cols = c("Hidden", "Active"),
               names_to = "Activity",
               values_to = "Count")

######### Plot crayfish activity timeseries (Figure S4) #########

cray_act |>
  ggplot(aes(x = Time,
             y = Count,
             color = Activity,
             group = interaction(Treatment, Activity),
             linetype = Treatment, shape = Treatment)) +
  stat_summary(fun = mean, geom = "line", linewidth = 0.5, alpha = 0.6,
               position = position_dodge(width = 0.3)) +
  stat_summary(fun = mean, geom = "point", size = 3,
               position = position_dodge(width = 0.3)) +
  stat_summary(fun.data = mean_cl_boot,
               #fun.args = list(mult = 1),
               geom = "errorbar", width = 0.2,
               position = position_dodge(width = 0.3)) +
  scale_shape_manual(values = c(16, 17)) +
  scale_color_manual(values = c("#009E73",
                     "#000000"))+
  #facet_wrap(~ Activity) +
  labs(y = "Mean activity count",
       x = "Experiment week")+
  theme_rsm() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey80", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold", size = 12),
    axis.line = element_line(color = "black"),
    strip.placement = "outside",
    axis.title.y = element_text(size = 12, color = "black", face = "bold"),
    strip.text.x = element_text(face = "bold", size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11, color = "black")
  )

####### Save the Crayfish activity Figure S4##########
ggsave("emerg_d/crayfish_activity.png",
       height = 5, width = 7,
       dpi = 300,
       units = "in")







