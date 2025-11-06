library(glmmTMB)
library(DHARMa)
library(emmeans)
library(car)

abund = readxl::read_excel("data/Emergence/Emergence for SIA.xlsx", sheet = "Abundance") |>
  pivot_longer(
    7:25,
    names_to = "Family",
    values_to = "abundance"
  )


bio = readxl::read_excel("data/Emergence/Emergence for SIA.xlsx", sheet = "Biomass (mg)") |> 
  pivot_longer(
    6:26,
    names_to = "Family",
    values_to = "biomass"
  )


#######Percentage of abundance of all emerging aquatic insects######
ab_perc = abund |>
  dplyr::group_by(Family) |>
  dplyr::summarise(total_abund = sum(abundance, na.rm = TRUE), .groups = "drop") |> 
  dplyr::mutate(s_abund = sum(total_abund, na.rm = TRUE)) |> 
  #dplyr::left_join(abund, by = c("Family")) |> 
  dplyr::mutate(abund_perc = total_abund/s_abund*100)

##################################

b = bio |>
  # dplyr::mutate(across(c(Treatment, Week), factor)) |> 
  # dplyr::mutate(Week = case_when(
  #   Week == 1 ~ "Wk 1",
  #   Week == 4 ~ "Wk 4",
  #   TRUE ~ Week
  # )) |>
  # dplyr::mutate(Treatment = factor(Treatment, levels = c("control", "alan", "cray", "alan+cray"))) |> 
  #dplyr::filter(!biomass %in% 0 & Family %in% "Chironomidae") |> 
  dplyr::group_by(Treatment, Week) |>
  rstatix::get_summary_stats(biomass, type = "mean_se") |> 
  ggplot(aes(x = Treatment, y = mean, color = Treatment)) +
  #geom_line(linewidth = 0.5) +  # Line plot
  geom_point( size = 1.5) +  # Add points
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.1) +  # Error bars
  #scale_color_brewer(palette = "Dark2")+
  #scale_color_manual(values = c("Wk 1" = "#E69F00", "Wk 4" = "#56B4E9"))+
  scale_x_discrete(labels = c("Control", "ALAN", "Crayfish", "ALAN + Crayfish"))+
  #ylim(0,35)+
  labs(tag = "b",
       x = "Treatment",
       y = "Mean total biomass") +
  facet_wrap(~Week)+
  theme_minimal() +
  theme(legend.position = "bottom")  # Rotate x-axis labels for readability
b

f_lab = c("Week 1" = "Week 1 (1 week after treatment start)", "Week 4" = "Week 4 (4 weeks after treatment start)")

sp_l = sp |> 
  mutate(letter = case_when(
    Family == "Tetragnatha" & Treatment == "Control" ~ "a",  # Assign 'a' to Control group in Family 2
    Family == "Tetragnatha" & Treatment == "Low flow" ~ "b",  # Assign 'b' to Treated group in Family 2
    TRUE ~ ""  # Other cases remain blank (no significant difference)
  ))

#saveRDS(abund_r, "data/Emergence/emergence.rds")

emerg_r = readRDS("data/Emergence/emergence.rds")

emerg_r |>
  dplyr::group_by(treat, week) |>
  rstatix::get_summary_stats(abund_rate, type = "mean_se") |>
  dplyr::mutate(letter = case_when(
    treat == "control" & week == "Week 1" ~ "a",  # Assign 'a' to Control group in Family 2
    treat == "control" & week == "Week 4" ~ "a",  # Assign 'b' to Treated group in Family 2
    treat == "alan+cray" & week == "Week 1" ~ "b",
    treat == "cray" & week == "Week 4" ~ "b",
    TRUE ~ ""
  )) |>
  ggplot(aes(x = treat, y = mean)) +
  geom_point(size = 1.5) +  # Add points
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.1) +  # Error bars
  geom_text(aes(label = letter, y = mean + se+3.5), size = 5, color = "black")+
  scale_x_discrete(labels = c("Control", "ALAN", "Crayfish", "ALAN + Crayfish"))+
  labs(x = "Treatment",
    y = expression(bold("Emergence flux (" * ind ~ m^{-2} ~ day^{-1} * ")"))) +
  facet_wrap(~week)+
  ggh4x::facet_wrap2(vars(week), labeller = as_labeller(f_lab))+
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

a = abund_r |>
  # dplyr::mutate(across(c(treat, week), factor)) |>
  # dplyr::mutate(week = case_when(
  #   week == "w1" ~ "Week 1",
  #   week == "w4" ~ "Week 4",
  #   TRUE ~ week
  # )) |> 
  # dplyr::mutate(treat = factor(treat, levels = c("control", "alan", "cray", "alan+cray"))) |> 
  #dplyr::filter(!abundance %in% 0 & Family %in% c("chiro", "tany")) |>
  #dplyr::filter(!abundance %in% 0 & Family %in% c("hepta", "other.ephe", "simu", "chiro", "tany")) |> 
  dplyr::group_by(treat, week) |>
  rstatix::get_summary_stats(abund_rate, type = "mean_se") |>
  dplyr::mutate(letter = case_when(
    treat == "control" & week == "Week 1" ~ "a",  # Assign 'a' to Control group in Family 2
    treat == "control" & week == "Week 4" ~ "a",  # Assign 'b' to Treated group in Family 2
    treat == "alan+cray" & week == "Week 1" ~ "b",
    treat == "cray" & week == "Week 4" ~ "b",
    TRUE ~ ""
  )) |>
  ggplot(aes(x = treat, y = mean)) +
  #geom_boxplot(size = 0.5, aes(group = treat)) +  # Line plot
  geom_point(size = 1.5) +  # Add points
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.1) +  # Error bars
  geom_text(aes(label = letter, y = mean + se+3.5), size = 5, color = "black")+
  #scale_color_manual(values = c("Wk 1" = "#E69F00", "Wk 4" = "#56B4E9"))+
  #scale_color_brewer(palette = "Dark2")+
  #scale_color_viridis_d(option = "C")+
  scale_x_discrete(labels = c("Control", "ALAN", "Crayfish", "ALAN + Crayfish"))+
  #ylim(0,100)+
  labs(#tag = "a",
    x = "Treatment",
    y = expression(bold("Emergence flux (" * ind ~ m^{-2} ~ day^{-1} * ")"))) +
  facet_wrap(~week)+
  ggh4x::facet_wrap2(vars(week), labeller = as_labeller(f_lab))+
  theme_rsm() +
  theme(
    #panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 11, family = "Helvetica"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 11),
    axis.title.x = element_text(face = "bold", size = 11, color = "black"),
    axis.line = element_line(color = "black"),
    strip.background = element_rect(fill = "grey90", color = "black"),  # box for facet labels
    strip.placement = "outside",  # place the day_night labels on the right
    strip.text.y.left = element_text(angle = 0, face = "bold"),  # day/night on the right
    strip.text.x = element_text(face = "bold", size = 11)  # treatment labels
    
  )

a

#######Save######
ggsave("col_output/emergence.png", a,
       height = 5,
       width = 7,
       units = "in",
       dpi = 600)


##
#Combining the two plot



(a | b) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom") 

mixture_crayfish

install.packages("glmmTMB","DHARMa","emmeans", dependencies=TRUE)





# Fit the model
abund_r = abund |> 
  dplyr::mutate(across(c(treat, week), factor)) |>
  dplyr::mutate(week = case_when(
    week == "w1" ~ "Week 1",
    week == "w4" ~ "Week 4",
    TRUE ~ week
  )) |> 
  dplyr::select(treat, week, flume, trap, Family, abundance, day) |> 
  dplyr::mutate(treat = factor(treat, levels = c("control", "alan", "cray", "alan+cray"))) |> 
  dplyr::mutate(week = factor(week, levels = c("Week 1", "Week 4"))) |> 
  dplyr::mutate(flume = factor(flume)) |>
  dplyr::mutate(day = factor(day)) |>
  dplyr::mutate(trap_day = 4) |>
  dplyr::mutate(area = 1) |>
  dplyr::mutate(abund_rate = abundance/(trap_day*area)) |>   # Convert
  #dplyr::filter(!abundance %in% 0)
  dplyr::filter(!abundance %in% 0 & Family %in% "chiro")
#"hepta", "other.ephe", "simu", "chiro", "tany"
# Fit the model



bio_r = bio |> 
  dplyr::mutate(across(c(Treatment, Week), factor)) |>
  dplyr::mutate(Week = case_when(
    Week == "1" ~ "Week 1",
    Week == "4" ~ "Week 4",
    TRUE ~ Week
  )) |> 
  #dplyr::select(treat, week, flume, trap, Family, abundance, day) |> 
  dplyr::mutate(Treatment = factor(Treatment, levels = c("control", "alan", "cray", "alan+cray"))) |> 
  dplyr::mutate(Week = factor(Week, levels = c("Week 1", "Week 4"))) |> 
  dplyr::mutate(Flume = factor(Flume)) |>
  dplyr::mutate(Day = factor(Day)) |>
  dplyr::mutate(trap_day = 4) |>
  dplyr::mutate(area = 1) |>
  dplyr::mutate(mass_flux = biomass/(trap_day*area)) |>   # Convert
  #dplyr::filter(!abundance %in% 0)
  dplyr::filter(!biomass %in% 0 & Family %in% "Chironomidae")




mod1 = glmmTMB(abund_rate ~ treat * week +
                (1|flume),
              family = Gamma(link = "log"),
              data = abund_r)
summary(mod1)
res <- residuals(mod1)
auto_cor <- acf(res, lag.max=40, plot=FALSE)
plot(auto_cor, main="ACF of Residuals", xlab="Lag", ylab="ACF")

Box.test(res, lag = 20, type = "Ljung-Box")

# Fit the model with AR(1) correlation structure

model <- glmmTMB(abund_rate ~ treat * week + (1 | flume)
                 + ar1(week + 0 |flume),
                 #ziformula = ~1,
                 dispformula = ~week,
                 data = abund_r,
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
e_results <- emmeans(model, ~ treat * week,
                     adjust = "bonferroni")
summary(e_results)

result = contrast(e_results, method = "pairwise", by = "week", adjust = "bonferroni")
# Plotting the results
openxlsx::write.xlsx(result, "data/Emergence/Emmean_output/emerg_emm.xlsx")

plot(emmeans_results, comparisons = TRUE)
