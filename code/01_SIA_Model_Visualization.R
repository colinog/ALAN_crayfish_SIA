##################################################################
# MixSIAR model visualization for spider and crayfish
##################################################################

source("code/00_setup.R")

for (file in list.files("data", pattern = "*.csv")) {
  assign(
    stringr::str_remove(file, ".csv"),
    read.csv(paste0("data/", file)) %>% 
      tibble::tibble() %>%
      recode_treatment()
  )
}

#### Spider --------------
source_spider_tef_corrected <- 
  source_spider %>%
  # adjust soruces by TEF
  dplyr::left_join(TEF_spider, by = dplyr::join_by(Family)) %>%
  dplyr::mutate(d15N = d15N + Meand15N,
                d13C = d13C + Meand13C,
                flume_nr = readr::parse_number(Flume)) %>%
  dplyr::group_by(Treatment) %>%
  dplyr::mutate(flume_nr = as.numeric(factor(flume_nr)))

sstcs <- # source_spider_tef_corrected_summary
  source_spider_tef_corrected %>%
  dplyr::group_by(ALAN, Crayfish, Family, flume_nr) %>%
  dplyr::summarise(
    dplyr::across(
      c(d15N, d13C),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd = ~sd(.x, na.rm = TRUE)
      )
    ), .groups = "drop"
  )

##### Spider polgyon Plot -----------------------

plt_tet = source_spider_tef_corrected %>%
  dplyr::rename(Taxa = Family) %>%
  dplyr::mutate(Treatment = factor(Treatment, levels = c( "Control","ALAN", "Crayfish", "ALAN + Crayfish"))) |> 
  dplyr::group_by(Treatment,Flume, Taxa) |>
  dplyr::summarise(
    d13C = mean(d13C, na.rm = TRUE),
    d15N = mean(d15N, na.rm = TRUE), .groups = "drop"
  ) |>
  ggplot(aes(x = d13C, y = d15N, color = Taxa, fill = Taxa)) +
  ggpubr::stat_chull(
    geom = "polygon",
    color = "transparent",
    fill = "#000000",
    alpha = 0.15
  ) +
  geom_point(
    size = 1.5,
    alpha = 0.8
  ) +
  geom_point(
    data = mixture_spider %>%
      dplyr::rename(Taxa = Family) %>%
      dplyr::mutate(Treatment = factor(Treatment, levels = c( "Control","ALAN", "Crayfish", "ALAN + Crayfish"))) |> 
      dplyr::group_by(Treatment,Flume, Taxa) |> 
      dplyr::summarise(
        d13C = mean(d13C, na.rm = TRUE),
        d15N = mean(d15N, na.rm = TRUE), .groups = "drop") |> 
      dplyr::mutate(flume_nr = readr::parse_number(Flume)) %>%
      dplyr::group_by(Treatment) %>%
      dplyr::mutate(flume_nr = as.numeric(factor(flume_nr))),
    color = "#000000",
    fill = "transparent",
    size = 2,
    alpha = 0.6
  ) +
  scale_color_manual(values = col_src_spider, breaks = names(col_src_spider)) +
  scale_fill_manual(values = col_src_spider, breaks = names(col_src_spider)) +
  ggh4x::facet_nested(Treatment~.) +
  #ggh4x::facet_nested(ALAN + Crayfish~"Replicate" ) +
  labs(
    x = "&delta;<sup>13</sup>C (&permil;)",
    y = "&delta;<sup>15</sup>N (&permil;)"
  ) +
  theme_rsm() +
  guides(color = guide_legend(nrow = 2)) +
  theme(
    strip.text.y = element_blank(),
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown(),
    legend.title.position = "top",
    legend.direction = "horizontal",
    legend.title = element_text(hjust = 0.5),
    legend.position = "top",
    legend.key.size = unit(0.5,"line"),
  )

ggsave(
  "output/spider/spider_polygon.png",
  plt_tet,
  height = 12,
  width = 16.5,
  units = "cm"
)

##### Run MixSIAR for spider --------
set.seed(123)
spider_models <-
  unique(mixture_spider$Treatment) %>%
  purrr::map(function(x) {
    run_spider_model(
      run = "short",
      treatment = x,
      mix = mixture_spider,
      src = source_spider,
      discr = TEF_spider
    )
  })

##########Read in pre-run saved spider models ############
spider_models = readRDS("rds_data/spider/spider_model.rds")

##### Spider model trace plots for model diagnostics --------------
spider_models %>%
  purrr::walk(
    function(x) {
      plt_temp <-
        x$draws %>%
        dplyr::mutate(iter = 1:n()) %>%
        tidyr::pivot_longer(-c(Treatment, iter)) %>%
        ggplot(aes(iter, value, color = name)) +
        geom_line(linewidth = 0.1) +
        scale_color_manual(values = col_src_spider) +
        facet_grid(name~Treatment) +
        theme_rsm()
      
      ggsave(
        paste0("output/spider_trace/spider_", unique(x$draws$Treatment), "_trace.png"),
        plt_temp,
        height = 10,
        width = 16.5,
        units = "cm"
      )
    }
  )

##### Spider Diet Proportion Plot ---------------------
plt_diet_spider <-
  dplyr::bind_rows(
    spider_models %>%
      purrr::map_df(~.x$draws %>% dplyr::mutate(iter = 1:n()))
  ) %>%
  dplyr::mutate(Aquatic = Chironomidae+Tanypodinae,
                Terrestrial = Cicadellidae+Curculionidae) %>%
  dplyr::select(Treatment:Terrestrial) %>%
  tidyr::pivot_longer(-c(Treatment, iter)) %>%
  dplyr::left_join(
    mixture_spider %>%
      dplyr::select(Treatment, ALAN, Crayfish) %>%
      dplyr::distinct(),
    by = "Treatment"
  ) %>%
  dplyr::mutate(Treatment = factor(Treatment, 
                                   levels = c("Control","ALAN","Crayfish", "ALAN + Crayfish"))) |> 
  dplyr::group_by(Treatment, name, iter) %>%
  dplyr::summarise(value = mean(value), .groups = "drop") %>%
  ggplot(aes(x = value, fill = name,  color = name)) +
  geom_density(alpha = 0.5) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_fill_manual(values = unname(col_hab)[c(4, 2)], guide = "legend", name = "Habitat",
                    limits = c("Terrestrial", "Aquatic")) +
  scale_color_manual(values = unname(col_hab)[c(3, 1)], guide = "legend", name = "Habitat",
                     limits = c("Terrestrial", "Aquatic")) +
  scale_x_continuous(breaks = 0:4/4, labels = scales::percent_format()) +
  guides(
    fill = guide_legend(override.aes = list(alpha = 0.5, color = NA))
  )+
  ggh4x::facet_nested(Treatment~1) +
  labs(y = NULL, x = "Dietary proportion") +
  theme_rsm() +
  theme(
    strip.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 33),
    axis.ticks.y = element_blank(),
    legend.key.size = unit(0.5,"line"),
    legend.title.position = "top",
    legend.direction = "horizontal",
    legend.title = element_text(hjust = 0.5),
    legend.position = "top"
  )

####### Patchwork  layout design -------------

design <- "AAABB"

############ Save combined spider isotope polygon and diet proportion plot -------------

ggsave(
  "output/spider/spider_poly_prop.png",
  plt_tet + 
    plt_diet_spider + 
    patchwork::plot_annotation(
      tag_levels = "a") +
    theme(plot.tag.position = c(0, 1),
          plot.tag = element_text(hjust = 0.5, vjust = 1)) +
    patchwork::plot_layout(design = design, guides = "keep") &
    theme(legend.position = "top",
          legend.box = "horizontal",
          legend.justification = "center",
          legend.box.spacing = unit(0.2, "cm"),
          legend.text = element_text(margin = margin(l = 8)),
          legend.spacing.x = unit(1.3, "cm"),
          plot.tag.position =  "topleft",
          plot.tag = element_text(face = "bold", size = 10)),
  height = 6,
  width = 7,
  units = "in",
  dpi = 600 
)

#####To estimate the median dietary contribution frome each source#####

#######The data frame in wide format ##########

sp_wide = dplyr::bind_rows(
  spider_models %>%
    purrr::map_df(~.x$draws %>% dplyr::mutate(iter = 1:n()))) %>%
  dplyr::mutate(Aquatic = Chironomidae+Tanypodinae,
                Terrestrial = Cicadellidae+Curculionidae) %>%
  dplyr::select(Treatment:Terrestrial) %>%
  tidyr::pivot_longer(-c(Treatment, iter)) %>%
  dplyr::left_join(
    mixture_spider %>%
      dplyr::select(Treatment,ALAN, Crayfish) %>%
      dplyr::distinct(),
    by = "Treatment"
  ) |> 
  select(iter, Treatment, name, value) %>%
  tidyr::pivot_wider(names_from = Treatment, values_from = value)

#######The data frame in long format ##########
##Treatment posterior draws - control draws
###combined treatment draws - individual treatment draws
con_long = sp_wide %>%
  dplyr::filter(!name %in% "Terrestrial") |>  # <-- Filter if needed
  mutate(
    ALAN_Control = `ALAN` - `Control`,
    Crayfish_Control = `Crayfish` - `Control`,
    `ALAN + Crayfish_Control` = `ALAN + Crayfish` - `Control`,
    `ALAN + Crayfish_ALAN` = `ALAN + Crayfish` - `ALAN`,
    `ALAN + Crayfish_Crayfish` = `ALAN + Crayfish` - `Crayfish`
  ) %>%
  select(iter, name, ends_with("_Control"), ends_with("_ALAN"), ends_with("_Crayfish")) %>%
  pivot_longer(
    cols = -c(iter, name),
    names_to = "contrast",
    values_to = "diff"
  )

# Bayes Factor data frame -----------------------------------------------

dat_con = data.frame(
  contrast = c(
    "ALAN + Crayfish - Crayfish",
    "ALAN + Crayfish - ALAN",
    "ALAN + Crayfish - Control",
    "Crayfish - Control",
    "ALAN - Control"),
  BF = c(0.31,0.50,1.87,11.1,6.71)
)

# Plot of the contrast ----------------------------------------------------
con_long$contrast <- gsub("_", " - ", con_long$contrast)

# Order in the plot -------------------------------------------------------

d_order <- c(
  "ALAN + Crayfish - Crayfish",
  "ALAN + Crayfish - ALAN",
  "ALAN + Crayfish - Control",
  "Crayfish - Control",
  "ALAN - Control")

con_long$contrast <- factor(con_long$contrast, levels = d_order)

####### MAP function to estimate the marginal a posteriori ##########

MAP_fun <- function(x) {
  d <- density(x, na.rm = TRUE)
  d$x[which.max(d$y)]
}

######## Contrast Ridge plot ##########
con_long |> 
  ggplot(aes(x = (diff*100), y = contrast, fill = contrast)) +
  ggridges::geom_density_ridges(scale = 1.0, alpha = 0.7) +
  #geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  stat_summary(fun = MAP_fun, geom = "point", color = "black", size = 2) +
  geom_text(
    data = dat_con,  # still use small summary frame
    aes(
      x = 10,            # fixed x position
      y = contrast,
      label = sprintf("BF = %.2f", BF)
    ),
    size = 3, hjust = -2, vjust = -0.5, fontface = "bold", color = "black"
  ) +
  scale_fill_manual(values = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"))+
  ggridges::theme_ridges() +
  labs(
    x = "Estimated difference in aquatic diet contribution (%)",
    y = "Treatment comparison")+
  theme_rsm() +
  theme(legend.position = "None",
        axis.title = element_text(size = 12, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black", face = "bold"),
        axis.title.x = element_text(margin = margin(t = 15)))

###### Save ridge plot of contrasts Figure 3 ##########
ggsave(
  "output/spider/est_diff_ridge.png",
  height = 5,
  width = 8,
  units = "in",
  dpi = 600)


contrasts <- sp_wide %>%
  mutate(
    ALAN_Control = `ALAN` - `Control`,
    Crayfish_Control = `Crayfish` - `Control`,
    ALANCray_Control = `ALAN + Crayfish` - `Control`,
    ALANCray_ALAN = `ALAN + Crayfish` - `ALAN`,
    ALANCray_Crayfish = `ALAN + Crayfish` - `Crayfish`
  ) %>%
  select(iter, name, ends_with("_Control"), ends_with("_ALAN"), ends_with("_Crayfish")) %>%
  pivot_longer(
    cols = -c(iter, name),
    names_to = "contrast",
    values_to = "diff"
  ) %>%
  group_by(name, contrast) %>%
  summarise(
    Median = median(diff),
    Mean = mean(diff),
    MAP = {
      d <- density(diff)
      d$x[which.max(d$y)]
    },
    # Lower95 = quantile(diff, probs = 0.025),
    # Upper95 = quantile(diff, probs = 0.975),
    p_g = mean(diff > 0),
    BF = p_g/(1-p_g),
    hdi_bounds = list(bayestestR::hdi(diff, ci = 0.95)),
    Lower95 = purrr::map_dbl(hdi_bounds, ~ .x$CI_low),
    Upper95 = purrr::map_dbl(hdi_bounds, ~ .x$CI_high),
    .groups = "drop"
  ) |>
  select(-hdi_bounds)

write.csv(contrasts, "col_output/spider.csv")



#### Invasive signal crayfish polygon and diet proportion plot --------------
crayfish_ <- source_crayfish_tef_corrected %>%
  dplyr::rename(Taxa = Family) %>%
  dplyr::filter(!Treatment %in% c("Control", "ALAN")) |> 
  dplyr::mutate(Treatment = factor(Treatment, levels = c( "Control","ALAN", "Crayfish", "ALAN + Crayfish"))) |> 
  dplyr::mutate(Taxa = factor(Taxa, levels = c("Chironomidae", "Gammaridae", "Simuliidae","Ephemeroptera"))) |> 
  dplyr::group_by(Treatment,Flume, Taxa) |>
  dplyr::summarise(
    d13C = mean(d13C, na.rm = TRUE),
    d15N = mean(d15N, na.rm = TRUE), .groups = "drop") |>
  ggplot(aes(x = d13C, y = d15N, color = Taxa, fill = Taxa)) +
  ggpubr::stat_chull(
    geom = "polygon",
    color = "transparent",
    fill = "#000000",
    alpha = 0.15
  ) +
  geom_point(
    size = 2,
    alpha = 1
  ) +
  geom_point(
    data = mixture_crayfish %>%
      dplyr::rename(Taxa = Family) %>%
      dplyr::mutate(Treatment = factor(Treatment, levels = c("Control","ALAN", "Crayfish", "ALAN + Crayfish"))) |> 
      dplyr::group_by(Treatment,Flume, Taxa) |>
      dplyr::summarise(
        d13C = mean(d13C, na.rm = TRUE),
        d15N = mean(d15N, na.rm = TRUE), .groups = "drop") |>
      dplyr::mutate(flume_nr = readr::parse_number(Flume))%>%
      dplyr::group_by(Treatment) %>%
      dplyr::mutate(flume_nr = as.numeric(factor(flume_nr))),
    color = "#000000",
    fill = "transparent",
    size = 2.5,
    alpha = 1
  ) +
  scale_color_manual(values = col_src_crayfish) +
  scale_fill_manual(values = col_src_crayfish) +
  ggh4x::facet_nested(Treatment ~.) +
  #ggh4x::facet_nested(ALAN + Crayfish ~ "Replicate" + flume_nr) +
  labs(
    x = "**&delta;<sup>13</sup>C (&permil;)**",
    y = "**&delta;<sup>15</sup>N (&permil;)**"
  ) +
  theme_rsm() +
  theme(
    strip.text.y = element_blank(),
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown(),
    legend.key.size = unit(0.5,"line"),
    legend.title.position = "top",
    legend.direction = "horizontal",
    legend.title = element_text(hjust = 0.5),
    legend.position = "top"
  )

####### Save crayfish isotope polygon plot -------------
ggsave(
  "output/crayfish/crayfish_.png",
  crayfish_,
  height = 12,
  width = 16.5,
  units = "cm"
)

##### MixSIAR --------
# only select flumes with crayfish
# Moderate confidence (numbers scaled to 10)
# set.seed(123)
# crayfish_models <-
#   unique(mixture_crayfish$Treatment[mixture_crayfish$Crayfish == "Crayfish"]) %>%
#   purrr::map(function(x) {
#     run_crayfish_model(
#       run = "long",
#       treatment = x,
#       mix = mixture_crayfish,
#       src = source_crayfish,
#       discr = TEF_crayfish
#     )
#   })
##########Read in the saved pre-run saved crayfish models ############

crayfish_models = readRDS("rds_data/crayfish/crayfish_models.rds")

##### Crayfish model trace plots for model diagnostics --------------

crayfish_models %>%
  purrr::walk(
    function(x) {
      plt_temp <-
        x$draws %>%
        dplyr::mutate(iter = 1:n()) %>%
        tidyr::pivot_longer(-c(Treatment, iter)) %>%
        ggplot(aes(iter, value, color = name)) +
        geom_line(linewidth = 0.1) +
        scale_color_manual(values = col_src_crayfish) +
        facet_grid(name~Treatment) +
        theme_rsm()
      
      ggsave(
        paste0("output/crayfish_trace/crayfish_", unique(x$draws$Treatment), "_trace.png"),
        plt_temp,
        height = 10,
        width = 16.5,
        units = "cm"
      )
    }
  )


##### Crayfish dietary proportion plot ---------------------
plt_diet_crayfish_full <-
  dplyr::bind_rows(
    crayfish_models %>%
      purrr::map_df(~.x$draws %>% dplyr::mutate(iter = 1:n()))
  ) %>%
  dplyr::mutate(`Chiro&Gamm` = Chironomidae+Gammaridae) %>%
  dplyr::select(-c(Chironomidae, Gammaridae)) %>%
  tidyr::pivot_longer(-c(Treatment, iter)) %>%
  tidyr::drop_na(value) %>%
  dplyr::left_join(
    mixture_crayfish %>%
      dplyr::select(Treatment, ALAN, Crayfish) %>%
      dplyr::distinct(),
    by = "Treatment"
  ) %>%
  dplyr::mutate(Treatment = factor(Treatment, 
                                   levels = c("Crayfish", "ALAN + Crayfish"))) |> 
  dplyr::group_by(Treatment, name, iter) %>%
  dplyr::summarise(value = mean(value), .groups = "drop") %>%
  dplyr::mutate(Treat = case_when(Treatment == "Crayfish" ~ "No ALAN",
                                  Treatment == "ALAN + Crayfish" ~ "ALAN",
                                  TRUE ~ NA_character_)) |> 
  dplyr::mutate(Treat = factor(Treat, 
                               levels = c("No ALAN", "ALAN"))) |>
  ggplot(aes(x = value, fill = name, color = name)) +
  geom_density(alpha = 0.8) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_color_manual(values = unname(col_src_crayfish)[c(2,4,1)], guide = "legend", name = "Taxa",
                     limits = c("Chiro&Gamm","Simuliidae","Ephemeroptera")) +
  scale_fill_manual(values = unname(col_src_crayfish)[c(2,4,1)], guide = "legend", name = "Taxa",
                    limits = c("Chiro&Gamm","Simuliidae","Ephemeroptera")) +
  scale_x_continuous(breaks = 0:4/4, labels = scales::percent_format()) +
  ggh4x::facet_nested(Treat~1) +
  labs(y = NULL, x = expression(bold("Dietary proportion"))) +
  theme_rsm() +
  theme(
    strip.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 33),
    axis.ticks.y = element_blank(),
    legend.key.size = unit(0.5,"line"),
    legend.title.position = "top",
    legend.direction = "horizontal",
    legend.title = element_text(hjust = 0.5),
    legend.position = "top"
  )

####### Save combined crayfish isotope polygon and diet proportion plot -------------
######The cr_p is the crayfish isotopic niche width plot saved separately ##########
######The cr_p code is in the crayfish isotopic niche width script ##########

ggsave(
  "output/crayfish/crayfish_poly_prop_niche.png",
  (crayfish_ | plt_diet_crayfish_full) / (cr_p) +
    patchwork::plot_annotation(
      tag_levels = 'a',
      theme = theme(
        plot.tag = element_text(face = "bold", size = 14)
      )) +
    theme(
      legend.position = "top",
      legend.title = element_text(face = "bold", size = 14, family = "sans"),
      legend.text = element_text(margin = margin(l = 2), size = 8, face = "bold"),
      legend.margin = margin(),
      legend.justification = "center",
      legend.box.spacing = unit(0.3, "cm")),
  height = 10,
  width = 9,
  units = "in",
  dpi = 600
)

#####To estimate the median dietary contribution frome each source#####
cr_sum = dplyr::bind_rows(
  crayfish_models %>%
    purrr::map_df(~.x$draws %>% dplyr::mutate(iter = 1:n()))) %>%
  dplyr::mutate(ChiroGamma = Chironomidae+Gammaridae) %>%
  dplyr::select(-c(Chironomidae, Gammaridae)) %>%
  tidyr::pivot_longer(-c(Treatment, iter)) %>%
  dplyr::left_join(
    mixture_spider %>%
      dplyr::select(Treatment,ALAN, Crayfish) %>%
      dplyr::distinct(),
    by = "Treatment"
  ) %>%
  dplyr::group_by(Treatment, name) %>%
  dplyr::summarise(
    Median = median(value),
    MAP = value[which.max(density(value)$y)], 
    hdi_bounds = list(hdi(value, ci = 0.95)),
    Lower95 = map_dbl(hdi_bounds, ~ .x$CI_low),
    Upper95 = map_dbl(hdi_bounds, ~ .x$CI_high),
    .groups = "drop") |> 
  select(-hdi_bounds)
