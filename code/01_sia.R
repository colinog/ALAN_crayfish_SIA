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
    )
  )

##### Polgyon Plot -----------------------
# plt_polygon_spider <- 
# source_spider_tef_corrected %>%
#   ggplot(aes(x = d13C, y = d15N, color = Family, fill = Family)) +
#   ggpubr::stat_chull(
#     geom = "polygon",
#     color = "transparent",
#     fill = "#000000",
#     alpha = 0.15
#   ) +
#   geom_point(
#     size = 1.5,
#     alpha = 0.6
#     ) +
#   # For mean with sd (not really valuable if plotted for each stream)
#   # ggstance::geom_errorbarh(
#   #   data = sstcs,
#   #   mapping = aes(xmin = d13C_mean - d13C_sd, 
#   #                 xmax = d13C_mean + d13C_sd,
#   #                 y = d15N_mean,
#   #                 color = Family),
#   #   inherit.aes = FALSE
#   # ) +
#   # geom_errorbar(
#   #   data = sstcs,
#   #   mapping = aes(ymin = d15N_mean - d15N_sd, 
#   #                 ymax = d15N_mean + d15N_sd,
#   #                 x = d13C_mean,
#   #                 color = Family),
#   #   inherit.aes = FALSE
#   # ) +
#   # geom_point(
#   #   data = sstcs,
#   #   mapping = aes(x = d13C_mean, y = d15N_mean),
#   #   size = 2
#   # ) +
#   geom_point(
#     data = mixture_spider %>%
#       dplyr::mutate(flume_nr = readr::parse_number(Flume)) %>%
#       dplyr::group_by(Treatment) %>%
#       dplyr::mutate(flume_nr = as.numeric(factor(flume_nr))),
#     color = "#000000",
#     fill = "transparent",
#     size = 2,
#     alpha = 0.6
#   ) +
#   scale_color_manual(values = col_src_spider, breaks = names(col_src_spider)) +
#   scale_fill_manual(values = col_src_spider, breaks = names(col_src_spider)) +
#   ggh4x::facet_nested(ALAN + Crayfish~"Replicate" + flume_nr) +
#   labs(
#     strip.text.y = element_blank(),
#     x = "&delta;<sup>13</sup>C (&permil;)",
#     y = "&delta;<sup>15</sup>N (&permil;)"
#   ) +
#   theme_rsm() +
#   theme(
#     strip.text.y = element_blank(),
#     axis.title.x = ggtext::element_markdown(),
#     axis.title.y = ggtext::element_markdown()
#   )
##############The above will not be used#######


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
  "output/spider.png",
  plt_tet,
  height = 12,
  width = 16.5,
  units = "cm"
)

##### MixSIAR --------
spider_models <-
  unique(mixture_spider$Treatment) %>%
  purrr::map(function(x) {
    run_spider_model(
      run = "long", 
      treatment = x,
      mix = mixture_spider,
      src = source_spider,
      discr = TEF_spider
    )
  })

saveRDS(spider_models, "rds/spider_model.rds")

spider_models = readRDS("rds/spider_model.rds")


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
        paste0("output/spider_", unique(x$draws$Treatment), "_trace.png"),
        plt_temp,
        height = 10,
        width = 16.5,
        units = "cm"
      )
    }
  )


#"legend", name = "Habitat",
#limits = c("Terrestrial", "Aquatic")
##### Diet Plot ---------------------
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
  dplyr::summarise(value = mean(value)) %>%
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

design <- "AAABB"

# plt_diet_spider <- plt_diet_spider +
#   theme(
#     legend.key.width = unit(0.8, "lines"),
#     legend.key.height = unit(0.5, "lines"),
#     legend.text = element_text(size = 8),   # optional
#     legend.title = element_text(size = 9)   # optional
#     )
# 
# library(gridExtra)
# 
# gridExtra::grid.arrange(plt_tet, plt_diet_spider, ncol = 2)

ggsave(
  "output/spider_full.png",
  #plt_polygon_spider 
  plt_tet + 
    plt_diet_spider + 
    patchwork::plot_annotation(
      tag_levels = "a") +
    theme(plot.tag.position = c(0, 1),       # x = 0 (left), y = 1 (top)
           plot.tag = element_text(hjust = 0.5, vjust = 1)) +
    patchwork::plot_layout(design = design, guides = "keep") &
    theme(legend.position = "top",
          legend.box = "horizontal",
          #legend.box.just = "left",
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
  

combined <- (plt_tet + 
               plt_diet_spider) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "top",
    legend.box.just = "left",
    legend.justification = "left",
    legend.box.spacing = unit(0.2, "cm")
  )

combined



ggsave(
  "output/spider_full.png",
  #plt_polygon_spider 
  plt_tet + 
    plt_diet_spider + 
    patchwork::plot_layout(design = design, guides = "collect") & 
    theme(legend.position = "top", 
          legend.margin = margin(), 
          legend.box.spacing = unit(0, "cm"),
          legend.justification = "left",
          legend.text = element_text(margin = margin(l = -5)),
          legend.key.width = unit(1.7, "lines"),
          legend.key.height = unit(1, "lines"),
          legend.key.size = unit(0.2, "cm")),
  height = 12,
  width = 16.5,
  units = "cm",
  dpi = 600
)


#####To estimate the median dietary contribution frome each source#####
#The estimation of the additive effect minus the combined effect:
##Combined effect (Proportion) - (ALAN + Crayfish individual)
##"≈ 0 → Additive, > 0 → Synergistic, < 0 → Antagonistic"

sp_wide <- sp_sum %>%
  select(iter, Treatment, name, value) %>%
  tidyr::pivot_wider(names_from = Treatment, values_from = value)


con_long = sp_wide %>%
  dplyr::filter(!name %in% "Terrestrial") |>  # <-- Filter if needed
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
  )

MAP_fun <- function(x) {
  d <- density(x, na.rm = TRUE)
  d$x[which.max(d$y)]
}

dat_con = data.frame(
  contrast = c(
    "ALANCray - Crayfish",
    "ALANCray - ALAN",
    "ALANCray - Control",
    "Crayfish - Control",
    "ALAN - Control"),
  BF = c(0.27,0.51,1.90,11.99,6.85)
)

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
    x = "Estimated difference in dietary proportion (%)",
    y = "Treatment comparison")+
  theme_rsm() +
  theme(legend.position = "None",
        axis.title = element_text(size = 12, color = "black", face = "bold"),
        axis.text = element_text(size = 11, color = "black", face = "bold"),
        axis.title.x = element_text(margin = margin(t = 15)))

ggsave(
  "output/est_diff_ridge.png",
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
    Lower95 = quantile(diff, probs = 0.025),
    Upper95 = quantile(diff, probs = 0.975),
    p_g = mean(diff > 0),
    BF = p_g/(1-p_g),
    # hdi_bounds = list(bayestestR::hdi(diff, ci = 0.95)),
    # Lower95 = purrr::map_dbl(hdi_bounds, ~ .x$CI_low),
    # Upper95 = purrr::map_dbl(hdi_bounds, ~ .x$CI_high),
    .groups = "drop"
  ) 

write.csv(contrasts, "col_output/spider.csv")


# Plot of the contrast ----------------------------------------------------
con_long$contrast <- gsub("_", " - ", con_long$contrast)


# Order in the plot -------------------------------------------------------

d_order <- c(
  "ALANCray - Crayfish",
  "ALANCray - ALAN",
  "ALANCray - Control",
  "Crayfish - Control",
  "ALAN - Control")

con_long$contrast <- factor(con_long$contrast, levels = d_order)


#################
contrasts |>
  dplyr::filter(!name %in% "Terrestrial") |>  # <-- Filter if needed
  # mutate(
  #   name = factor(name, levels = c("Aquatic", "Terrestrial"))
  # ) |> 
  ggplot(aes(x = contrast, y = (MAP*100))) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = (Lower95*100), ymax = (Upper95*100)), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  coord_flip() +
  facet_wrap(~ name, scales = "free_y") +  # <-- Facet by habitat
  labs(
    y = "Estimated MAP difference (%)",
    x = "Treatment differences"
  ) +
  theme_rsm()+
  theme(
    #strip.text.y = element_text(angle = 0),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"), 
    axis.title.y = element_text(),
    axis.title = element_text(size = 11, color = "black", face = "bold"),
    axis.text = element_text(size = 11, color = "black", face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(size = 11, color = "black", face = "bold"),
    legend.title = element_text(size = 11, face = "bold"),
    plot.tag = element_text(size = 11, face = "bold"),
    strip.text = element_text(face = "bold", size = 11),
    plot.margin = margin(t = 10, r = 20, b = 10, l = 10), 
    strip.background = element_rect(fill = "grey90", color = "black")
  )

#################Save the figure#####
ggsave(
  "output/est_MAP_diff.png",
  height = 12,
  width = 20,
  units = "cm",
  dpi = 600
)


########################




p_diff <- sp_wide %>%
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
    p_g = mean(diff > 0),
    BF = p_g/(1-p_g),
    .groups = "drop"
  )

write.csv(p_diff, "col_output/spider_p_diff.csv")

sp_sum = dplyr::bind_rows(
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
  ) %>%
  dplyr::group_by(Treatment, name) %>%
  dplyr::summarise(
    Median = median(value),
    Mean = mean(value),
    MAP = value[which.max(density(value)$y)], 
    Lower95 = quantile(value,probs = 0.025),
    Upper95 = quantile(value,probs = 0.975),
    .groups = "drop"
    # hdi_l = list(bayestestR::hdi(value, ci = 0.95)),
    # Lower95 = purrr::map_dbl(hdi_l, ~ .x$CI_low),
    # Upper95 = purrr::map_dbl(hdi_l, ~ .x$CI_high)
    ) 

#library(bayestestR)

# Save as tab-delimited
write.table(sp_sum, file = "col_output/tet_spsum_new.txt", sep = "\t", row.names = FALSE)


# cray_sc = source_crayfish |> 
#   dplyr::filter(!Treatment %in% c("ALAN", "Control"))
# 
# write.csv(cray_sc, "data/cray_sc.csv", row.names = FALSE)

#### Crayfish --------------
source_crayfish_tef_corrected <- 
  source_crayfish %>%
  # adjust soruces by TEF
  dplyr::left_join(TEF_crayfish, by = dplyr::join_by(Family)) %>%
  dplyr::mutate(d15N = d15N + Meand15N,
                d13C = d13C + Meand13C,
                flume_nr = readr::parse_number(Flume)) %>%
  dplyr::group_by(Treatment) %>%
  dplyr::mutate(flume_nr = as.numeric(factor(flume_nr)))

sctcs <- # source_crayfish_tef_corrected_summary
  source_crayfish_tef_corrected %>%
  dplyr::group_by(ALAN, Crayfish, Family, flume_nr) %>%
  dplyr::summarise(
    dplyr::across(
      c(d15N, d13C),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd = ~sd(.x, na.rm = TRUE) 
      )
    )
  )

##### Polgyon Plot -----------------------
# plt_polygon_crayfish <- 
#   source_crayfish_tef_corrected %>%
#   ggplot(aes(x = d13C, y = d15N, color = Family, fill = Family)) +
#   ggpubr::stat_chull(
#     geom = "polygon",
#     color = "transparent",
#     fill = "#000000",
#     alpha = 0.15
#   ) +
#   geom_point(
#     size = 1.5,
#     alpha = 0.6
#   ) +
#   geom_point(
#     data = mixture_crayfish %>%
#       dplyr::mutate(flume_nr = readr::parse_number(Flume)) %>%
#       dplyr::group_by(Treatment) %>%
#       dplyr::mutate(flume_nr = as.numeric(factor(flume_nr))),
#     color = "#000000",
#     fill = "transparent",
#     size = 2,
#     alpha = 0.6
#   ) +
#   scale_color_manual(values = col_src_crayfish) +
#   scale_fill_manual(values = col_src_crayfish) +
#   ggh4x::facet_nested(ALAN + Crayfish~"Replicate" + flume_nr) +
#   labs(
#     strip.text.y = element_blank(),
#     x = "&delta;<sup>13</sup>C (&permil;)",
#     y = "&delta;<sup>15</sup>N (&permil;)"
#   ) +
#   theme_rsm() +
#   theme(
#     strip.text.y = element_blank(),
#     axis.title.x = ggtext::element_markdown(),
#     axis.title.y = ggtext::element_markdown()
#   )

########New plot with mean ######
crayfish_ <- source_crayfish_tef_corrected %>%
  dplyr::rename(Taxa = Family) %>%
  dplyr::filter(!Treatment %in% c("Control", "ALAN")) |> 
  dplyr::mutate(Treatment = factor(Treatment, levels = c( "Control","ALAN", "Crayfish", "ALAN + Crayfish"))) |> 
  dplyr::mutate(Taxa = factor(Taxa, levels = c("Ephemeroptera", "Simuliidae", "Chironomidae", "Gammaridae"))) |> 
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

ggsave(
  "output/crayfish_.png",
  crayfish_,
  height = 12,
  width = 16.5,
  units = "cm"
)

##### MixSIAR --------
# only select flumes with crayfish
# Moderate confidence (numbers scaled to 10)

crayfish_models <-
  unique(mixture_crayfish$Treatment[mixture_crayfish$Crayfish == "Crayfish"]) %>%
  purrr::map(function(x) {
    run_crayfish_model(
      run = "long",
      treatment = x,
      mix = mixture_crayfish,
      src = source_crayfish,
      discr = TEF_crayfish
    )
  })

saveRDS(crayfish_models, "rds/new_crayfish_model_full.rds")
#saveRDS(crayfish_models, "rds/new_crayfish_model_cl.rds")
#crayfish_models_sex = readRDS("rds/crayfish_model_sex.rds")



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
        paste0("output/crayfish_cl", unique(x$draws$Treatment), "_trace.png"),
        plt_temp,
        height = 10,
        width = 16.5,
        units = "cm"
      )
    }
  )


##### Diet Plot ---------------------
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
  dplyr::summarise(value = mean(value)) %>%
  dplyr::mutate(Treat = case_when(Treatment == "Crayfish" ~ "No ALAN",
                               Treatment == "ALAN + Crayfish" ~ "ALAN",
                               TRUE ~ NA_character_)) |> 
  dplyr::mutate(Treat = factor(Treat, 
                                   levels = c("No ALAN", "ALAN"))) |>
  # # dplyr::bind_rows(
  # #   tibble::tibble(
  # #     ALAN = c("No ALAN", "ALAN"),
  # #     Crayfish = "No Crayfish",
  # #     name = "ChiroGamma",
  # #     value = NA
  # #   )
  # ) %>%
  ggplot(aes(x = value, fill = name, color = name)) +
  geom_density(alpha = 0.8) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_color_manual(values = unname(col_src_crayfish)[c(1,4,2,3)], guide = "legend", name = "Taxa",
                     limits = c("Ephemeroptera", "Simuliidae", "Chiro&Gamm")) +
  scale_fill_manual(values = unname(col_src_crayfish)[c(1,4,2,3)], guide = "legend", name = "Taxa",
                    limits = c("Ephemeroptera", "Simuliidae", "Chiro&Gamm")) +
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

plt_list = list(crayfish_ + plt_diet_crayfish_full + cr_p )

design_c <- "AAABBCC"

plt_diet_crayfish <- plt_diet_crayfish_full + theme(
  strip.text.y = element_blank()
)

wrap_plots(plt_list, design = design_c, guides = "keep") +
  plot_annotation(
    tag_levels = list(c("(a)", "(b)", "(c)")),
    theme = theme(
      plot.tag.position = c(0, 1),        # normalized top-left
      plot.tag = element_text(face = "bold", size = 10)
    )
  ) & 
  theme(
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 9, family = "Helvetica"),
    legend.text = element_text(margin = margin(l = 2), size = 8, family = "Helvetica", face = "bold"),
    legend.margin = margin(),
    legend.justification = "center",
    legend.box.spacing = unit(0.3, "cm")
  )

design_t = "
AABB
CCDD
"

ggsave(
  "output/crayfish_full.png",
  (crayfish_ | plt_diet_crayfish_full) / (cr_p) +
  patchwork::plot_annotation(
    tag_levels = 'a',
    #tag_levels = list(c("(a)", "(b)", "(c)")),
    theme = theme(
      plot.tag = element_text(face = "bold", size = 14)
    )) +
  #patchwork::plot_layout(ncol = 1, nrow = 2, guides = "keep") & 
  theme(
    #plot.tag = element_text(face = "bold", size = 14),
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 14, family = "sans"),
    legend.text = element_text(margin = margin(l = 2), size = 8, face = "bold"),
    legend.margin = margin(),
    legend.justification = "center",
    legend.box.spacing = unit(0.3, "cm")),
  height = 10,
  width = 9,
  units = "in",
  dpi = 300
)

dev.off()



ggsave(
  "output/crayfish_full.png",
  (crayfish_ | plt_diet_crayfish_full) / cr_p +
    patchwork::plot_annotation(tag_levels = list(c("(a)", "(b)", "(c)"))) +
    theme(plot.tag.position = c(0, 1),       # x = 0 (left), y = 1 (top)
           plot.tag = element_text(hjust = 0.1, vjust = 0.5))+
    #patchwork::plot_layout(design = design_c, guides = "keep") & 
    theme(legend.position = "top",
          #legend.title = element_text(face = "bold", size = 9, family = "Helvetica"),
          #legend.text = element_text(margin = margin(l = 2),size = 8, family = "Helvetica", face = "bold"),
          #legend.margin = margin(),
          legend.justification = "center",
          legend.box.spacing = unit(0.3, "cm"),
          plot.tag.position =  "topleft",
          plot.tag = element_text(face = "bold",size = 10)
          ),
  height = 10,
  width = 9,
  units = "in",
  dpi = 600
)



str(crayfish_models, max.level = 2)

library(bayestestR)

cr_sum = dplyr::bind_rows(
  crayfish_models %>%
    purrr::map_df(~.x$draws %>% dplyr::mutate(iter = 1:n()))) %>%
  dplyr::mutate(ChiroGamma = Chironomidae+Gammaridae) %>%
  dplyr::select(-c(Chironomidae, Gammaridae)) %>%
  #dplyr::select(Treatment:Terrestrial) %>%
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
    # hdi_bounds = list(hdi(value, ci = 0.95)),
    # Lower95 = map_dbl(hdi_bounds, ~ .x$CI_low),
    # Upper95 = map_dbl(hdi_bounds, ~ .x$CI_high),
    Lower95 = quantile(value, 0.025),
    Upper95 = quantile(value, 0.975),
    .groups = "drop")


write.table(cr_sum, file = "col_output/cr_new_clsum.txt", sep = "\t", row.names = FALSE)


# ggplot(cr_sum, aes(x = name, y = m_value, fill = interaction(Treatment))) +
#   geom_bar(stat = "identity", position = "dodge") +
#   geom_errorbar(aes(ymin = hdi_lower, ymax = hdi_upper), width = 0.2, position = position_dodge(0.9)) +
#   labs(y = "Mean Dietary Proportion", x = "Source") +
#   theme_minimal()


cl_g = ggplot(sp_sum, aes(x = name, y = m_value, color = Treatment)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = hdi_lower, ymax = hdi_upper),
                width = 0.1,
                position = position_dodge(width = 0.5)) +
  ylab("Proportional Contribution") +
  theme_rsm()

sex_g
cl_g





source_cray = source_crayfish %>%
  dplyr::filter(Treatment %in% c("Crayfish", "ALAN + Crayfish")) |> 
  dplyr::group_by(Family,Treatment, Flume) %>%
  dplyr::summarize(n = n(),
            Meand15N = mean(d15N, na.rm = T),
            SDd15N = sd(d15N, na.rm = T),
            Meand13C = mean(d13C, na.rm = T),
            SDd13C = sd(d13C, na.rm = T), .groups = "drop")

write.csv(source_cray, "data/source_cray.csv", row.names = F)

mix_cr <- load_mix_data(filename = "data/mixture_crayfish.csv", 
                     iso_names = c("d13C","d15N"), 
                     factors   = c("Flume","Treatment"), 
                     fac_random = c(TRUE,FALSE), 
                     fac_nested = c(FALSE,FALSE), 
                     cont_effects = "cp_mm")

source_cr <- load_source_data(filename = "data/source_cray.csv",
                           source_factors = NULL, 
                           conc_dep = FALSE, 
                           data_type = "means", 
                           mix_cr)


discr <- load_discr_data(filename="data/TEF_crayfish.csv", mix_cr)

# Write the JAGS model file
model_filename <- "MixSIAR_model_test.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix_cr, source_cr)


jags.1 <- run_model(run="test", mix_cr, source_cr, discr, model_filename, resid_err=T, process_err=T)

dev.off()
plot_prior(alpha.prior = c(1, 0.5,0.5,0.5), source_cr, plot_save_pdf=FALSE)

jags_1 = "MixSIAR_model_crayfish.txt"

names(jags.1)

ab = jags.1$BUGSoutput$sims.matrix |>
  as.data.frame() |> 
  dplyr::select(dplyr::starts_with("ilr.cont1"))


ab$`ilr.cont1[1]` +
  ggplot2::theme(legend.position="right")




aaa = plot_continuous_var(jags.1, mix_cr, source_cr, output_options, alphaCI = 0.05, 
                    exclude_sources_below = 0.1)

output_options <- list(summary_save = TRUE,
                       summary_name = "summary_statistics",
                       sup_post = FALSE,
                       plot_post_save_pdf = TRUE,
                       plot_post_name = "posterior_density",
                       sup_pairs = FALSE,
                       plot_pairs_save_pdf = TRUE,
                       plot_pairs_name = "pairs_plot",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "diagnostics",
                       indiv_effect = FALSE,
                       plot_post_save_png = FALSE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)
