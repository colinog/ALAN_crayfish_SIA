# You could also create a lockfile from your old R library path, 
# and then restore based on that. 
# For example, something like this could work:
# renv::snapshot(
#   library  = "~/r/win-library/3.6",
#   lockfile = "r-36.lock",
#   type     = "simple"
# )
# 
# renv::restore(
#   lockfile = "r-36.lock"
# )

# spider_mod <-
#   unique(mixture_spider$Flume) %>%
#   purrr::map(function(x) {
#     run_spider_model(
#       run = "short",
#       flume = x,
#       mix = mixture_spider,
#       src = source_spider,
#       discr = TEF_spider
#     )
#   })




plt_tet = source_spider_tef_corrected %>%
  dplyr::mutate(Treatment = factor(Treatment, levels = c("ALAN + Crayfish", "ALAN", "Crayfish", "Control"))) |> 
  dplyr::group_by(Treatment,Flume, Family) |>
  dplyr::summarise(
    d13C = mean(d13C, na.rm = TRUE),
    d15N = mean(d15N, na.rm = TRUE), .groups = "drop"
  ) |>
  ggplot(aes(x = d13C, y = d15N, color = Family, fill = Family)) +
  ggpubr::stat_chull(
    geom = "polygon",
    color = "transparent",
    fill = "#000000",
    alpha = 0.15
  ) +
  geom_point(
    size = 1.5,
    alpha = 0.6
  ) +
  geom_point(
    data = mixture_spider %>%
      dplyr::mutate(Treatment = factor(Treatment, levels = c("ALAN + Crayfish", "ALAN", "Crayfish", "Control"))) |> 
      dplyr::group_by(Treatment,Flume, Family) |> 
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
  ggh4x::facet_nested("Treatment") +
  #ggh4x::facet_nested(ALAN + Crayfish~"Replicate" ) +
  labs(
    strip.text.y = element_blank(),
    x = "&delta;<sup>13</sup>C (&permil;)",
    y = "&delta;<sup>15</sup>N (&permil;)"
  ) +
  theme_rsm() +
  theme(
    strip.text.y = element_blank(),
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown()
  )

ggsave(
  "output/spider_.png",
  #plt_polygon_spider 
  plt_tet + plt_diet_tet + 
    patchwork::plot_layout(design = design, guides = "collect") & 
    theme(legend.position = "top", 
          legend.margin = margin(), 
          legend.box.spacing = unit(0, "cm")),
  height = 12,
  width = 16.5,
  units = "cm",
  dpi = 600
)

#######Crayfish ######
crayfish_ <- source_crayfish_tef_corrected %>%
  dplyr::mutate(Treatment = factor(Treatment, levels = c("ALAN + Crayfish", "ALAN", "Crayfish", "Control"))) |> 
  dplyr::group_by(Treatment,Flume, Family) |>
  dplyr::summarise(
    d13C = mean(d13C, na.rm = TRUE),
    d15N = mean(d15N, na.rm = TRUE), .groups = "drop") |>
  ggplot(aes(x = d13C, y = d15N, color = Family, fill = Family)) +
  ggpubr::stat_chull(
    geom = "polygon",
    color = "transparent",
    fill = "#000000",
    alpha = 0.15
  ) +
  geom_point(
    size = 1.5,
    alpha = 0.6
  ) +
  geom_point(
    data = mixture_crayfish %>%
      dplyr::mutate(Treatment = factor(Treatment, levels = c("ALAN + Crayfish", "ALAN", "Crayfish", "Control"))) |> 
      dplyr::group_by(Treatment,Flume, Family) |>
      dplyr::summarise(
        d13C = mean(d13C, na.rm = TRUE),
        d15N = mean(d15N, na.rm = TRUE), .groups = "drop") |>
      dplyr::mutate(flume_nr = readr::parse_number(Flume))%>%
      dplyr::group_by(Treatment) %>%
      dplyr::mutate(flume_nr = as.numeric(factor(flume_nr))),
    color = "#000000",
    fill = "transparent",
    size = 2,
    alpha = 0.6
  ) +
  scale_color_manual(values = col_src_crayfish) +
  scale_fill_manual(values = col_src_crayfish) +
  ggh4x::facet_nested(Treatment ~.) +
  #ggh4x::facet_nested(ALAN + Crayfish ~ "Replicate" + flume_nr) +
  labs(
    strip.text.y = element_blank(),
    x = "&delta;<sup>13</sup>C (&permil;)",
    y = "&delta;<sup>15</sup>N (&permil;)"
  ) +
  theme_rsm() +
  theme(
    strip.text.y = element_blank(),
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown()
  )


combined_plot =  crayfish_+ plt_diet_crayfish +
  patchwork::plot_layout(design = design, guides = "collect")
  
combined_plot = combined_plot & theme(legend.position = "top",
                                      legend.margin = margin(), 
                                      legend.box.spacing = unit(0, "cm"))

ggsave("output/crayfish_.png", height = 12, width = 16.5,units = "cm",dpi = 600)
  

ggsave(
  "output/crayfish_.png",
  crayfish_ + plt_diet_crayfish + 
    patchwork::plot_layout(design = design, guides = "collect") & 
    theme(legend.position = "top", 
          legend.margin = margin(), 
          legend.box.spacing = unit(0, "cm")),
  height = 12,
  width = 16.5,
  units = "cm",
  dpi = 600
)  
  