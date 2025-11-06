run_crayfish_mod <- 
  function(
    run = "short",
    treatment,
    mix,
    src,
    discr
  ) {
    cat(paste0("##### ", Sys.time(), ": Initializing model for ", treatment, " #####\n\n"))
    
    mix %>%
      dplyr::filter(Treatment == treatment) %>%
      write.csv("data/temp/mixture_crayfish_filtered.csv")
    
    mix_crayfish <- load_mix_data("data/temp/mixture_crayfish_filtered.csv",
                                  iso_names = c("d13C", "d15N"),
                                  factors = NULL,
                                  fac_random = NULL,
                                  fac_nested = NULL,
                                  cont_effects = NULL)
    
    src_temp <- 
      src %>%
      dplyr::filter(Treatment == treatment) 
    
    write.csv(src_temp, "data/temp/source_crayfish_means.csv", row.names = src_temp$Family)
    
    src_crayfish <- load_source_data("data/temp/source_crayfish_means.csv",
                                     source_factors = NULL,
                                     conc_dep = FALSE,
                                     data_type = "raw",
                                     mix_crayfish)
    
    TEF_crayfish %>%
      dplyr::filter(Family %in% src_crayfish$source_names) %>%
      dplyr::select(-Family) %>%
      write.csv("data/temp/TEF_crayfish_filtered.csv")
    
    discr_crayfish <- load_discr_data("data/temp/TEF_crayfish_filtered.csv",
                                      mix_crayfish)
    
    write_JAGS_model(
      filename = "MixSIAR_model_crayfish.txt",
      resid_err = TRUE,
      process_err = TRUE,
      mix = mix_crayfish, 
      source = src_crayfish
    )
    
    crayfish_model <-
      run_model(
        run = run,
        mix = mix_crayfish, source = src_crayfish, discr = discr_crayfish,
        model_filename = "MixSIAR_model_crayfish.txt",
        process_err = TRUE,
        resid_err = TRUE
      )
    
    crayfish_output <- 
      crayfish_model$BUGSoutput$sims.matrix %>%
      as.data.frame() %>%
      dplyr::select(dplyr::starts_with("p.global"))
    
    colnames(crayfish_output) <- src_crayfish$source_names
    
    return(
      list(
        model = crayfish_model,
        draws = crayfish_output %>% dplyr::mutate(Treatment = treatment)
      )
    )
  }


cray_mod1 <-
  #unique(mixture_crayfish$Treatment) %>%
  unique(mixture_crayfish$Treatment[mixture_crayfish$Crayfish == "Crayfish"]) |> 
  purrr::map(function(x) {
    run_crayfish_mod(
      run = "test",
      treatment = x,
      mix = mixture_crayfish,
      src = source_crayfish,
      discr = TEF_crayfish
    )
  })


#####Trace plot crayfish model######
plt_diet_cr <-
  dplyr::bind_rows(
    cray_mod %>%
      purrr::map_df(~.x$draws %>% dplyr::mutate(iter = 1:n()))
  ) %>%
  dplyr::mutate(ChiroGamma = Chironomidae+Gammaridae) %>%
  dplyr::select(-c(Chironomidae, Gammaridae)) %>%
  tidyr::pivot_longer(-c(Treatment, iter)) %>%
  tidyr::drop_na(value) %>%
  dplyr::left_join(
    mixture_crayfish %>%
      dplyr::select(Treatment, ALAN, Crayfish) %>%
      dplyr::distinct(),
    by = "Treatment"
  ) %>%
  dplyr::group_by(ALAN, Crayfish, name, iter) %>%
  dplyr::summarise(value = mean(value), .groups = "drop") %>%
  dplyr::bind_rows(
    tibble::tibble(
      ALAN = c("No ALAN", "ALAN"),
      Crayfish = "No Crayfish",
      name = "ChiroGamma",
      value = NA
    )
  ) %>%
  ggplot(aes(x = value, fill = name, color = name)) +
  geom_density(alpha = 0.5) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0,NA)) +
  scale_color_manual(values = unname(col_src_crayfish)[c(3,1,4)], guide = "none") +
  scale_fill_manual(values = unname(col_src_crayfish)[c(2,1,4)], guide = "none") +
  scale_x_continuous(breaks = 0:4/4, labels = scales::percent_format()) +
  ggh4x::facet_nested(ALAN+Crayfish~1) +
  labs(y = NULL, x = "Dietary proportion") +
  theme_rsm() +
  theme(
    strip.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 33),
    axis.ticks.y = element_blank()
  )

design <- "AAAAB"

ggsave(
  "output/crayfish_t.png",
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

crayfish_models

cr_sum = dplyr::bind_rows(
  crayfish_models %>%
    purrr::map_df(~.x$draws %>% dplyr::mutate(iter = 1:n()))) %>%
  dplyr::mutate(ChiroGamma = Chironomidae+Gammaridae) %>%
  dplyr::select(-c(Chironomidae, Gammaridae)) %>%
  #dplyr::select(Treatment:Terrestrial) %>%
  tidyr::pivot_longer(-c(Flume, iter)) %>%
  dplyr::left_join(
    mixture_spider %>%
      dplyr::select(Flume,ALAN, Crayfish) %>%
      dplyr::distinct(),
    by = "Flume"
  ) %>%
  dplyr::group_by(ALAN, name) %>%
  dplyr::summarise(m_value = median(value),
                   min_v = min(value),
                   max_v = max(value))
