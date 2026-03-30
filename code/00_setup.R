################################################################
# Setup script for the manuscript:
# Artificial light at night and invasive signal crayfish alter aquatic-terrestrial food webs
################################################################

# Load libraries required ----------------------------------------------------------
library(tidyverse)
library(Hmisc)
library(ggstance)
library(MixSIAR)
library(SIBER)
library(ggthemes)
library(ggh4x)
library(ggtext)
library(patchwork)
library(bayestestR)
library(ggplotify)


########Define colors for sources in spider polygon plot --------------------------
 col_src_spider <- c(
   Chironomidae = "#4f6a8f",
   Tanypodinae = "#3197bc",
   Cicadellidae = "#3E9C94",
   Curculionidae = "#7A9E24")
 
########Define colors and themes for spider and crayfish diet plot -------------------
col_hab = c(
  Chironomidae = "#56B1F7",
  Tanypodinae = "#56B1F7",
  Cicadellidae = "#66A61E",
  Curculionidae = "#66A61E")

col_src_crayfish <- ggthemes::tableau_color_pal()(4)
names(col_src_crayfish) <- c("Ephemeroptera", "Chironomidae", "Gammaridae", "Simuliidae")


recode_treatment <- function(data) {
  if ("Treatment" %in% colnames(data)) {
    data %>%
      dplyr::mutate(
        ALAN = ifelse(
          stringr::str_detect(Treatment, "ALAN"),
          "ALAN",
          "No ALAN"
        ),
        Crayfish = ifelse(
          stringr::str_detect(Treatment, "Crayfish"),
          "Crayfish",
          "No Crayfish"
        )
      ) %>%
      dplyr::relocate(ALAN, Crayfish, .after = Treatment)
  } else {
    data
  }
}

theme_rsm <- function() {
  ggthemes::theme_few() +
    theme(
      text = element_text(size = 11, family = "Helvetica"),
      legend.position = "bottom",
      legend.title = element_text(size = 11, face = "bold"),
      legend.margin = margin(),
      legend.box.margin = margin(),
      strip.text = element_text(face = "bold", margin = margin(1,1,1,1, unit = "mm")),
      panel.grid.major = element_line(linewidth = 0.2, color = "#eeeeee"),
      panel.grid.minor = element_line(linewidth = 0.1, color = "#eeeeee")
    )
}

run_spider_model <- 
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
      write.csv("data/temp/mixture_spider_filtered.csv")
    
    mix_spider <- load_mix_data("data/temp/mixture_spider_filtered.csv",
                                  iso_names = c("d13C", "d15N"),
                                  factors = c("Flume"),
                                  fac_random = c(TRUE),
                                  fac_nested = c(FALSE),
                                  cont_effects = NULL)
    
    src_temp <- 
      src %>%
      dplyr::filter(Treatment == treatment) 
    
    write.csv(src_temp, "data/temp/source_spider_means.csv", row.names = src_temp$Family)
    
    src_spider <- load_source_data("data/temp/source_spider_means.csv",
                                     source_factors = NULL,
                                     conc_dep = FALSE,
                                     data_type = "raw",
                                     mix_spider)
    
    TEF_spider %>%
      dplyr::filter(Family %in% src_spider$source_names) %>%
      dplyr::select(-Family) %>%
      write.csv("data/temp/TEF_spider_filtered.csv")
    
    discr_spider <- load_discr_data("data/temp/TEF_spider_filtered.csv",
                                      mix_spider)
    
    write_JAGS_model(
      filename = "MixSIAR_model_spider.txt",
      resid_err = TRUE,
      process_err = TRUE,
      mix = mix_spider, 
      source = src_spider
    )
    
    spider_model <-
      run_model(
        run = run,
        mix = mix_spider, source = src_spider, discr = discr_spider,
        model_filename = "MixSIAR_model_spider.txt",
        process_err = TRUE,
        resid_err = TRUE
      )
    
    spider_output <- 
      spider_model$BUGSoutput$sims.matrix %>%
      as.data.frame() %>%
      dplyr::select(dplyr::starts_with("p.global"))
    
    colnames(spider_output) <- src_spider$source_names
    
    return(
      list(
        model = spider_model,
        draws = spider_output %>% dplyr::mutate(Treatment = treatment)
      )
    )
  }


run_crayfish_model <- 
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
                                  factors = c("Flume"),
                                  fac_random = c(TRUE),
                                  fac_nested = c(FALSE),
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

