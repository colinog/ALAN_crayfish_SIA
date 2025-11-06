###########################################################################
# Spider and Crayfish Niche Width Analysis Using SIBER
###########################################################################
# Spider Niche Calculation ------------------------------------------------
source_sample = bind_rows(source_spider) |> 
  mutate(across(Treatment, ~factor(.x, levels = c("Control", "ALAN", "Crayfish", "ALAN + Crayfish"))))

mix_sample = bind_rows(mixture_spider) |> 
  mutate(across(Treatment, ~factor(.x, levels = c("Control", "ALAN", "Crayfish", "ALAN + Crayfish"))))



siber_source <- subset(source_sample, 
                       Treatment == "Control" | 
                         Treatment == "ALAN" |
                         Treatment == "Crayfish"|
                         Treatment == "ALAN + Crayfish",
                       select = c("d15N", "d13C", "Treatment", "Family"))

siber_mix <- subset(mix_sample, 
                    Treatment == "Control" | 
                      Treatment == "ALAN" |
                      Treatment == "Crayfish"|
                      Treatment == "ALAN + Crayfish",
                    select = c("d15N", "d13C", "Treatment", "Family"))
######### Rename the columns ############
siber_s = siber_source |>
  dplyr::mutate(community1 = 1) |> 
  dplyr::rename("iso1" = "d13C","iso2" = "d15N", "group" = "Treatment",
                "community" = "community1") |>
  dplyr::select(iso1, iso2, group, community)

siber_m = siber_mix |>
  dplyr::rename("iso1" = "d13C","iso2" = "d15N", "group" = "Treatment",
                "community" = "Family") |>
  dplyr::select(iso1, iso2, group, community)


# Convert to data frame ---------------------------------------------------

siber_s = as.data.frame(siber_s)
siber_m = as.data.frame(siber_m)


# Convert to SIBER object -------------------------------------------------

siber_ob_s <- createSiberObject(siber_s)
siber_ob_m <- createSiberObject(siber_m)

#########SEA#####
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3




ellipses_s <- siberMVN(siber_ob_s, parms, priors) #run the model
ellipses_m <- siberMVN(siber_ob_m, parms, priors)

#group_names <- sub("Tetragnathidae\\.", "", colnames(ellipses.posterior))

SEA_B_s <- siberEllipses(ellipses_s)
SEA_B_m <- siberEllipses(ellipses_m)


SEA_B_s <- SEA_B_s[, c(4, 3, 2, 1)]
SEA_B_m <- SEA_B_m[, c(4, 3, 2, 1)] #spider order


SEA_B_norm = SEA_B_m/SEA_B_s

###################
my_clrs <- matrix(c("#08306b", "#6baed6", "#deebf7",  # Column 1: Light to Dark
                    "#08306b", "#6baed6", "#deebf7",  # Column 2: Light to Dark
                    "#08306b", "#6baed6", "#deebf7",  # Column 3: Light to Dark
                    "#08306b", "#6baed6", "#deebf7"), # Column 4: Light to Dark
                  nrow = 3, ncol = 4)



#"#1f78b4", "#ff7f00", "#33a02c"

# To save the figure ------------------------------------------------------

####Estimate the probability of mean differences between treatment and conttrol######
SEA_B_norm[,1]

diff_A <-  SEA_B_norm[,2] - SEA_B_norm[,1]
prob_A <- mean(diff_A > 0)
BF_A <- prob_A/(1-prob_A)

diff_C <-  SEA_B_norm[,3] - SEA_B_norm[,1]
prob_C <- mean(diff_C > 0)
BF_C <- prob_C/(1-prob_C)

diff_AC <-  SEA_B_norm[,4] - SEA_B_norm[,1]
prob_AC <- mean(diff_AC > 0)
BF_AC <- prob_AC/(1-prob_AC)



sp_pt =ggplotify::as.ggplot(~siberDensityPlot(SEA_B_norm *100, 
                                   xlab = expression(bold("Treatment")),
                                   ylab = expression(bold("Standard ellipse area (%)")),
                                   #ylab = expression(bold("Standard ellipse area " ('\u2030' ^2) )),
                                   ylab.line = 2.7,
                                   #font.axis = 2,
                                   bty = "L",
                                   las = 2,
                                   ylims = c(0,25),
                                   ct = "mean",
                                   clr = my_clrs,
                                   scl = 1,
                                   #xticklabels = c("Crayfish","ALAN + Crayfish"), #Crayfish
                                   xticklabels = c("Control","ALAN","Crayfish","ALAN + Crayfish"), 
                                   #main = "SIBER ellipses of each Treatment"
                                   prn = TRUE,
                                   probs = probs
)+ 
  points(1:ncol(SEA_B_norm*100), group_norm[3,]*100, col="red", pch = "x", lwd = 2))

# Save the figure (Figure 3) ---------------------------------------------------------
ggsave("col_output/spider_SEA.png",
       plot = sp_pt,
       bg = "white",
       width = 6,
       height = 5,
       units = "in",
       dpi = 600)

############################################################################

# Crayfish Niche Calculation ------------------------------------------------

cray_niche = bind_rows(mixture_crayfish) |> 
  mutate(across(Treatment, ~factor(.x, levels = c("Crayfish", "ALAN + Crayfish"))))

###########################################################################
s_ncray = bind_rows(source_crayfish) |> 
  mutate(across(Treatment, ~factor(.x, levels = c("Crayfish", "ALAN + Crayfish"))))

m_ncray = bind_rows(mixture_crayfish) |> 
  mutate(across(Treatment, ~factor(.x, levels = c("Crayfish", "ALAN + Crayfish"))))



siber_scray <- subset(s_ncray,
                      Treatment == "Crayfish"|
                        Treatment == "ALAN + Crayfish",
                      select = c("d15N", "d13C", "Treatment", "Family"))

siber_mcray <- subset(m_ncray,
                      Treatment == "Crayfish"|
                        Treatment == "ALAN + Crayfish",
                      select = c("d15N", "d13C", "Treatment", "Family"))
######### Rename the columns ############
siber_sc = siber_scray |>
  dplyr::mutate(community1 = 1) |> 
  dplyr::rename("iso1" = "d13C","iso2" = "d15N", "group" = "Treatment",
                "community" = "community1") |>
  dplyr::select(iso1, iso2, group, community)

siber_mc = siber_mcray |>
  dplyr::rename("iso1" = "d13C","iso2" = "d15N", "group" = "Treatment",
                "community" = "Family") |>
  dplyr::select(iso1, iso2, group, community)


# Convert to data frame ---------------------------------------------------

siber_sc = as.data.frame(siber_sc)
siber_mc = as.data.frame(siber_mc)


# Convert to SIBER object -------------------------------------------------

siber_ob_sc <- createSiberObject(siber_sc)
siber_ob_mc <- createSiberObject(siber_mc)

#########SEA#####
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3




ellipses_sc <- siberMVN(siber_ob_sc, parms, priors) #run the model
ellipses_mc <- siberMVN(siber_ob_mc, parms, priors)

#group_names <- sub("Tetragnathidae\\.", "", colnames(ellipses.posterior))

SEA_B_sc <- siberEllipses(ellipses_sc)
SEA_B_mc <- siberEllipses(ellipses_mc)


SEA_B_sc <- SEA_B_sc[, c(1, 2)]
SEA_B_mc <- SEA_B_mc[, c(2, 1)] #crayfish column re-order

##Normalized niche Bayesian Stable Isotope Ellipse Area (SEA_b)

SEA_B_normc = SEA_B_mc/SEA_B_sc

###################

mc_clrs <- matrix(c("#8B5E3C", "#C4A484", "#F5F5DC",  # Column 1: Light to Dark
                    "#8B5E3C", "#C4A484", "#F5F5DC"),  # Column 2: Light to Dark
                  nrow = 3, ncol = 2)

#"#1f78b4", "#ff7f00", "#33a02c"

# To save the figure ------------------------------------------------------

####Estimate the probability of mean differences between treatment and conttrol######
SEA_B_normc[,1]

diff_A <-  SEA_B_norm[,2] - SEA_B_norm[,1]
prob_A <- mean(diff_A > 0)

diff_C <-  SEA_B_norm[,3] - SEA_B_norm[,1]
prob_C <- mean(diff_C > 0)

diff_AC <-  SEA_B_norm[,4] - SEA_B_norm[,1]
prob_AC <- mean(diff_AC > 0)

cr_p = as.ggplot(~siberDensityPlot(SEA_B_normc*100, 
                                   xlab = expression(bold("Treatment")),
                                   ylab = expression(bold("Standard Ellipse Area (%)")),
                                   #ylab = expression(bold("Standard ellipse area " ('\u2030' ^2) )),
                                   ylab.line = 2.6,
                                   #font.axis = 2,
                                   bty = "L",
                                   las = 2,
                                   ylims = c(0,10),
                                   ct = "mean",
                                   clr = mc_clrs,
                                   scl = 1,
                                   #xticklabels = c("Crayfish","ALAN + Crayfish"), #Crayfish
                                   xticklabels = c("No ALAN","ALAN"), 
                                   #main = "SIBER ellipses of each Treatment"
                                   prn = TRUE,
                                   probs = probs
)+
  points(1:ncol(SEA_B_normc*100), group_normc[3,]*100, col="red", pch = "x", lwd = 2))
