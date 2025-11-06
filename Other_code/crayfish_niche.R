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




library(ggplotify)

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

mtext("(c)", side = 3, line = 1.2, at = par("usr")[1] + 0.01 * diff(par("usr")[1:2]), cex = 1.2, font = 2)


library(ggplotify)
library(patchwork)

p_overlay <- as.ggplot(expression({
  points(1:ncol(SEA_B_normc*100), group_normc[3,]*100, col="red", pch="x", lwd=2)
}))

p_main + p_overlay


cr_p +
  geom_point(
    aes(x = 1:ncol(SEA_B_normc*100), y = group_normc[3,]*100),
    color = "red",
    shape = 4,    # shape 4 is "x"
    size = 3,     # adjust for lwd equivalent
    stroke = 1.5
  )


# Now manually add the bold x-axis labels
# axis(1, 
#      at = 1:4, 
#      labels = c("Control", "ALAN", "Crayfish", "ALAN + Crayfish"), 
#      font = 2, 
#      cex.axis = 0.9)  # Adjust size if needed


group_sc <- groupMetricsML(siber_ob_sc)##Group matrix by treatment (4)
group_mc <- groupMetricsML(siber_ob_mc)
group_sc_ML <- group_sc[, c(1,2)]
group_mc_ML <- group_mc[, c(2,1)]

group_normc = group_mc_ML/group_sc_ML



#group.ML1 <- group.ML[, c(4, 3, 2, 1)] #Change the arrangement of the group matrix

# Add red x's for the ML estimated SEA-c
points(1:ncol(SEA_B_normc*100), group_normc[3,]*100, col="red", pch = "x", lwd = 2)

# Then reopen a device to save it
png("col_output/new_cray_SEA.png", width = 16, height = 10, units = "cm", res = 600)

dev.off()


##################################Crayfish trophic niche plot #################













cray_niche |> 
  ggplot( aes(x = d13C, y = d15N, colour = Family)) +
  geom_point(alpha = 0.7, size=2) +
  facet_grid(. ~ Treatment) +
  theme_bw() +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)")))

# Assuming your dataset has columns 'd15N', 'd13C', and 'Treatment'
cray_niche |> 
  ggplot( aes(x = d13C, y = d15N, color = Family)) +
  geom_point() +  # Plot points
  stat_ellipse(aes(group = Family), position = "identity", level = 0.95) +  # Plot ellipses
  facet_wrap(~ Treatment) +
  theme_minimal()  # A clean theme

#############Resouce SEA########
res_cray = source_crayfish |> 
  dplyr::filter(Treatment %in% c("Crayfish", "ALAN + Crayfish")) |>
  dplyr::rename("iso1" = "d13C","iso2" = "d15N", "group" = "Treatment",
                "community" = "Family") |>
  dplyr::select(iso1, iso2, group, community)

table(res_cray1$group)
res_cray1 = as.data.frame(res_cray)

res_gp = createSiberObject(res_cray1)

plotSiberObject(res_gp,
                ax.pad = 2, 
                hulls = F, community.hulls.args, 
                ellipses = T, group.ellipses.args,
                group.hulls = F, group.hull.args,
                bty = "L",
                #col = gp_cray,
                #pch = c_pchs,
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030'))


###Rename the isogroup for the crayfish niche by treatment
cray_gp = cray_niche |>
  dplyr::rename("iso1" = "d13C","iso2" = "d15N", "group" = "Treatment",
                "community" = "Family") |>
  dplyr::select(iso1, iso2, group, community)

cray_gp1 = as.data.frame(cray_gp)
head(cray_gp1)
colnames(cray_gp1)
table(cray_gp1$community)

cray_gp_eg <- createSiberObject(cray_gp1)

community.hulls.args <- list(col = 1, lty = 0, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)
group.hull.args      <- list(lty = 2, col = "grey20")

########Estimating the overlap percentage ######
###################
# Convert the SIBER object to a data frame
#df <- as.data.frame(siber.g.example$original.data)


par(mfrow=c(1,1))
#color_blind_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cray_gp_eg$group.names
cray_leg <- c("Crayfish", "ALAN + Crayfish")
gp_cray <- palette(c("#E69F00", "#56B4E9"))
c_pchs <- c(3,4)

#gp = c("#E69F00", "#56B4E9", "#009E73", "#D55E00")
plotSiberObject(cray_gp_eg,
                ax.pad = 2, 
                hulls = F, community.hulls.args, 
                ellipses = F, group.ellipses.args,
                group.hulls = T, group.hull.args,
                bty = "L",
                col = gp_cray,
                pch = c_pchs,
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030'))

legend("topright", cray_leg, col = gp_cray, pch = c_pchs,  bty = "n")

#I decided to use ggplot for the convex hull as the original
#plot function was unable to distinguish the group by color

######################
# Extract original data from SIBER object
df_c <- as.data.frame(cray_gp_eg$original.data)

# # Rename columns (adjust if needed)
# colnames(df) <- c("d13C", "d15N", "Community", "Group")
# 
# # Define colors and shapes manually
#group_colors <- c("#E7298A", "#7570B3", "#1B9E77", "#D95F02")  # Modify as needed
cray_clr <- c("#E69F00", "#0072B2")
c_pchs <- c(16, 17)  # Different point shapes per group

# Function to compute convex hulls for each group
cray_hull <- function(df_c) {
  df_c %>%
    group_by(group) %>%
    slice(chull(iso1, iso2))  # Compute convex hull points
}

# Get hulls for each group
crayfish_hull <- cray_hull(df_c)

dev.off()
# Plot the group hulls showing the crayfish niche under the two treatment
ggplot(df_c, aes(x = iso1, y = iso2, color = factor(group), shape = factor(group))) +
  geom_point(size = 3, alpha = 0.8) +  # Add points
  #stat_ellipse(level = 0.95, linetype = "dashed")+
  geom_polygon(data = crayfish_hull, aes(group = group, fill = factor(group)), alpha = 0.1) +  # Group hulls
  scale_color_manual(values = cray_clr,
                     labels = c("Crayfish" = "no-ALAN", "ALAN + Crayfish" = "ALAN")) +  # Custom colors
  scale_fill_manual(values = cray_clr,
                    labels = c("Crayfish" = "no-ALAN", "ALAN + Crayfish" = "ALAN")) +  # Match fill colors to points
  scale_shape_manual(values = c_pchs,
                     labels = c("Crayfish" = "no-ALAN", "ALAN + Crayfish" = "ALAN")) +  # Custom shapes
  labs(
    x = expression(delta^13*C~("\u2030")), 
    y = expression(delta^15*N~("\u2030")),
    color = "Treatment",
    shape = "Treatment",
    fill = "Treatment"
  ) +
  theme_classic(base_size = 14) +  # Clean theme with readable font size
  theme(
    legend.position = "bottom",  # Move legend for clarity
    legend.title = element_text(face = "bold")  # Bold legend titles
    #panel.border = element_rect(fill = NA, color = "black", linewidth = 1)  # Add a border
  )

ggsave(
  "col_output/cray_hull.png",
  height = 12,
  width = 16.5,
  units = "cm"
)



# plotGroupEllipses(siber.g.example, n = 100, p.interval = 0.95, ci.mean = T,
#                   lty = 1, lwd = 2)

# community.ML <- communityMetricsML(cray_g_eg) ##For community matrix of which we have one
# print(community.ML)

group.ML_cray <- groupMetricsML(cray_gp_eg)##Group matrix by treatment (4)
print(group.ML_cray)


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

post_ellipse <- siberMVN(cray_gp_eg, parms, priors) #run the model

#group_names <- sub("Tetragnathidae\\.", "", colnames(ellipses.posterior))

SEA.cray <- siberEllipses(post_ellipse)

SEA.cray <- SEA.cray[, c(2, 1)] ##Change the column order

SEA_df = data.frame(SEA.cray)

names(SEA_df)[1] <- "Crayfish"  # changes the 2nd column name

names(SEA_df)[2] <- "ALAN + Crayfish"

#########Estimate of poesterior mean and confidence intervals#######

#group1_SEA <- SEA_df[, 2]  # adjust for your group

# Posterior mean
posterior_mean <- mean(SEA_df$Crayfish)

# 95% credible interval
ci <- quantile(SEA_df$Crayfish, probs = c(0.025, 0.975))


###############




##Changing the density plot colour


clr = t(col2rgb(c("red", "green", "blue", "orange")))
probs = c(95, 75, 50)
#length(probs) == nrow(clr)

base_plot <- recordPlot({})

dev.off()
#png("col_output/b_panel.png", width = 6, height = 5, units = "in", res = 300)

siberDensityPlot(SEA_df, 
                 xlab = "Treatment",
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 ylims = c(0,0.5),
                 ct = "median",
                 clr = my_clrs,
                 scl = 1,
                 #clr = clrs,
                 #xticklabels = c("Crayfish","ALAN + Crayfish"), #Crayfish
                 xticklabels = c("Crayfish","ALAN + Crayfish"), 
                 #main = "SIBER ellipses of each Treatment"
                 prn = TRUE,
                 probs = probs,
                 #clr = clr
                 
)

group.MLc <- group.ML_cray[, c(2, 1)] #Change the arrangement of the group matrix

# Add red x's for the ML estimated SEA-c
points(1:ncol(SEA.cray), group.MLc[3,], col="red", pch = "x", lwd = 2) 
dev.off()
###The Standard ellipses area included the small sample size corrected (SEAc)
# Get group names and their numeric index
cray_gp_t = cray_gp
cray_gp_t$community_id <- as.numeric(as.factor(cray_gp_t$community))
cray_gp_t$group_id <- as.numeric(interaction(cray_gp_t$community, cray_gp_t$group, drop = TRUE))

cray_gpc = cray_gp_t |> 
  select(-c(group, community)) |> 
  rename(group = group_id, community = community_id)

cray_gpc = as.data.frame(cray_gpc)
m_cray = createSiberObject(cray_gpc)

ellipse1 = as.character(c(1.1))
ellipse2 = as.character(c(1.2))

class(ellipse1)
str(m_cray)
overlap <- maxLikOverlap(ellipse1, ellipse2,m_cray)
print(overlap)

# # Suppose you want to compare group 1 and 2 in community 1
# ####Estimate the percentage overlap#####
# area.1 <- 1.580332
# area.2 <- 2.003313
# n_overlap <- 1.579251
# 
# p_overlap_1in2 <- (n_overlap / area.1) * 100
# print(p_overlap_1in2) # ~99.93%
# pe_overlap_2in1 <- (n_overlap / area.2) * 100  
# print(p_overlap_2in1)# ~78.82%

###########new from Andrew#####
# the overlap betweeen the corresponding 95% prediction ellipses is given by:
ellipse95.overlap <- maxLikOverlap(ellipse1, ellipse2, m_cray, 
                                   p.interval = 0.95, n = 1000)

# so in this case, the overlap as a proportion of the non-overlapping area of 
# the two ellipses, would be
prop.95.over <- ellipse95.overlap[3] / (ellipse95.overlap[2] + 
                                          ellipse95.overlap[1] -
                                          ellipse95.overlap[3])


#legend("topright", legend = "SEAc", col = "red", pch = "x", lwd = 0, bty = "0", cex = 0.5)

cr1 = image_read("col_output/cray_hull.png")
cr2 = image_read("col_output/b_panel.png")

# Resize second to match height of first
plt1_info = image_info(cr1)
img2_res <- image_resize(cr2, geometry = paste0("x", plt1_info$height))

# Combine horizontally
#combined <- image_append(c(img1, img2_resized))

plt_comb = image_append(c(cr1, img2_res))


plt_ann = image_annotate(
  plt_comb,
  "(a)",
  size = 50,
  font = "Arial-Bold", 
  gravity = "northwest", 
  location = "+20+20", 
  color = "black"
)

# Add (b) on the right — shift x position to right half

half_width <- plt1_info$width
plt_ann1 <- image_annotate(
  plt_ann, 
  "(b)", 
  size = 50,
  font = "Arial-Bold",
  gravity = "northwest", 
  location = paste0("+", half_width + 20, "+20"), 
  color = "black"
)

#img12_res = image_resize(plt_ann1, "x10000")

# Save to file
image_write(plt_ann1, path = "col_output/anotated_cray_plot.png", 
            format = "png", quality = 100)





#####End of crayfish #####
library(patchwork)
cowplot::plot_grid(p2, base_plot,
                   labels = 'AUTO',
                   hjust = 0, vjust = 1)


png("col_output/b_panel.png", width = 16, height = 12, units = "cm", res = 300)




# Define custom axis labels using plotmath
xlab_custom <- bquote(mu[delta] ~ "(" * .(unique(df.par$neutron)) * ") " * .(unique(df.par$element)))
ylab_custom <- bquote(Sigma[delta])  # example, can be customized similarly

# Plot
niche.par.plot(df.par,
               col = clrs,
               plot.index = 1,
               xlab = xlab_custom,
               ylab = ylab_custom)


##############Niche Rover#########
library(nicheROVER)
names(mixture_crayfish)
nr_df = mixture_crayfish |> 
  dplyr::select(c(Treatment,ALAN, d13C, d15N, cp_mm, sex )) |> 
  mutate(tr_sex = str_c(ALAN,sex, sep = "_"))
  
aggregate(nr_df[2:3], df[3], mean)


cr_par <- nr_df %>% 
  split(.$tr_sex) %>% 
  map(~ select(., d13C, d15N)) %>% 
  map(~niw.post(nsamples = nsamples, X = .))

df_mu <- map(cr_par, pluck, 1) %>% 
  imap(~ as_tibble(.x) %>% 
         mutate( 
           metric = "mu", 
           tr_sex = .y
         )
  ) %>%
  bind_rows() %>% 
  mutate(
    tr_sex = factor(tr_sex, 
                     levels = c("No ALAN_male", "No ALAN_female", "ALAN_male", "ALAN_female"))
  ) %>% 
  group_by(tr_sex) %>% 
  mutate(
    sample_number = 1:1000
  ) %>% 
  ungroup() |> 
  rename(Treatment = tr_sex)

df_mu_long <- df_mu %>% 
  pivot_longer(cols = -c(metric, Treatment, sample_number), 
               names_to = "isotope", 
               values_to = "mu_est") %>% 
  mutate(
    element = case_when(
      isotope == "d15N" ~ "N",
      isotope == "d13C" ~ "C",
    ), 
    neutron = case_when(
      isotope == "d15N" ~ 15,
      isotope == "d13C" ~ 13,
    ) 
  )


df_sigma <- map(cr_par, pluck, 2) %>%
  imap(~ as_tibble(.x) %>%
         mutate(
           metric = "sigma",
           id = c("d15N", "d13C"),
           tr_sex = .y
         )
  ) %>%
  bind_rows() %>%
  pivot_longer(cols = -c("id", "tr_sex", "metric"),
               names_to = "isotope",
               values_to = "post_sample"
  )  %>%
  separate(isotope, into = c("isotopes", "sample_number"), sep = "\\.") |> 
  rename(Treatment = tr_sex)

df_sigma_cn <- df_sigma %>%
  filter(id != isotopes)


posterior_plots <- df_mu_long %>%
  split(.$isotope) %>%
  imap(
    ~ ggplot(data = ., aes(x = mu_est)) +
      geom_density(aes(fill = Treatment), alpha = 0.5) +
      scale_fill_viridis_d(begin = 0.25, end = 0.75,
                           option = "D", name = "Treatment by sex") +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.title.x =  element_markdown(),
            axis.title.y =  element_markdown(),
            legend.position = "none"
      ) +
      labs(
        x = paste("\u00b5<sub>\U03B4</sub>", "<sub><sup>",
                  unique(.$neutron), "</sup></sub>",
                  "<sub>",unique(.$element), "</sub>", sep = ""),
        y = paste0("p(\u00b5 <sub>\U03B4</sub>","<sub><sup>",
                   unique(.$neutron), "</sub></sup>",
                   "<sub>",unique(.$element),"</sub>",
                   " | X)"), sep = "")
  )

posterior_plots$d15N +
  posterior_plots$d13C +
  theme(legend.position = "right") 

df_sigma_cn <- df_sigma_cn %>%
  mutate(
    element_id = case_when(
      id == "d15N" ~ "N",
      id == "d13C" ~ "C",
    ),
    neutron_id = case_when(
      id == "d15N" ~ 15,
      id == "d13C" ~ 13,
    ),
    element_iso = case_when(
      isotopes == "d15N" ~ "N",
      isotopes == "d13C" ~ "C",
    ),
    neutron_iso = case_when(
      isotopes == "d15N" ~ 15,
      isotopes == "d13C" ~ 13,
    )
  )

sigma_plots <- df_sigma_cn %>%
  group_split(id, isotopes) %>%
  imap(
    ~ ggplot(data = ., aes(x = post_sample)) +
      geom_density(aes(fill = sex), alpha = 0.5) +
      scale_fill_viridis_d(begin = 0.25, end = 0.75,
                           option = "D", name = "Sex") +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.title.x =  element_markdown(),
            axis.title.y =  element_markdown(),
            legend.position = "none"
      ) +
      labs(
        x = paste("\U03A3","<sub>\U03B4</sub>",
                  "<sub><sup>", unique(.$neutron_id), "</sub></sup>",
                  "<sub>",unique(.$element_id),"</sub>"," ",
                  "<sub>\U03B4</sub>",
                  "<sub><sup>", unique(.$neutron_iso), "</sub></sup>",
                  "<sub>",unique(.$element_iso),"</sub>", sep = ""),
        y = paste("p(", "\U03A3","<sub>\U03B4</sub>",
                  "<sub><sup>", unique(.$neutron_id), "</sub></sup>",
                  "<sub>",unique(.$element_id),"</sub>"," ",
                  "<sub>\U03B4</sub>",
                  "<sub><sup>", unique(.$neutron_iso), "</sub></sup>",
                  "<sub>",unique(.$element_iso),"</sub>", " | X)", sep = ""),
      )
  )

sigma_plots[[1]] + 
  theme(legend.position = c(0.6, 0.85))






nsamples <- 1e3
system.time({
  df.par <- tapply(1:nrow(df), nr_df$cp_mm,
                   function(ii) niw.post(nsamples = nsamples, X = nr_df[ii,2:3]))
})

# various parameter plots
clrs <- c("#E69F00", "#0072B2", "#009E73", "#D55E00") # colors for each species

# mu1 (del15N), mu2 (del13C), and Sigma12
par(mar = c(4, 4, .5, .1)+.1, mfrow = c(1,3))
niche.par.plot(df.par, col = clrs, plot.index = 1)
niche.par.plot(df.par, col = clrs, plot.index = 2)
niche.par.plot(df.par, col = clrs, plot.index = 1:2)
legend("topright", legend = names(df.par), fill = clrs, cex = 0.6, box.lty = 0.6)

#########################
# all mu (del15N, del13C, del34S)
niche.par.plot(df.par, col = clrs, plot.mu = TRUE, plot.Sigma = FALSE)
legend("topleft", legend = names(df.par), fill = clrs, cex = 0.5)




###########################
par(mar = c(4.2, 4.2, 2, 1)+.1)
niche.par.plot(df.par, col = clrs, plot.mu = TRUE, plot.Sigma = TRUE)
legend("topright", legend = names(df.par), fill = clrs, cex = 0.5)


clrs <- c("black", "red", "blue", "orange") # colors for each species
nsamples <- 10
df.par <- tapply(1:nrow(df), df$group,
                 function(ii) niw.post(nsamples = nsamples, X = df[ii,1:2]))

# format data for plotting function
df.data <- tapply(1:nrow(df), df$group, function(ii) X = df[ii,1:2])

niche.plot(niche.par = df.par, niche.data = df.data, pfrac = .05,
           iso.names = expression(delta^{15}*N, delta^{13}*C, delta^{34}*S),
           col = clrs, xlab = expression("Isotope Ratio (per mil)"))


# niche overlap plots for 95% niche region sizes
nsamples <- 1000
df.par <- tapply(1:nrow(df), df$group,
                 function(ii) niw.post(nsamples = nsamples, X = df[ii,1:2]))

# Overlap calculation.  use nsamples = nprob = 10000 (1e4) for higher accuracy.
# the variable over.stat can be supplied directly to the overlap.plot function

over.stat <- overlap(df.par, nreps = nsamples, nprob = 1e3, alpha = c(.95, 0.99))

#The mean overlap metrics calculated across iteratations for both niche 
#region sizes (alpha = .95 and alpha = .99) can be calculated and displayed in an array.
over.mean <- apply(over.stat, c(1:2,4), mean)*100
round(over.mean, 2)

over.cred <- apply(over.stat*100, c(1:2, 4), quantile, prob = c(.025, .975), na.rm = TRUE)
round(over.cred[,,,1]) # display alpha = .95 niche region

# Overlap plot.Before you run this, make sure that you have chosen your 
#alpha level.
clrs <- c("black", "red", "blue", "orange") # colors for each species
over.stat <- overlap(df.par, nreps = nsamples, nprob = 1e3, alpha = .95)
overlap.plot(over.stat, col = clrs, mean.cred.col = "turquoise", equal.axis = TRUE,
             xlab = "Overlap Probability (%) -- Niche Region Size: 95%")

# posterior distribution of (mu, Sigma) for each species
nsamples <- 1000
df.par <- tapply(1:nrow(df_c), df_c$group,
                 function(ii) niw.post(nsamples = nsamples, X = df[ii,1:2]))

# posterior distribution of niche size by species
df.size <- sapply(df.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = .95)
})

# point estimate and standard error
rbind(est = colMeans(df.size),
      se = apply(df.size, 2, sd))

dev.off()
# boxplots
clrs <- c("black", "red", "blue", "orange") # colors for each species
boxplot(df.size, col = clrs, pch = 16, cex = .5,
        ylab = "Niche Size", xlab = "Species")
