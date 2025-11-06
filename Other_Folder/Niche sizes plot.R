packageVersion("SIBER")
citation("SIBER")

source_sample = bind_rows(source_spider) |> 
  mutate(across(Treatment, ~factor(.x, levels = c("Control", "ALAN", "Crayfish", "ALAN + Crayfish"))))

mix_sample = bind_rows(mixture_spider) |> 
  mutate(across(Treatment, ~factor(.x, levels = c("Control", "ALAN", "Crayfish", "ALAN + Crayfish"))))

#str(all_sample)


siber.biplots<-ggplot(all_sample, aes(x = d13C, y = d15N, colour = Family)) +
  geom_point(alpha = 0.7, size=2) +
  facet_grid(. ~ Treatment) +
  theme_bw() +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)")))
siber.biplots

# Precompute means and standard errors by group
summary_data <- all_sample %>%
  group_by(Treatment, Family) %>%
  summarise(
    mean_d13C = mean(d13C, na.rm = TRUE),
    mean_d15N = mean(d15N, na.rm = TRUE),
    se_d13C = sd(d13C, na.rm = TRUE) / sqrt(n()),
    se_d15N = sd(d15N, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

ggplot(summary_data, aes(x = mean_d13C, y = mean_d15N, colour = Family)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_d15N - se_d15N, ymax = mean_d15N + se_d15N), width = 0.1) +
  geom_errorbarh(aes(xmin = mean_d13C - se_d13C, xmax = mean_d13C + se_d13C), height = 0.1) +
  facet_grid(. ~ Treatment) +
  theme_bw() +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)")))



# siber.group<-subset(all_sample, 
#                     Treatment == "ALAN + Crayfish" | 
#                       Treatment == "Crayfish")[,c("d15N", "d13C", "Treatment", "Family")]

siber.group <- subset(all_sample, 
                        Treatment == "Control" | 
                        Treatment == "ALAN" |
                        Treatment == "Crayfish"|
                        Treatment == "ALAN + Crayfish",
                      select = c("d15N", "d13C", "Treatment", "Family"))


sum(is.na(siber.group$d15N))  # Count missing values in d15N
sum(is.na(siber.group$d13C))  # Count missing values in d13C
head(siber.group)

siber.group.biplots<-ggplot(siber.group, aes(x = d13C, y = d15N, colour = Family, group = Treatment)) +
  geom_point(alpha = 0.7, size=2) +
  facet_grid(. ~ Treatment) +
  theme_bw() +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)")))
siber.group.biplots



# Assuming your dataset has columns 'd15N', 'd13C', and 'Treatment'
ggplot(data = siber.group, aes(x = d13C, y = d15N, color = Family)) +
  geom_point() +  # Plot points
  stat_ellipse(aes(group = Family), position = "identity", level = 0.95) +  # Plot ellipses
  facet_wrap(~ Treatment) +
  theme_minimal()  # A clean theme



#siber.group$Treatment <- factor(siber.group$Treatment)

siber.group.biplots +
  geom_point(aes(color = Family)) +  # Ensure points are colored by Treatment
  stat_ellipse(aes(color = Treatment), position = "identity", level = 0.4, linewidth = 2) +  # 40% ellipse
  stat_ellipse(aes(color = Treatment), position = "identity", level = 0.95, linewidth = 1)   # 95% ellipse


siber.group.biplots+
  stat_ellipse(aes(group = Treatment),position="identity", level=0.4, linewidth=2)+
  stat_ellipse(aes(group = Treatment),position="identity", level=0.95, linewidth=1)

siber.group1 = siber.group |>
  dplyr::rename("iso1" = "d13C","iso2" = "d15N", "group" = "Treatment",
         "community" = "Family") |>
  dplyr::select(iso1, iso2, group, community)

siber.group2 = as.data.frame(siber.group1)
head(siber.group2)
colnames(siber.group2)
table(siber.group2$community)

siber.g.example <- createSiberObject(siber.group2)

community.hulls.args <- list(col = 1, lty = 0, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)
group.hull.args      <- list(lty = 2, col = "grey20")

###################
# Convert the SIBER object to a data frame
#df <- as.data.frame(siber.g.example$original.data)


par(mfrow=c(1,1))
#color_blind_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
siber.g.example$group.names
legend <- c("Control", "ALAN", "Crayfish", "ALAN + Crayfish")
gp <- palette(c("#E69F00", "#56B4E9", "#009E73", "#D55E00"))
pchs <- c(1,2,3,4)

#gp = c("#E69F00", "#56B4E9", "#009E73", "#D55E00")
plotSiberObject(siber.g.example,
                ax.pad = 2, 
                hulls = F, community.hulls.args, 
                ellipses = T, group.ellipses.args,
                group.hulls = T, group.hull.args,
                bty = "L",
                col = gp,
                pch = pchs,
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030'))

legend("topright", legend, col = colours,pch = pchs,  bty = "n")


siber.group2 %>% ggplot(aes(iso1, iso2, fill = group, colour = group))+
  geom_point(size=3)+
  stat_ellipse()+
  xlab(expression({delta}^13*C~'\u2030'))+
  ylab(expression({delta}^15*N~'\u2030'))

##############Estimating overlap###########
sp_gp = siber.group2
sp_gp$community_id <- as.numeric(as.factor(sp_gp$community))
sp_gp$group_id <- as.numeric(interaction(sp_gp$community, sp_gp$group, drop = TRUE))

sp_gpc = sp_gp |> 
  select(-c(group, community)) |> 
  rename(group = group_id, community = community_id) |> 
  select(iso1, iso2, group, community)

sp_gpc = as.data.frame(sp_gpc)
m_spider = createSiberObject(sp_gpc)

ellipse1 = as.character(c(1.1))
ellipse2 = as.character(c(1.2))
ellipse3 = as.character(c(1.3))
ellipse4 = as.character(c(1.4))

class(ellipse1)
str(m_spider)
sp_overlap <- maxLikOverlap(ellipse1, ellipse2, m_spider)
print(overlap)

# # Suppose you want to compare group 1 and 2 in community 1
# ####Estimate the percentage overlap#####
# area.1 <- 1.580332
# area.2 <- 2.003313
# n_overlap <- 1.579251
# 
# #############
# p_overlap_1in2 <- (n_overlap / area.1) * 100
# print(p_overlap_1in2) # ~99.93%
# pe_overlap_2in1 <- (n_overlap / area.2) * 100  
# print(p_overlap_2in1)# ~78.82%

##########new from Andrew#####
# the overlap betweeen the corresponding 95% prediction ellipses is given by:
#compare 1 = Control ellipse to the treatment ellipses: where 2 = ALAN, 3 = Crayfish, 4 = ALAN + Crayfish
ellipse95.overlap1_2 <- maxLikOverlap(ellipse1, ellipse2, m_spider, 
                                   p.interval = 0.95, n = 1000)

ellipse95.overlap1_3 <- maxLikOverlap(ellipse1, ellipse3, m_spider, 
                                      p.interval = 0.95, n = 1000)

ellipse95.overlap1_4 <- maxLikOverlap(ellipse1, ellipse4, m_spider, 
                                      p.interval = 0.95, n = 1000)

# so in this case, the overlap as a proportion of the non-overlapping area of 
# the two ellipses, would be
prop.95.over12 <- ellipse95.overlap1_2[3] / (ellipse95.overlap1_2[2] + 
                                          ellipse95.overlap1_2[1] -
                                          ellipse95.overlap1_2[3])
print(prop.95.over12) # 41.09%

prop.95.over13 <- ellipse95.overlap1_3[3] / (ellipse95.overlap1_3[2] + 
                                            ellipse95.overlap1_3[1] -
                                            ellipse95.overlap1_3[3])
print(prop.95.over13) # 42.17%

prop.95.over14 <- ellipse95.overlap1_4[3] / (ellipse95.overlap1_4[2] + 
                                            ellipse95.overlap1_4[1] -
                                            ellipse95.overlap1_4[3])
print(prop.95.over14) # 38.47%


######################
# Extract original data from SIBER object
df <- as.data.frame(siber.g.example$original.data)

# # Rename columns (adjust if needed)
# colnames(df) <- c("d13C", "d15N", "Community", "Group")
# 
# # Define colors and shapes manually
#group_colors <- c("#E7298A", "#7570B3", "#1B9E77", "#D95F02")  # Modify as needed
group_colors <- c("#E69F00", "#0072B2", "#009E73", "#D55E00")
pchs <- c(16, 17, 18, 15)  # Different point shapes per group

# Function to compute convex hulls for each group
compute_hull <- function(df) {
  df %>%
    group_by(group) %>%
    slice(chull(iso1, iso2))  # Compute convex hull points
}

# Get hulls for each group
hulls <- compute_hull(df)

library(nicheROVER)

aggregate(df[1:2], df[3], mean)

nsamples <- 1e3
system.time({
  df.par <- tapply(1:nrow(df), df$group,
                     function(ii) niw.post(nsamples = nsamples, X = df[ii,1:2]))
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






# Plot using ggplot2
ggplot(df, aes(x = iso1, y = iso2, color = factor(group), shape = factor(group))) +
  geom_point(size = 3, alpha = 0.8) +  # Add points
  #stat_ellipse(level = 0.95, linetype = "dashed")+
  geom_polygon(data = hulls, aes(group = group, fill = factor(group)), alpha = 0.1) +  # Group hulls
  scale_color_manual(values = group_colors) +  # Custom colors
  scale_fill_manual(values = group_colors) +  # Match fill colors to points
  scale_shape_manual(values = pchs) +  # Custom shapes
  labs(
    x = expression(delta^13*C~("\u2030")), 
    y = expression(delta^15*N~("\u2030")),
    color = "Treatment",
    shape = "Treatment",
    fill = "Treatment"
  ) +
  theme_classic() +  # Clean theme with readable font size
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "bottom",  # Move legend for clarity
    legend.title = element_text(face = "bold"),  # Bold legend titles
    #panel.border = element_rect(fill = NA, color = "black", linewidth = 1)  # Add a border
  )

ggsave(
  "col_output/tet_hull.png",
  height = 12,
  width = 16,
  units = "cm"
)


plotGroupEllipses(siber.g.example, n = 100, p.interval = 0.95, ci.mean = T,
                  lty = 1, lwd = 2)

community.ML <- communityMetricsML(siber.g.example) ##For community matrix of which we have one
print(community.ML)

str(siber.g.example)
group.ML <- groupMetricsML(siber.g.example)##Group matrix by treatment (4)
print(group.ML)


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

ellipses.posterior <- siberMVN(siber.g.example, parms, priors) #run the model

#group_names <- sub("Tetragnathidae\\.", "", colnames(ellipses.posterior))

SEA.B <- siberEllipses(ellipses.posterior)

colnames(SEA.B)

SEA.B <- SEA.B[, c(4, 3, 2, 1)] #spider order

#SEA.B <- SEA.B[, c(2, 1)]

head(SEA.B)

# c_value <- group.ML[3, c("Tetragnathidae.Control", "Tetragnathidae.ALAN", 
#                                     "Tetragnathidae.Crayfish", "Tetragnathidae.ALAN + Crayfish")]
# ##Changing the density plot colour

# colors = viridis::viridis(9)
# clrs = matrix(rep(colors, times = ncol(SEA.B)), nrow = 9, ncol = ncol(SEA.B))

library(magick)

plt1 = image_read("col_output/tet_hull.png")
plt2 = image_read("col_output/spider_SEA.png")

# Resize second to match height of first
img2_res <- image_resize(plt2, geometry = paste0("x", plt1_info$height))

# Combine horizontally
#combined <- image_append(c(img1, img2_resized))

plt_comb = image_append(c(plt1, img2_res))


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
plt1_info = image_info(plt1)
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
image_write(plt_ann1, path = "col_output/anotated_plot.png", 
            format = "png", quality = 100)


my_clrs <- matrix(c("#E69F00", "#D55E00", "#009E73",  # Column 1: Light to Dark
                    "#E69F00", "#D55E00", "#009E73",  # Column 2: Light to Dark
                    "#E69F00", "#D55E00", "#009E73",  # Column 3: Light to Dark
                    "#E69F00", "#D55E00", "#009E73"), # Column 4: Light to Dark
                  nrow = 3, ncol = 4)







siberDensityPlot(SEA.B, 
                 xlab = "Treatment",
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 ylims = c(0,2),
                 ct = "median",
                 clr = my_clrs,
                 scl = 1,
                 #xticklabels = c("Crayfish","ALAN + Crayfish"), #Crayfish
                 xticklabels = c("Control","ALAN","Crayfish","ALAN + Crayfish"), 
                 #main = "SIBER ellipses of each Treatment"
                 prn = TRUE,
                 probs = probs
)

group.ML1 <- group.ML[, c(4, 3, 2, 1)] #Change the arrangement of the group matrix

# Add red x's for the ML estimated SEA-c
points(1:ncol(SEA.B), group.ML1[3,], col="red", pch = "x", lwd = 2)

?siberDensityPlot

# Then reopen a device to save it
png("col_output/spider_SEA.png", width = 16, height = 12, units = "cm", res = 300)

dev.off()


#png("col_output/spider_SEA.png", width = 16.5, height = 12, units = "cm", res = 300)
#dev.off()



#legend("topright", legend = "SEAc", col = "red", pch = "x", lwd = 0, bty = "0", cex = 0.5)

#####End of Spider #####
####Estimating the treatment overlap in percentage ######

overlaps <- maxLikOverlap(
  SEA.B,   # group means
  siber_object$ML.cov   # group covariances
)














#points(1:ncol(SEA.B), c_value, col="red", pch = "x", lwd = 2)
############################ SIBER#################
# Perform MANOVA
manova_results <- manova(cbind(d15N, d13C) ~ Treatment, data = all_sample)

# Print MANOVA summary
summary(manova_results)
print(manova_results)

# Run ANOVA for each isotope separately
aov_d15N <- aov(d15N ~ Treatment , data = all_sample)
aov_d13C <- aov(d13C ~ Treatment, data = all_sample)

# Tukey HSD test for pairwise differences
tukey_d15N <- TukeyHSD(aov_d15N)
tukey_d13C <- TukeyHSD(aov_d13C)

# Print results
print(tukey_d15N)
print(tukey_d13C)


# Create a SIBER object from your dataset
s_all = all_sample |> 
  filter(Treatment %in% c("Crayfish", "ALAN", "Control", "ALAN + Crayfish"))
siber.data <- data.frame(
  iso1 = all_sample$d13C,  # First isotope
  iso2 = all_sample$d15N,  # Second isotope
  group = as.numeric(as.factor(all_sample$Family)),  # Convert family to numeric groups
  community = as.numeric(as.factor(all_sample$Treatment))  # Assign a single community if working in one ecosystem
)

# Convert to SIBER object
siber.obj <- createSiberObject(siber.data)

ellipses.posterior <- siberMVN(siber.obj, parms, priors) #run the model

SEA.B <- siberEllipses(ellipses.posterior)

#Normalize SEA by mean of SEAb of all treatment for easy comparability between treatment.
m_SEA <- mean(SEA.B)  
nor_SEA.B <- SEA.B / m_SEA  

###################################

str(SEA.B)

group_names <- c("ALAN", "ALAN + Crayfish", "Control", "Crayfish")
colnames(SEA.B) <- group_names


siberDensityPlot(SEA.B, 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 ylims = c(0,2.5)
                 #main = "SIBER ellipses on each group"
)

# Define custom group names
group_names <- c("ALAN", "ALAN + Crayfish", "Control", "Crayfish")  # Replace with actual names

# Add custom x-axis labels
axis(1, at = 1:4, labels = group_names, las = 1)

# Calculate Bayesian Standard Ellipse Area (SEAb)
SEA <- siberEllipses(ellipses.posterior)

# Plot SEAb ellipses
plotSiberObject(siber.obj, ax.pad = 2, hulls = FALSE, ellipses = TRUE, group.hulls = TRUE)



####################
# of the three plotting functions.
community.hulls.args <- list(col = "black", lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)
group.hull.args      <- list(lty = 2, col = "grey20")


any(is.na(siber.obj$ML.mu))   # Should return FALSE
any(is.infinite(siber.obj$ML.mu))  # Should return FALSE
det_cov <- det(siber.obj$ML.mu)
print(siber.obj$ML.mu)
str(siber.obj)
siber.obj$ML.mu
# plot the raw data

par(mfrow=c(1,1))
plotSiberObject(siber.obj,
                ax.pad = 2, 
                hulls = TRUE, community.hulls.args, 
                ellipses = FALSE, group.ellipses.args,
                group.hulls = FALSE, group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030')
)
plotGroupEllipses(siber.obj, n = 1000, p.interval = 0.05,
                  ci.mean = F, lty = 1, lwd = 2)

############
# extract the posterior means

mu.post <- extractPosteriorMeans(siber.g.example, ellipses.posterior)

# calculate the corresponding distribution of layman metrics
layman.B <- bayesianLayman(mu.post)


# --------------------------------------
# Visualise the first community
# --------------------------------------

# drop the 3rd column of the posterior which is TA using -3.
siberDensityPlot(layman.B[[1]], 
                 xticklabels = colnames(layman.B[[1]]), 
                 bty="L", ylim = c(0,2))

# add the ML estimates (if you want). Extract the correct means 
# from the appropriate array held within the overall array of means.
comm1.layman.ml <- laymanMetrics(siber.obj$ML.mu[[1]][1,1,],
                                 siber.obj$ML.mu[[1]][1,2,]
)

# again drop the 3rd entry which relates to TA
points(1:5, comm1.layman.ml$metrics[-3], 
       col = "red", pch = "x", lwd = 2)


# --------------------------------------
# Visualise the second community
# --------------------------------------
siberDensityPlot(layman.B[[1]][ , -3], 
                 xticklabels = colnames(layman.B[[1]][ , -3]), 
                 bty="L", ylim = c(0,1))

# add the ML estimates. (if you want) Extract the correct means 
# from the appropriate array held within the overall array of means.
comm2.layman.ml <- laymanMetrics(siber.obj$ML.mu[[2]][1,1,],
                                 siber.obj$ML.mu[[2]][1,2,]
)
points(1:5, comm2.layman.ml$metrics[-3], 
       col = "red", pch = "x", lwd = 2)


# --------------------------------------
# Alternatively, pull out TA from both and aggregate them into a 
# single matrix using cbind() and plot them together on one graph.
# --------------------------------------

# go back to a 1x1 panel plot
par(mfrow=c(1,1))

# Now we only plot the TA data. We could address this as either
# layman.B[[1]][, "TA"]
# or
# layman.B[[1]][, 3]
siberDensityPlot(cbind(layman.B[[1]][ , "TA"], 
                       layman.B[[1]][ , "TA"]),
                 xticklabels = c("Community 1", "Community 2"), 
                 bty="L", ylim = c(0, 10),
                 las = 1,
                 ylab = "TA - Convex Hull Area",
                 xlab = "")


# We can take the code from just above for TA and simply 
# swap out the metric of choice for whatever we want. 
# Here I do this for NND and take care to change the 
# y axis labels too.
siberDensityPlot(cbind(#layman.B[[3]][ , "NND"], 
                       layman.B[[1]][ , "NND"],
                       #layman.B[[2]][ , "NND"],
                       layman.B[[1]][ , "NND"]),
                 #xticklabels = c("Control", "ALAN", "ALAN + Crayfish", "Crayfish"), 
                 bty="L", ylim = c(1, 2.5),
                 las = 1,
                 ylab = "NND (nearest neighbour distance)",
                 xlab = "")

siberDensityPlot(cbind(layman.B[[1]][ , "dY_range"], 
                       layman.B[[1]][ , "dY_range"],
                       layman.B[[1]][ , "dY_range"],
                       layman.B[[1]][ , "dY_range"]),
                 xticklabels = c("Control", "ALAN", "ALAN + Crayfish", "Crayfish"), 
                 bty="L", ylim = c(0, 5),
                 las = 1,
                 ylab = "dX",
                 xlab = "")
dNr1.lt.dNr2 <- sum(layman.B[[3]][,"dY_range"] < 
                      layman.B[[4]][,"dY_range"]) / 
  length(layman.B[[3]][,"dY_range"])

print(dNr1.lt.dNr2)

siberDensityPlot(cbind(layman.B[[1]][,"TA"], layman.B[[2]][,"TA"]),
                 xticklabels = c("Community 1", "Community 2"), 
                 bty="L", ylim = c(0,20),
                 las = 1,
                 ylab = "TA - Convex Hull Area",
                 xlab = "")
