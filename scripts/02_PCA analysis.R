# ------------------------------------------------------------------------------
# Andrea G. Castillo
# Global PCA and visualization of environmental space
#
# Chapter 3:
# "Where is the guanaco irreplaceable? Lineage-level niche structure and
# functional redundancy in wild South American camelids"
#
# This script constructs a common environmental PCA using the accessible
# environmental space (M) of all six groups and projects occurrence records
# onto the resulting PCA space.
#
# Date: 14 September 2026

#Load required packages 
library(ggplot2)
library(dplyr)

#Define environmental variables -----------------------------------------------
vars.env <- c(
  "bio6",
  "bio14",
  "cmi_max")

#Read occurrence environmental matrices --------------------------------------
#Lama guanicoe
cac.occs <- read.csv(
  "data/environmental_matrices/cacsilensis_occurrences_lonlat_bio_cmimax.csv")
lhyb.occs <- read.csv(
  "data/environmental_matrices/lhybrid_occurrences_lonlat_bio_cmimax.csv")
gua.occs <- read.csv(
  "data/environmental_matrices/guanicoe_occurrences_lonlat_bio_cmimax.csv")

#Vicugna vicugna
men.occs <- read.csv(
  "data/environmental_matrices/mensalis_occurrences_lonlat_bio_cmimax.csv")
vhyb.occs <- read.csv(
  "data/environmental_matrices/vhybrid_occurrences_lonlat_bio_cmimax.csv")
vic.occs <- read.csv(
  "data/environmental_matrices/vicugna_occurrences_lonlat_bio_cmimax.csv")


#Remove occurrence records with missing environmental values -----------------
cac.occs <- cac.occs[
  complete.cases(cac.occs[, vars.env]),]
lhyb.occs <- lhyb.occs[
  complete.cases(lhyb.occs[, vars.env]),]
gua.occs <- gua.occs[
  complete.cases(gua.occs[, vars.env]),]
men.occs <- men.occs[
  complete.cases(men.occs[, vars.env]),]
vhyb.occs <- vhyb.occs[
  complete.cases(vhyb.occs[, vars.env]),]
vic.occs <- vic.occs[
  complete.cases(vic.occs[, vars.env]),]


#Read background environmental matrices --------------------------------------
# Lama guanicoe
cac.bg <- read.csv(
  "data/environmental_matrices/cacsilensis_10K_Mpoints_lonlat_bio_cmimax.csv")
lhyb.bg <- read.csv(
  "data/environmental_matrices/lhybrid_10K_Mpoints_lonlat_bio_cmimax.csv")
gua.bg <- read.csv(
  "data/environmental_matrices/guanicoe_10K_Mpoints_lonlat_bio_cmimax6.csv")

# Vicugna vicugna
men.bg <- read.csv(
  "data/environmental_matrices/mensalis_5K_Mpoints_lonlat_bio_cmimax.csv")
vhyb.bg <- read.csv(
  "data/environmental_matrices/vhybrid_5K_Mpoints_lonlat_bio_cmimax.csv")
vic.bg <- read.csv(
  "data/environmental_matrices/vicugna_5K_Mpoints_lonlat_bio_cmimax.csv")


#Remove background records with missing environmental values ------------------
cac.bg <- cac.bg[
  complete.cases(cac.bg[, vars.env]),]
lhyb.bg <- lhyb.bg[
  complete.cases(lhyb.bg[, vars.env]),]
gua.bg <- gua.bg[
  complete.cases(gua.bg[, vars.env]),]
men.bg <- men.bg[
  complete.cases(men.bg[, vars.env]),]
vhyb.bg <- vhyb.bg[
  complete.cases(vhyb.bg[, vars.env]),]
vic.bg <- vic.bg[
  complete.cases(vic.bg[, vars.env]),]


# Construct a global environmental PCA -----------------------------------------
# Combine the accessible environmental space (M) of all six groups.
# Only environmental variables are included in the PCA.

env.global <- rbind(
  cac.bg[, vars.env],
  lhyb.bg[, vars.env],
  gua.bg[, vars.env],
  men.bg[, vars.env],
  vhyb.bg[, vars.env],
  vic.bg[, vars.env])

#Perform PCA using centered and standardized environmental variables.
pca.global <- prcomp(
  env.global,
  center = TRUE,
  scale. = TRUE)


#Examine explained variance and variable loadings.
summary(pca.global)
pca.global$rotation

# Project occurrence records onto the global PCA space -------------------------

cac.occ.pca <- predict(
  pca.global,
  cac.occs[, vars.env])

lhyb.occ.pca <- predict(
  pca.global,
  lhyb.occs[, vars.env])

gua.occ.pca <- predict(
  pca.global,
  gua.occs[, vars.env])

men.occ.pca <- predict(
  pca.global,
  men.occs[, vars.env])

vhyb.occ.pca <- predict(
  pca.global,
  vhyb.occs[, vars.env])

vic.occ.pca <- predict(
  pca.global,
  vic.occs[, vars.env])

# Combine PC1 and PC2 scores for all groups ------------------------------------
pca.global.df <- bind_rows(
  
  data.frame(
    PC1 = cac.occ.pca[, 1],
    PC2 = cac.occ.pca[, 2],
    Group = "L. g. cacsilensis"),
  
  data.frame(
    PC1 = lhyb.occ.pca[, 1],
    PC2 = lhyb.occ.pca[, 2],
    Group = "Guanaco contact-zone"),
  
  data.frame(
    PC1 = gua.occ.pca[, 1],
    PC2 = gua.occ.pca[, 2],
    Group = "L. g. guanicoe"),
  
  data.frame(
    PC1 = men.occ.pca[, 1],
    PC2 = men.occ.pca[, 2],
    Group = "V. v. mensalis"),
  
  data.frame(
    PC1 = vhyb.occ.pca[, 1],
    PC2 = vhyb.occ.pca[, 2],
    Group = "Vicuña contact-zone"),
  
  data.frame(
    PC1 = vic.occ.pca[, 1],
    PC2 = vic.occ.pca[, 2],
    Group = "V. v. vicugna"))

# Define group order ------------------------------------------------------------
pca.global.df$Group <- factor(
  pca.global.df$Group,
  levels = c(
    "L. g. cacsilensis",
    "Guanaco contact-zone",
    "L. g. guanicoe",
    "V. v. mensalis",
    "Vicuña contact-zone",
    "V. v. vicugna"))


# Define color palette ----------------------------------------------------------
lineage_colors <- c(
  "L. g. cacsilensis" = "#4a0c6b",
  "Guanaco contact-zone" = "#781c6d",
  "L. g. guanicoe" = "#a52c60",
  "V. v. mensalis" = "#cf4446",
  "Vicuña contact-zone" = "#ed6925",
  "V. v. vicugna" = "#fb9b06")

# Calculate group centroids -----------------------------------------------------
centroids <- pca.global.df %>%
  group_by(Group) %>%
  summarise(
    PC1 = mean(PC1),
    PC2 = mean(PC2),
    .groups = "drop")


# Extract percentage of variance explained -------------------------------------
pca_summary <- summary(pca.global)
PC1_var <- pca_summary$importance[2, 1] * 100
PC2_var <- pca_summary$importance[2, 2] * 100

# Plot global environmental PCA -------------------------------------------------
p_pca <- ggplot(
  pca.global.df,
  aes(
    x = PC1,
    y = PC2,
    colour = Group,
    fill = Group
  )
) +
  
  # 95% confidence ellipses
  stat_ellipse(
    geom = "polygon",
    level = 0.95,
    alpha = 0.18,
    linewidth = 0.9
  ) +
  
  # Group centroids
  geom_point(
    data = centroids,
    aes(
      x = PC1,
      y = PC2,
      colour = Group,
      fill = Group
    ),
    size = 4,
    shape = 21,
    stroke = 1,
    inherit.aes = FALSE
  ) +
  
  scale_colour_manual(
    values = lineage_colors
  ) +
  
  scale_fill_manual(
    values = lineage_colors
  ) +
  
  labs(
    x = paste0(
      "PC1 (",
      round(PC1_var, 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(PC2_var, 1),
      "%)"
    ),
    colour = "Group",
    fill = "Group"
  ) +
  
  theme_bw(
    base_size = 11
  ) +
  
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(
      colour = "black"
    ),
    legend.position = "right",
    legend.title = element_text(
      face = "bold"
    )
  )


# Display figure ----------------------------------------------------------------

p_pca