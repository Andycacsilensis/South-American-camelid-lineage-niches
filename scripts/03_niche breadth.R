# ------------------------------------------------------------------------------
# Andrea G. Castillo
# Niche breadth analysis using convex hulls in global PCA space
#
# Chapter 3:
# "Where is the guanaco irreplaceable? Lineage-level niche structure and
# functional redundancy in wild South American camelids"
#
# Niche breadth is estimated as the area of the convex hull occupied by each
# group in the first two axes of a global environmental PCA. Rarefaction is
# performed using 100 occurrence records per group and 999 iterations.
#
# Date: 14 September 2026

#Load required packages --------------------------------------------------------
library(dplyr)
library(ggplot2)
library(ggtext)
library(grid)
library(geometry) #convex hull analysis

#Define environmental variables -----------------------------------------------
vars.env <- c(
  "bio6",
  "bio14",
  "cmi_max")

#Read occurrence environmental matrices --------------------------------------

# Lama guanicoe
cac.occs <- read.csv(
  "data/environmental_matrices/cacsilensis_occurrences_lonlat_bio_cmimax.csv")
lhyb.occs <- read.csv(
  "data/environmental_matrices/lhybrid_occurrences_lonlat_bio_cmimax.csv")
gua.occs <- read.csv(
  "data/environmental_matrices/guanicoe_occurrences_lonlat_bio_cmimax.csv")

# Vicugna vicugna
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
  "data/environmental_matrices/guanicoe_10K_Mpoints_lonlat_bio_cmimax.csv")

# Vicugna vicugna
men.bg <- read.csv(
  "data/environmental_matrices/mensalis_5K_Mpoints_lonlat_bio_cmimax.csv")
vhyb.bg <- read.csv(
  "data/environmental_matrices/vhybrid_5K_Mpoints_lonlat_bio_cmimax.csv")
vic.bg <- read.csv(
  "data/environmental_matrices/vicugna_5K_Mpoints_lonlat_bio_cmimax.csv")

#Remove background records with missing environmental values -----------------
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


#Construct global environmental space -----------------------------------------
#Combine the accessible environmental space (M) of all six groups.
env.global <- rbind(
  cac.bg[, vars.env],
  lhyb.bg[, vars.env],
  gua.bg[, vars.env],
  men.bg[, vars.env],
  vhyb.bg[, vars.env],
  vic.bg[, vars.env])

#Perform global PCA ------------------------------------------------------------
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

# Retain PC1 and PC2 ------------------------------------------------------------
sp.cac <- cac.occ.pca[, 1:2]
sp.lhyb <- lhyb.occ.pca[, 1:2]
sp.gua <- gua.occ.pca[, 1:2]

sp.men <- men.occ.pca[, 1:2]
sp.vhyb <- vhyb.occ.pca[, 1:2]
sp.vic <- vic.occ.pca[, 1:2]

# Check PCA projections ---------------------------------------------------------
# Confirm that no missing values are present in the projected coordinates.
missing.pca <- sapply(
  list(
    cacsilensis = sp.cac,
    guanaco_contact = sp.lhyb,
    guanicoe = sp.gua,
    mensalis = sp.men,
    vicuna_contact = sp.vhyb,
    vicugna = sp.vic
  ),
  function(x) sum(!complete.cases(x)))

missing.pca
#All groups should return zero missing values.

#Calculate raw niche breadth ---------------------------------------------------
#Function to calculate the area of a two-dimensional convex hull.

hull_area <- function(xy) {
  xy <- as.matrix(xy[, 1:2])
  hull <- chull(xy)
  hull <- c(hull, hull[1])
  x <- xy[hull, 1]
  y <- xy[hull, 2]
  0.5 * abs(
    sum(
      x[-1] * y[-length(y)] -
        x[-length(x)] * y[-1]))}

#Store occurrence PCA coordinates in a single list.
pca.occ.list <- list(
  cacsilensis = sp.cac,
  guanaco_contact = sp.lhyb,
  guanicoe = sp.gua,
  mensalis = sp.men,
  vicuna_contact = sp.vhyb,
  vicugna = sp.vic)

#Calculate raw convex-hull area for each group.
raw.breadth <- sapply(
  pca.occ.list,
  hull_area)

raw.breadth

#Summarize raw niche breadth ---------------------------------------------------
raw.breadth.df <- data.frame(
  Group = names(pca.occ.list),
  n = sapply(pca.occ.list, nrow),
  Raw_breadth = raw.breadth)
raw.breadth.df

#Export raw niche breadth results.
write.csv(raw.breadth.df, "output/niche_breadth/niche_breadth_raw.csv",
  row.names = FALSE)


#Rarefaction of niche breadth --------------------------------------------------

#Note: Each group is randomly subsampled to 100 occurrence records.
#Convex-hull area is calculated for each subsample over 999 iterations.
rarefy_hull <- function(
    pca.coords,
    n.sample = 100,
    n.rep = 999
) {
  if (nrow(pca.coords) < n.sample) {
    stop(
      "n.sample is larger than the number of available occurrences.")}                                
  replicate(
    n.rep,
    {   ids <- sample(
        seq_len(nrow(pca.coords)),
        size = n.sample,
        replace = FALSE)
      samp <- pca.coords[
        ids,
        1:2,
        drop = FALSE]
      hull_area(samp)})}


# Set seed to ensure reproducibility of random subsampling.
set.seed(123)

# Run rarefaction for all six groups.
rare.results <- lapply(
  pca.occ.list,
  rarefy_hull,
  n.sample = 100,
  n.rep = 999)


#Validate rarefaction results --------------------------------------------------
#NOTE: A convex hull calculated from a subset of points cannot exceed the convex
#hull calculated from the complete occurrence dataset. Therefore, all rarefied 
#estimates should be less than or equal to the corresponding raw niche breadth.

validation <- sapply(
  names(rare.results),
  function(g) {
    all(
      rare.results[[g]] <=
        raw.breadth[g] + 1e-12)})

validation

# Compare the largest rarefied hull with the corresponding full convex hull.
# Difference should be <= 0 for all groups.
validation.detail <- data.frame(
  Group = names(rare.results),
  Raw = as.numeric(
    raw.breadth[names(rare.results)]),
  Max_rarefied = sapply(
    rare.results,
    max))
validation.detail$Difference <-
  validation.detail$Max_rarefied -
  validation.detail$Raw

validation.detail

# Summarize rarefied niche breadth ---------------------------------------------
breadth.rarefied <- data.frame(
  Group = names(rare.results),
  n_original = sapply(
    pca.occ.list,
    nrow),
  n_rarefied = 100,
  Raw_breadth = as.numeric(
    raw.breadth[names(rare.results)]),
  Rarefied_mean = sapply(
    rare.results,
    mean),
  Lower_95 = sapply(
    rare.results,
    quantile,
    probs = 0.025),
  Upper_95 = sapply(
    rare.results,
    quantile,
    probs = 0.975))

#Calculate the proportion of occurrence records retained after rarefaction
#and the proportion of raw convex-hull area retained by the rarefied estimate.

breadth.rarefied <- breadth.rarefied %>%
  mutate(
    Records_retained = n_rarefied / n_original,
    Area_retained = Rarefied_mean / Raw_breadth)
breadth.rarefied

#Export rarefied niche breadth results.
write.csv(breadth.rarefied,
  "output/niche_breadth/niche_breadth_rarefied_n100_999iter.csv", row.names = FALSE)

#Plot niche breadth results ----------------------------------------------------
breadth.plot <- breadth.rarefied %>%
  mutate(
    Group = recode(
      Group,
      "cacsilensis" =
        "<i>L. g. cacsilensis</i>",
      "guanaco_contact" =
        "Guanaco contact-zone",
      "guanicoe" =
        "<i>L. g. guanicoe</i>",
      "mensalis" =
        "<i>V. v. mensalis</i>",
      "vicuna_contact" =
        "Vicuña contact-zone",
      "vicugna" =
        "<i>V. v. vicugna</i>"),
    
    Species = case_when(
      Group %in% c(
        "<i>L. g. cacsilensis</i>",
        "Guanaco contact-zone",
        "<i>L. g. guanicoe</i>"
      ) ~ "Guanaco",
      TRUE ~ "Vicuña"))


# Define group order.
breadth.plot$Group <- factor(
  breadth.plot$Group,
  levels = c(
    "<i>L. g. cacsilensis</i>",
    "Guanaco contact-zone",
    "<i>L. g. guanicoe</i>",
    "<i>V. v. mensalis</i>",
    "Vicuña contact-zone",
    "<i>V. v. vicugna</i>"))

breadth.plot$Species <- factor(
  breadth.plot$Species,
  levels = c(
    "Guanaco",
    "Vicuña"))

# Define group colors.
group_colors <- c(
  "<i>L. g. cacsilensis</i>" = "#4a0c6b",
  "Guanaco contact-zone" = "#781c6d",
  "<i>L. g. guanicoe</i>" = "#a52c60",
  "<i>V. v. mensalis</i>" = "#cf4446",
  "Vicuña contact-zone" = "#ed6925",
  "<i>V. v. vicugna</i>" = "#fb9b06")

# Define horizontal positions for sample-size labels and plot limits.
n_position <- max(breadth.plot$Upper_95) * 1.18
x_max <- max(breadth.plot$Upper_95) * 1.35

# Create niche breadth figure.
p_breadth <- ggplot(
  breadth.plot,
  aes(
    x = Rarefied_mean,
    y = Group,
    colour = Group)) +
  #95% percentile interval
  geom_errorbarh(
    aes(
      xmin = Lower_95,
      xmax = Upper_95),
    height = 0,
    linewidth = 0.8) +

  #Mean rarefied niche breadth
  geom_point(
    size = 4) +
  
  # Original occurrence sample size
  geom_text(
    aes(
      x = n_position,
      label = paste0(
        "n = ",
        n_original)),
    hjust = 0,
    size = 3.6) +
  
  # One panel per species
  facet_grid(
    Species ~ .,
    scales = "free_y",
    space = "free_y") +
  
  scale_colour_manual(
    values = group_colors) +
  
  scale_x_continuous(
    limits = c(
      0,
      x_max),
    expand = c(
      0,
      0)) +
  
  labs(
    x = expression(
      "Rarefied niche breadth (PCA-space units"^2 * ")"),
    y = NULL) +
  
  theme_bw(
    base_size = 11) +
  
  theme(legend.position = "none",
    
    #Vertical grid lines only
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(
      linewidth = 0.3,
      colour = "grey85"),
    
    #Italicize taxonomic names using markdown
    axis.text.y = ggtext::element_markdown(
      size = 10,
      colour = "black"),
    
    axis.text.x = element_text(
      colour = "black"),
    
    axis.title.x = element_text(
      size = 11,
      margin = margin(
        t = 8)),
    
    # Species labels
    strip.background = element_rect(
      fill = "white",
      colour = "black"),
    
    strip.text.y = element_text(
      face = "bold",
      size = 11),
    
    strip.placement = "outside",
    
    panel.spacing.y = unit(
      0.15,
      "cm"))


# Display figure.

p_breadth