# ------------------------------------------------------------------------------
# Andrea G. Castillo
# Niche overlap analysis using "ecospat" package (Di Cola et al., 2017) 
#
# Chapter 3:
# "Where is the guanaco irreplaceable? Lineage-level niche structure and
# functional redundancy in wild South American camelids"
#
# Pairwise niche overlap is quantified in a common environmental PCA space
# using Schoener's D and Warren's I (URL:https://doi.org/10.1111/ecog.02671)
#
# Date: 14 September 2026

# Load required packages --------------------------------------------------------
library(ecospat)
library(dplyr)
library(ggplot2)
library(grid)

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

#Organize occurrence and background data --------------------------------------
occ.list <- list(
  cacsilensis = cac.occs,
  guanaco_contact = lhyb.occs,
  guanicoe = gua.occs,
  mensalis = men.occs,
  vicuna_contact = vhyb.occs,
  vicugna = vic.occs)

bg.list <- list(
  cacsilensis = cac.bg,
  guanaco_contact = lhyb.bg,
  guanicoe = gua.bg,
  mensalis = men.bg,
  vicuna_contact = vhyb.bg,
  vicugna = vic.bg)

#Remove records with missing environmental values -----------------------------
occ.list <- lapply(
  occ.list,
  function(x) {
    x[
      complete.cases(x[, vars.env]),]})

bg.list <- lapply(
  bg.list,
  function(x) {
    x[
      complete.cases(x[, vars.env]),]})

# Check sample sizes after removing missing values ------------------------------
sapply(occ.list, nrow)
sapply(bg.list, nrow)


#Construct global environmental PCA -------------------------------------------
#Combine the accessible environmental space (M) of all six groups.
env.global <- do.call(
  rbind,
  lapply(
    bg.list,
    function(x) x[, vars.env]))

#Perform PCA using centered and standardized environmental variables.
pca.global <- prcomp(
  env.global,
  center = TRUE,
  scale. = TRUE)

#Examine explained variance and variable loadings.
summary(pca.global)
pca.global$rotation

#Project backgrounds and occurrences onto the global PCA space ----------------
bg.pca <- lapply(
  bg.list,
  function(x) {
    predict(
      pca.global,
      newdata = x[, vars.env])})

occ.pca <- lapply(
  occ.list,
  function(x) {
    predict(
      pca.global,
      newdata = x[, vars.env])})

#Check dimensions of projected datasets.
lapply(bg.pca, dim)
lapply(occ.pca, dim)

#Build ecospat environmental grids --------------------------------------------
#Note: For each pairwise comparison, the environmental space is defined by the
#combined background of both groups. Group-specific environmental availability
#is represented by the corresponding background within this common space.
build_z <- function(
    group1,
    group2,
    bg.pca,
    occ.pca,
    R = 100) {
  
  glob <- rbind(
    bg.pca[[group1]][, 1:2],
    bg.pca[[group2]][, 1:2])
  
  z1 <- ecospat.grid.clim.dyn(
    glob = glob,
    glob1 = bg.pca[[group1]][, 1:2],
    sp = occ.pca[[group1]][, 1:2],
    R = R )
  
  z2 <- ecospat.grid.clim.dyn(
    glob = glob,
    glob1 = bg.pca[[group2]][, 1:2],
    sp = occ.pca[[group2]][, 1:2],
    R = R)
  
  list(
    z1 = z1,
    z2 = z2)}

#Calculate pairwise niche overlap ---------------------------------------------
#Niche overlap is quantified using Schoener's D and Warren's I.

calc_overlap <- function(
    group1,
    group2,
    bg.pca,
    occ.pca,
    R = 100
) {
  
  z <- build_z(
    group1 = group1,
    group2 = group2,
    bg.pca = bg.pca,
    occ.pca = occ.pca,
    R = R)
  
  ov <- ecospat.niche.overlap(
    z$z1,
    z$z2,
    cor = TRUE)
  
  data.frame(
    group1 = group1,
    group2 = group2,
    D = ov$D,
    I = ov$I)}

#Run all pairwise comparisons --------------------------------------------------
groups <- c(
  "cacsilensis",
  "guanaco_contact",
  "guanicoe",
  "mensalis",
  "vicuna_contact",
  "vicugna")

pairs <- combn(
  groups,
  2,
  simplify = FALSE)

overlap.results <- do.call(
  rbind,
  lapply(
    pairs,
    function(pair) {
      
      calc_overlap(
        group1 = pair[1],
        group2 = pair[2],
        bg.pca = bg.pca,
        occ.pca = occ.pca,
        R = 100)}))

overlap.results

#Check overlap results ---------------------------------------------------------
range(overlap.results$D)
range(overlap.results$I)

overlap.results %>%
  arrange(D)

# Export pairwise niche overlap results -----------------------------------------

write.csv(overlap.results,"output/niche_overlap/niche_overlap_D_I.csv",
  row.names = FALSE)


# Prepare Schoener's D matrix for plotting -------------------------------------
# Define group order.

group_order <- c(
  "cacsilensis",
  "guanaco_contact",
  "guanicoe",
  "mensalis",
  "vicuna_contact",
  "vicugna")

# Define group labels.
group_labels_expr <- expression(
  italic(L.~g.~cacsilensis),
  "Guanaco contact-zone",
  italic(L.~g.~guanicoe),
  italic(V.~v.~mensalis),
  "Vicuña contact-zone",
  italic(V.~v.~vicugna))

# Make pairwise Schoener's D results symmetric.
overlap.symmetric <- bind_rows(
  overlap.results %>%
    transmute(
      row = group1,
      col = group2,
      D = D),
  
  overlap.results %>%
    transmute(
      row = group2,
      col = group1,
      D = D))

# Build complete Schoener's D matrix.
overlap.matrix <- expand.grid(
  row = group_order,
  col = group_order,
  stringsAsFactors = FALSE
) %>%
  left_join(
    overlap.symmetric,
    by = c(
      "row",
      "col")
  ) %>%
  mutate(
    row_id = match(
      row,
      group_order),
    col_id = match(
      col,
      group_order ))


#Retain lower triangle only.
overlap.lower <- overlap.matrix %>%
  filter(
    row_id > col_id
  ) %>%
  mutate(
    x = col_id,
    y = 7 - row_id )

# Create diagonal values --------------------------------------------------------
# Each group has complete niche overlap with itself (D = 1).
diagonal.df <- data.frame(
  x = 1:6,
  y = 6:1,
  D = 1)

# Plot Schoener's D overlap matrix ----------------------------------------------
p_overlap <- ggplot(
  overlap.lower,
  aes(
    x = x,
    y = y,
    fill = D)) +
  
  # Lower triangle: pairwise Schoener's D
  geom_tile(
    width = 1,
    height = 1,
    colour = "white",
    linewidth = 0.8) +
  
  geom_text(
    aes(
      label = sprintf(
        "%.2f",
        D),
      colour = D > 0.30),
    size = 4,
    show.legend = FALSE) +
  
  # Diagonal: D = 1
  # Plotted separately so that diagonal values do not affect the
  # colour scale used for pairwise comparisons.
  geom_tile(
    data = diagonal.df,
    aes(
      x = x,
      y = y),
    inherit.aes = FALSE,
    width = 1,
    height = 1,
    fill = "black",
    colour = "white",
    linewidth = 0.8) +
  
  geom_text(
    data = diagonal.df,
    aes(
      x = x,
      y = y,
      label = "1.00"),
    inherit.aes = FALSE,
    colour = "white",
    size = 4) +
  
  # Grayscale for pairwise Schoener's D
  scale_fill_gradient(
    low = "white",
    high = "black",
    limits = c(
      0,
      0.50),
    breaks = seq(
      0,
      0.50,
      0.10),
    name = "Schoener's D") +
  
  # Switch numeric labels to white on darker cells.
  scale_colour_manual(
    values = c(
      "FALSE" = "black",
      "TRUE" = "white")) +
  
  # Axes
  scale_x_continuous(
    breaks = 1:6,
    labels = group_labels_expr,
    limits = c(
      0.5,
      8.9),
    expand = c(
      0,
      0)) +
  
  scale_y_continuous(
    breaks = 1:6,
    labels = rev(
      group_labels_expr),
    limits = c(
      0.5,
      6.5),
    expand = c(
      0,
      0)) +
  
  # Guanaco bracket
  annotate(
    "segment",
    x = 6.65,
    xend = 6.65,
    y = 3.55,
    yend = 6.45,
    linewidth = 0.7) +
  
  annotate(
    "segment",
    x = 6.45,
    xend = 6.65,
    y = 6.45,
    yend = 6.45,
    linewidth = 0.7) +
  
  annotate(
    "segment",
    x = 6.45,
    xend = 6.65,
    y = 3.55,
    yend = 3.55,
    linewidth = 0.7) +
  
  annotate(
    "text",
    x = 6.85,
    y = 5,
    label = "Guanaco",
    hjust = 0,
    size = 4) +
  
  annotate(
    "text",
    x = 6.85,
    y = 4.65,
    label = "Lama guanicoe",
    hjust = 0,
    fontface = "italic",
    size = 3.7) +
  
  # Vicuña bracket
  annotate(
    "segment",
    x = 6.65,
    xend = 6.65,
    y = 0.55,
    yend = 3.45,
    linewidth = 0.7) +
  
  annotate(
    "segment",
    x = 6.45,
    xend = 6.65,
    y = 3.45,
    yend = 3.45,
    linewidth = 0.7) +
  
  annotate(
    "segment",
    x = 6.45,
    xend = 6.65,
    y = 0.55,
    yend = 0.55,
    linewidth = 0.7) +
  
  annotate(
    "text",
    x = 6.85,
    y = 2,
    label = "Vicuña",
    hjust = 0,
    size = 4) +
  
  annotate(
    "text",
    x = 6.85,
    y = 1.65,
    label = "Vicugna vicugna",
    hjust = 0,
    fontface = "italic",
    size = 3.7) +
  
  # Appearance
  coord_fixed(
    ratio = 1,
    clip = "off") +
  
  labs(
    x = NULL,
    y = NULL) +
  
  guides(
    fill = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = unit(
        5,
        "cm"),
      barwidth = unit(
        0.6,
        "cm"))) +
  
  theme_minimal(
    base_size = 11) +
  
  theme(
    panel.grid = element_blank(),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      colour = "black",
      size = 10),
    
    axis.text.y = element_text(
      colour = "black",
      size = 10),
    
    axis.ticks = element_blank(),
    
    legend.position = "right",
    
    legend.title = element_text(
      size = 11),
    
    legend.text = element_text(
      size = 9),
    
    plot.margin = margin(
      t = 10,
      r = 10,
      b = 10,
      l = 10))


# Display figure ----------------------------------------------------------------

p_overlap