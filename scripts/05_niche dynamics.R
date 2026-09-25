# ------------------------------------------------------------------------------
# Andrea G. Castillo
# Niche dynamics analysis using "ecospat" package (Di Cola et al., 2017)
#
# Chapter 3:
# "Where is the guanaco irreplaceable? Lineage-level niche structure and
# functional redundancy in wild South American camelids"
#
# Expansion, stability, and unfilling are quantified between guanaco and
# vicuña groups in a common environmental PCA space. Guanaco groups are used
# as the reference (z1) and vicuña groups as the focal group (z2).
# (URL:https://doi.org/10.1111/ecog.02671)
#
# The main analysis is restricted to analogous environmental conditions.
# A sensitivity analysis subsequently evaluates the effect of including
# non-analogous environmental conditions.
#
# Date: 14 September 2026

#Load required packages --------------------------------------------------------
library(ecospat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)

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

#Check sample sizes after removing missing values.
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


#Define functions for niche dynamics ------------------------------------------
#NOTE: Build ecospat environmental grids for each pairwise comparison.
#z1 = guanaco group (reference)
#z2 = vicuña group (focal)

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
    R = R)
  
  z2 <- ecospat.grid.clim.dyn(
    glob = glob,
    glob1 = bg.pca[[group2]][, 1:2],
    sp = occ.pca[[group2]][, 1:2],
    R = R)
  
  list(
    z1 = z1,
    z2 = z2)}

#Calculate expansion, stability, and unfilling under analogous
#environmental conditions.
calc_dyn_analogue <- function(
    group1,
    group2,
    bg.pca,
    occ.pca,
    R = 100) {
  
  pair <- build_z(
    group1 = group1,
    group2 = group2,
    bg.pca = bg.pca,
    occ.pca = occ.pca,
    R = R)
  
  dyn <- ecospat.niche.dyn.index(
    z1 = pair$z1,
    z2 = pair$z2,
    intersection = 0)
  
  data.frame(
    Expansion = as.numeric(
      dyn$dynamic.index.w["expansion"]),
    Stability = as.numeric(
      dyn$dynamic.index.w["stability"]),
    Unfilling = as.numeric(
      dyn$dynamic.index.w["unfilling"]))}


#Define pairwise comparisons ---------------------------------------------------
#All comparisons are directional:
#group1 = guanaco reference (z1)
#group2 = vicuña focal group (z2).
comparisons <- data.frame(
  group1 = c(
    "cacsilensis",
    "cacsilensis",
    "cacsilensis",
    "guanaco_contact",
    "guanaco_contact",
    "guanaco_contact",
    "guanicoe",
    "guanicoe",
    "guanicoe"),
  
  group2 = c(
    "mensalis",
    "vicuna_contact",
    "vicugna",
    "mensalis",
    "vicuna_contact",
    "vicugna",
    "mensalis",
    "vicuna_contact",
    "vicugna"),
  
  Comparison = c(
    "cacsilensis-mensalis",
    "cacsilensis-vicuna_contact",
    "cacsilensis-vicugna",
    "guanaco_contact-mensalis",
    "guanaco_contact-vicuna_contact",
    "guanaco_contact-vicugna",
    "guanicoe-mensalis",
    "guanicoe-vicuna_contact",
    "guanicoe-vicugna"))

#Run main niche dynamics analysis ---------------------------------------------
dynamic.results <- lapply(
  seq_len(nrow(comparisons)),
  function(i) {
    result <- calc_dyn_analogue(
      group1 = comparisons$group1[i],
      group2 = comparisons$group2[i],
      bg.pca = bg.pca,
      occ.pca = occ.pca,
      R = 100 )
    data.frame(
      Comparison = comparisons$Comparison[i],
      result)}
) %>%
  bind_rows()

dynamic.results

#Export main niche dynamics results -------------------------------------------
write.csv( dynamic.results,"output/niche_dynamics/niche_dynamics_analogue.csv",
  row.names = FALSE)

#Prepare data for plotting -----------------------------------------------------
dyn.plot <- dynamic.results %>%
  mutate(
    Comparison = recode(
      Comparison,
      "cacsilensis-mensalis" =
        "L. g. cacsilensis – V. v. mensalis",
      "cacsilensis-vicuna_contact" =
        "L. g. cacsilensis – Vicuña contact-zone",
      "cacsilensis-vicugna" =
        "L. g. cacsilensis – V. v. vicugna",
      "guanaco_contact-mensalis" =
        "Guanaco contact-zone – V. v. mensalis",
      "guanaco_contact-vicuna_contact" =
        "Guanaco contact-zone – Vicuña contact-zone",
      "guanaco_contact-vicugna" =
        "Guanaco contact-zone – V. v. vicugna",
      "guanicoe-mensalis" =
        "L. g. guanicoe – V. v. mensalis",
      "guanicoe-vicuna_contact" =
        "L. g. guanicoe – Vicuña contact-zone",
      "guanicoe-vicugna" =
        "L. g. guanicoe – V. v. vicugna"))

#Correct negligible floating-point deviations so that expansion + stability
#sums exactly to 1 for plotting.
dyn.plot <- dyn.plot %>%
  mutate(
    ES_sum = Expansion + Stability,
    Expansion = Expansion / ES_sum,
    Stability = Stability / ES_sum
  ) %>%
  select(
    -ES_sum)

# Define comparison order -------------------------------------------------------
comparison_order <- c(
  "L. g. cacsilensis – V. v. mensalis",
  "L. g. cacsilensis – Vicuña contact-zone",
  "L. g. cacsilensis – V. v. vicugna",
  "Guanaco contact-zone – V. v. mensalis",
  "Guanaco contact-zone – Vicuña contact-zone",
  "Guanaco contact-zone – V. v. vicugna",
  "L. g. guanicoe – V. v. mensalis",
  "L. g. guanicoe – Vicuña contact-zone",
  "L. g. guanicoe – V. v. vicugna")

comparison_labels_expr <- expression(
  italic(L.~g.~cacsilensis) ~ "\u2013" ~ italic(V.~v.~mensalis),
  italic(L.~g.~cacsilensis) ~ "\u2013" ~ "Vicuña contact-zone",
  italic(L.~g.~cacsilensis) ~ "\u2013" ~ italic(V.~v.~vicugna),
  "Guanaco contact-zone" ~ "\u2013" ~ italic(V.~v.~mensalis),
  "Guanaco contact-zone" ~ "\u2013" ~ "Vicuña contact-zone",
  "Guanaco contact-zone" ~ "\u2013" ~ italic(V.~v.~vicugna),
  italic(L.~g.~guanicoe) ~ "\u2013" ~ italic(V.~v.~mensalis),
  italic(L.~g.~guanicoe) ~ "\u2013" ~ "Vicuña contact-zone",
  italic(L.~g.~guanicoe) ~ "\u2013" ~ italic(V.~v.~vicugna))

#ggplot draws the first factor level at the bottom, so the order is reversed.
dyn.plot <- dyn.plot %>%
  mutate(
    Comparison = factor(
      Comparison,
      levels = rev(
        comparison_order)))

# Convert expansion and stability to long format -------------------------------
dyn.long <- dyn.plot %>%
  select(
    Comparison,
    Expansion,
    Stability
  ) %>%
  pivot_longer(
    cols = c(
      Expansion,
      Stability),
    names_to = "Component",
    values_to = "Value"
  ) %>%
  mutate(
    Component = factor(
      Component,
      levels = c(
        "Expansion",
        "Stability")))

# Panel A: Expansion and stability ---------------------------------------------
p_ES <- ggplot(
  dyn.long,
  aes(
    x = Value,
    y = Comparison,
    fill = Component)) +
  
  geom_col(
    width = 0.72,
    position = position_stack(
      reverse = TRUE)) +
  
  geom_text(
    aes(
      label = ifelse(
        Value >= 0.025,
        sprintf(
          "%.2f",
          Value),
        "")),
    position = position_stack(
      vjust = 0.5,
      reverse = TRUE),
    colour = "white",
    size = 3.4) +
  
  scale_fill_manual(
    values = c(
      "Expansion" = "#3b528b",
      "Stability" = "#440154"),
    breaks = c(
      "Expansion",
      "Stability")) +
  
  scale_y_discrete(
    labels = rev(
      comparison_labels_exp)) +
  
  scale_x_continuous(
    breaks = c(
      0,
      0.25,
      0.50,
      0.75,
      1),
    labels = percent_format(
      accuracy = 1),
    expand = c(
      0,
      0)) +
  
  coord_cartesian(
    xlim = c(
      0,
      1)) +
  
  labs(
    x = NULL,
    y = NULL,
    fill = NULL,
    title = "Expansion and Stability") +
  
  theme_bw(
    base_size = 11) +
  
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    
    axis.text.y = element_text(
      colour = "black",
      size = 9),
    
    axis.text.x = element_text(
      colour = "black"),
    
    axis.ticks.y = element_blank(),
    
    plot.title = element_text(
      face = "bold",
      hjust = 0.5),
    
    legend.position = "bottom",
    
    plot.margin = margin(
      t = 5,
      r = 12,
      b = 5,
      l = 5))

#Panel B: Unfilling ------------------------------------------------------------
p_U <- ggplot(
  dyn.plot,
  aes(
    x = Unfilling,
    y = Comparison)) +
  
  geom_segment(
    aes(
      x = 0,
      xend = Unfilling,
      yend = Comparison),
    colour = "grey55",
    linewidth = 0.65) +
  
  geom_point(
    size = 3.2,
    colour = "black") +
  
  geom_text(
    aes(
      label = sprintf(
        "%.2f",
        Unfilling)),
    nudge_x = 0.018,
    hjust = 0,
    size = 3.4) +
  
  scale_x_continuous(
    limits = c(
      0,
      0.68),
    breaks = c(
      0,
      0.2,
      0.4,
      0.6),
    labels = percent_format(
      accuracy = 1),
    expand = c(
      0,
      0)) +
  
  labs(
    x = "Unfilling",
    y = NULL,
    title = "Unfilling") +
  
  theme_bw(
    base_size = 11) +
  
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    
    axis.text.x = element_text(
      colour = "black"),
    
    plot.title = element_text(
      face = "bold",
      hjust = 0.5),
    
    plot.margin = margin(
      t = 5,
      r = 5,
      b = 5,
      l = 12))

#Combine panels ---------------------------------------------------------------
p_dynamic <-
  p_ES +
  plot_spacer() +
  p_U +
  plot_layout(
    widths = c(
      1.45,
      0.10,
      1)) +
  plot_annotation(
    title = "Niche dynamics between guanaco and vicuña groups",
    theme = theme(
      plot.title = element_text(
        size = 13,
        face = "bold")))

# Display figure ---------------------------------------------------------------
p_dynamic

# Sensitivity analysis: non-analogous environments -----------------------------

#This analysis evaluates whether expansion, stability, and unfilling estimates
#change when environmental categories outside the analogous environmental
#space are included.
#
#The standard ecospat dynamic indices (intersection = 0) represent the
#analogue-space estimates. Full-space estimates are reconstructed from
#category_quantity by including both analogous (A) and non-analogous (NA)
#environmental categories.
calc_dyn_sensitivity <- function(
    group1,
    group2,
    bg.pca,
    occ.pca,
    R = 100) {
  
  pair <- build_z(
    group1 = group1,
    group2 = group2,
    bg.pca = bg.pca,
    occ.pca = occ.pca,
    R = R)
  
  dyn <- ecospat.niche.dyn.index(
    z1 = pair$z1,
    z2 = pair$z2,
    intersection = 0)

  cq <- dyn$category_quantity
  
  #Expansion including analogous and non-analogous environments.
  expansion_full <-
    (cq["z2_only_A"] +
        cq["z2_only_NA"]) /
    (cq["z2_only_A"] +
        cq["z2_only_NA"] +
        cq["z2_z1"])

  # Stability including analogous and non-analogous environments.
  stability_full <-
    cq["z2_z1"] /
    (cq["z2_only_A"] +
        cq["z2_only_NA"] +
        cq["z2_z1"])
  
  # Unfilling including analogous and non-analogous environments.
  unfilling_full <-
    (cq["z1_only_A"] +
        cq["z1_only_NA"]) /
    (cq["z1_only_A"] +
        cq["z1_only_NA"] +
        cq["z1_z2"])
  
  data.frame(
    Expansion_analogue = as.numeric(
      dyn$dynamic.index.w["expansion"]),
    Stability_analogue = as.numeric(
      dyn$dynamic.index.w["stability"]),
    Unfilling_analogue = as.numeric(
      dyn$dynamic.index.w["unfilling"]),
    Expansion_full = as.numeric(
      expansion_full),
    Stability_full = as.numeric(
      stability_full),
    Unfilling_full = as.numeric(
      unfilling_full))}

#Run sensitivity analysis for all pairwise comparisons ------------------------
sensitivity.results <- lapply(
  seq_len(nrow(comparisons)),
  function(i) {
    result <- calc_dyn_sensitivity(
      group1 = comparisons$group1[i],
      group2 = comparisons$group2[i],
      bg.pca = bg.pca,
      occ.pca = occ.pca,
      R = 100)
    
    data.frame(
      Comparison = comparisons$Comparison[i],
      result)}) %>%
  bind_rows()

# Calculate differences between full-space and analogue-space estimates.
sensitivity.results <- sensitivity.results %>%
  mutate(
    Delta_Expansion =
      Expansion_full -
      Expansion_analogue,
    Delta_Stability =
      Stability_full -
      Stability_analogue,
    Delta_Unfilling =
      Unfilling_full -
      Unfilling_analogue)

sensitivity.results

#Create rounded table for inspection ------------------------------------------
sensitivity.rounded <- sensitivity.results
sensitivity.rounded[, -1] <- round(
  sensitivity.rounded[, -1],3)

sensitivity.rounded

#Export sensitivity analysis --------------------------------------------------

write.csv(sensitivity.results, "output/niche_dynamics/niche_dynamics_sensitivity.csv",
  row.names = FALSE)