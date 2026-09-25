# ------------------------------------------------------------------------------
# Andrea G. Castillo
# Niche similarity analyses using "humboldt" package (Brown & Carnaval, 2019)
#
# Chapter 3:
# "Where is the guanaco irreplaceable? Lineage-level niche structure and
# functional redundancy in wild South American camelids"
#
# Pairwise niche equivalency and background similarity tests among guanaco
# and vicuña lineages and contact-zone populations (https://doi.org/10.21425/F5FBG44158)
#
# Date: 14 September 2026

#Install the necessary packages-------------------------------------------------
install.packages("devtools")
install_github("jasonleebrown/humboldt")

#Load required packages --------------------------------------------------------
library(humboldt)
library(dplyr)
library(ggplot2)

#Environmental variables used in all analyses ----------------------------------
vars_env <- c("bio6", "bio14", "cmi_max")


#1. SET Humboldt PARAMETERS ----------------------------------------------------
#NOTE:These parameters are defined here to ensure that the same settings are used
#consistently across all pairwise comparisons.
#The standard analysis uses 100 randomisations. The number of randomisations
#can be increased when greater precision around a test statistic or P-value
#is required (e.g. rep=1000)
n_rep <- 100
env_resolution <- 0.04491576
kernel_smoothing <- 1
grid_resolution <- 100
density_threshold <- 0.001
n_cores <- 1


#2.READ DATA -------------------------------------------------------------------
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

#3. PREPARE DATA FOR HUMBOLDT --------------------------------------------------
#Humboldt requires occurrence data to contain:
#  sp = group identifier
# x  = longitude
#  y  = latitude
#
#The original occurrence matrices use decimalLon and decimalLat. These columns
#are renamed here so that the original input files do not need to be modified
#manually before running the analysis.

prepare_humboldt_occ <- function(data, group_name) {
  data %>%
    mutate(sp = group_name) %>%
    rename(
      x = decimalLon,
      y = decimalLat
    ) %>%
    select(sp, x, y, everything())}

cac.occs <- prepare_humboldt_occ(cac.occs, "cacsilensis")
lhyb.occs <- prepare_humboldt_occ(lhyb.occs, "lhyb")
gua.occs <- prepare_humboldt_occ(gua.occs, "guanicoe")
men.occs <- prepare_humboldt_occ(men.occs, "mensalis")
vhyb.occs <- prepare_humboldt_occ(vhyb.occs, "vhyb")
vic.occs <- prepare_humboldt_occ(vic.occs, "vicugna")

# Remove records with incomplete environmental information
remove_env_na <- function(data) {
  data[complete.cases(data[, vars_env]), ]}

cac.occs <- remove_env_na(cac.occs)
lhyb.occs <- remove_env_na(lhyb.occs)
gua.occs <- remove_env_na(gua.occs)
men.occs <- remove_env_na(men.occs)
vhyb.occs <- remove_env_na(vhyb.occs)
vic.occs <- remove_env_na(vic.occs)

Mcac.env <- remove_env_na(Mcac.env)
Mlhyb.env <- remove_env_na(Mlhyb.env)
Mgua.env <- remove_env_na(Mgua.env)
Mmen.env <- remove_env_na(Mmen.env)
Mvhyb.env <- remove_env_na(Mvhyb.env)
Mvic.env <- remove_env_na(Mvic.env)

# Check retained sample sizes
data.frame(
  group = c(
    "cacsilensis",
    "lhyb",
    "guanicoe",
    "mensalis",
    "vhyb",
    "vicugna"),
  occurrences = c(
    nrow(cac.occs),
    nrow(lhyb.occs),
    nrow(gua.occs),
    nrow(men.occs),
    nrow(vhyb.occs),
    nrow(vic.occs)),
  background = c(
    nrow(Mcac.env),
    nrow(Mlhyb.env),
    nrow(Mgua.env),
    nrow(Mmen.env),
    nrow(Mvhyb.env),
    nrow(Mvic.env)))


#4. ORGANISE GROUP DATA --------------------------------------------------------
#Using named lists prevents ambiguity between the guanaco and vicuña
#contact-zone populations and avoids repeatedly defining the same objects.
occ_list <- list(
  cacsilensis = cac.occs,
  lhyb = lhyb.occs,
  guanicoe = gua.occs,
  mensalis = men.occs,
  vhyb = vhyb.occs,
  vicugna = vic.occs)

env_list <- list(
  cacsilensis = Mcac.env,
  lhyb = Mlhyb.env,
  guanicoe = Mgua.env,
  mensalis = Mmen.env,
  vhyb = Mvhyb.env,
  vicugna = Mvic.env)


#5. DEFINE PAIRWISE COMPARISONS -----------------------------------------------
#Six within-species comparisons and nine between-species comparisons.

comparisons <- data.frame(
  group1 = c(
    # Guanaco
    "cacsilensis",
    "cacsilensis",
    "guanicoe",
    
    # Vicuña
    "mensalis",
    "mensalis",
    "vicugna",
    
    # Between species
    "cacsilensis",
    "cacsilensis",
    "cacsilensis",
    "lhyb",
    "lhyb",
    "lhyb",
    "guanicoe",
    "guanicoe",
    "guanicoe"),
  
  group2 = c(
    # Guanaco
    "guanicoe",
    "lhyb",
    "lhyb",
    
    # Vicuña
    "vicugna",
    "vhyb",
    "vhyb",
    
    # Between species
    "mensalis",
    "vhyb",
    "vicugna",
    "mensalis",
    "vhyb",
    "vicugna",
    "mensalis",
    "vhyb",
    "vicugna"),
  
  stringsAsFactors = FALSE)

comparisons


#6. FUNCTION FOR PAIRWISE HUMBOLDT ANALYSES ------------------------------------
run_humboldt_comparison <- function(group1, group2) {

  message(
    "\n------------------------------------------------------------\n",
    "Running Humboldt analysis: ", group1, " vs ", group2,
    "\n------------------------------------------------------------")
  
  # --------------------------------------------------------------------------
  # A. Build environmental space
  # --------------------------------------------------------------------------
  #
  # The suffix "_NOT" is retained to preserve consistency with the terminology
  # used during the original analyses and manuscript preparation.
  
  g2e_NOT <- humboldt.g2e(
    env1 = env_list[[group1]],
    env2 = env_list[[group2]],
    sp1 = occ_list[[group1]],
    sp2 = occ_list[[group2]],
    
    reduce.env = 0,
    reductype = "PCA",
    non.analogous.environments = "YES",
    env.trim = FALSE,
    
    e.var = 3:5,
    col.env = 3:5,
    
    env.reso = env_resolution,
    kern.smooth = kernel_smoothing,
    R = grid_resolution,
    
    run.silent = FALSE)

  # --------------------------------------------------------------------------
  # B. Extract scores from environmental space
  # --------------------------------------------------------------------------
  scores.env1 <- g2e_NOT$scores.env1[, 1:2]
  scores.env2 <- g2e_NOT$scores.env2[, 1:2]
  
  scores.env12 <- rbind(
    scores.env1,
    scores.env2)
  
  scores.sp1 <- g2e_NOT$scores.sp1[, 1:2]
  scores.sp2 <- g2e_NOT$scores.sp2[, 1:2]
  
  # --------------------------------------------------------------------------
  # C. Construct gridded niche spaces
  # --------------------------------------------------------------------------
  
  z.sp1 <- humboldt.grid.espace(
    scores.env12,
    scores.env1,
    scores.sp1,
    kern.smooth = kernel_smoothing,
    R = grid_resolution)
  
  z.sp2 <- humboldt.grid.espace(
    scores.env12,
    scores.env2,
    scores.sp2,
    kern.smooth = kernel_smoothing,
    R = grid_resolution)
  
  # Environmental availability grids
  # Retained because these objects were part of the original analytical
  # workflow and may be useful for diagnostics.
  
  z.env1 <- humboldt.grid.espace(
    scores.env12,
    scores.env1,
    scores.env1,
    kern.smooth = kernel_smoothing,
    R = grid_resolution)
  
  z.env2 <- humboldt.grid.espace(
    scores.env12,
    scores.env2,
    scores.env2,
    kern.smooth = kernel_smoothing,
    R = grid_resolution)
  
  # --------------------------------------------------------------------------
  # D. Observed niche similarity
  # --------------------------------------------------------------------------
  
  niche.sim.correct <- humboldt.niche.similarity(
    z.sp1,
    z.sp2,
    correct.env = TRUE,
    nae = "YES")
  
  niche.sim.uncorrect <- humboldt.niche.similarity(
    z.sp1,
    z.sp2,
    correct.env = FALSE,
    nae = "YES")
  
  # --------------------------------------------------------------------------
  # E. Niche equivalency test
  # --------------------------------------------------------------------------
  equivalency <- humboldt.equivalence.stat(
    z1 = z.sp1,
    z2 = z.sp2,
    
    rep = n_rep,
    
    correct.env = TRUE,
    kern.smooth = kernel_smoothing,
    nae = "YES",
    thresh.espace.z = density_threshold,
    
    run.silent.equ = FALSE,
    ncores = n_cores)
  
  # --------------------------------------------------------------------------
  # F. Background similarity tests
  # --------------------------------------------------------------------------
  # Background tests are directional:
  # sim.dir = 1 : group1 -> group2
  # sim.dir = 2 : group2 -> group1
  
  background_1_to_2 <- humboldt.background.stat(
    g2e = g2e_NOT,
    
    rep = n_rep,
    sim.dir = 1,
    
    env.reso = env_resolution,
    kern.smooth = kernel_smoothing,
    R = grid_resolution,
    
    correct.env = FALSE,
    thresh.espace.z = density_threshold,
    force.equal.sample = FALSE,
    
    run.silent.bak = FALSE,
    ncores = n_cores)
  
  background_2_to_1 <- humboldt.background.stat(
    g2e = g2e_NOT,
    
    rep = n_rep,
    sim.dir = 2,
    
    env.reso = env_resolution,
    kern.smooth = kernel_smoothing,
    R = grid_resolution,
    
    correct.env = FALSE,
    thresh.espace.z = density_threshold,
    force.equal.sample = FALSE,
    
    run.silent.bak = FALSE,
    ncores = n_cores)
  
  # --------------------------------------------------------------------------
  # G. Return complete analysis
  # --------------------------------------------------------------------------
  list(
    group1 = group1,
    group2 = group2,
    
    g2e_NOT = g2e_NOT,
    
    z.sp1 = z.sp1,
    z.sp2 = z.sp2,
    
    z.env1 = z.env1,
    z.env2 = z.env2,
    
    niche.sim.correct = niche.sim.correct,
    niche.sim.uncorrect = niche.sim.uncorrect,
    
    equivalency = equivalency,
    
    background_1_to_2 = background_1_to_2,
    background_2_to_1 = background_2_to_1)}


#7. RUN ALL PAIRWISE COMPARISONS ----------------------------------------------
humboldt_results <- vector(
  mode = "list",
  length = nrow(comparisons))

names(humboldt_results) <- paste(
  comparisons$group1,
  comparisons$group2,
  sep = "_vs_")

for (i in seq_len(nrow(comparisons))) {
  g1 <- comparisons$group1[i]
  g2 <- comparisons$group2[i]
  humboldt_results[[i]] <- run_humboldt_comparison(
    group1 = g1,
    group2 = g2)}

#8. EXTRACT PAIRWISE NICHE SIMILARITY RESULTS ----------------------------------
#NOTE: Observed niche similarity is extracted directly from each analysis object.
#These values are subsequently used to construct the results table and figure,
#avoiding manually entered D or I values.

pairwise_similarity <- do.call(
  rbind,
  lapply(humboldt_results, function(x) {
    
    data.frame(
      group1 = x$group1,
      group2 = x$group2,
      
      D = x$niche.sim.correct$D,
      I = x$niche.sim.correct$I,
      
      D_uncorrected = x$niche.sim.uncorrect$D,
      I_uncorrected = x$niche.sim.uncorrect$I,

      stringsAsFactors = FALSE)}))

rownames(pairwise_similarity) <- NULL

pairwise_similarity


# 9. EXTRACT DIRECTIONAL BACKGROUND TEST RESULTS -------------------------------
background_results <- do.call(
  rbind,
  lapply(humboldt_results, function(x) {
    
    data.frame(
      comparison = paste(x$group1, x$group2, sep = " - "),
    
      direction = c(
        paste(x$group1, "->", x$group2),
        paste(x$group2, "->", x$group1)),
      
      observed_D = c(
        x$background_1_to_2$obs$D,
        x$background_2_to_1$obs$D),
      
      observed_I = c(
        x$background_1_to_2$obs$I,
        x$background_2_to_1$obs$I),
      
      p_D = c(
        x$background_1_to_2$p.D,
        x$background_2_to_1$p.),
      
      p_I = c(
        x$background_1_to_2$p.I,
        x$background_2_to_1$p.I),
      
      stringsAsFactors = FALSE)}))

rownames(background_results) <- NULL

background_results


#10. INSPECT EQUIVALENCY TESTS -------------------------------------------------
# Equivalency-test objects are retained in full within humboldt_results.
# Example:
# humboldt_results[["cacsilensis_vs_guanicoe"]]$equivalency
# All equivalency tests use n_rep = 100. If a particular result requires
# greater numerical precision, n_rep can be increased and that comparison
# rerun separately.

equivalency_results <- lapply(
  humboldt_results,
  function(x) x$equivalency)

names(equivalency_results)


#11. NDT ANALYSIS: L. g. guanicoe vs V. v. mensalis ----------------------------
#This additional analysis follows the NDT configuration used for the
#guanicoe-mensalis comparison in the original workflow.

gua_men_NDT <- humboldt.g2e(
  env1 = Mgua.env,
  env2 = Mmen.env,
  
  sp1 = gua.occs,
  sp2 = men.occs,
  
  reduce.env = 2,
  reductype = "PCA",
  non.analogous.environments = "NO",
  env.trim = FALSE,
  
  e.var = 3:5,
  col.env = 3:5,
  
  env.reso = env_resolution,
  kern.smooth = kernel_smoothing,
  R = grid_resolution,
  
  run.silent = FALSE)

#Extract scores from the shared environmental space
scores.env1.NDT <- gua_men_NDT$scores.env1[, 1:2]
scores.env2.NDT <- gua_men_NDT$scores.env2[, 1:2]

scores.env12.NDT <- rbind(
  scores.env1.NDT,
  scores.env2.NDT)

scores.sp1.NDT <- gua_men_NDT$scores.sp1[, 1:2]
scores.sp2.NDT <- gua_men_NDT$scores.sp2[, 1:2]

#Construct gridded niche spaces
z.gua.NDT <- humboldt.grid.espace(
  scores.env12.NDT,
  scores.env1.NDT,
  scores.sp1.NDT,
  kern.smooth = kernel_smoothing,
  R = grid_resolution)

z.men.NDT <- humboldt.grid.espace(
  scores.env12.NDT,
  scores.env2.NDT,
  scores.sp2.NDT,
  kern.smooth = kernel_smoothing,
  R = grid_resolution)

#Observed similarity
niche.sim.NDT <- humboldt.niche.similarity(
  z.gua.NDT,
  z.men.NDT,
  correct.env = TRUE,
  nae = "NO")

niche.sim.NDT$D
niche.sim.NDT$I

#Equivalency test
equiv.gua_men_NDT <- humboldt.equivalence.stat(
  z1 = z.gua.NDT,
  z2 = z.men.NDT,
  rep = n_rep,
  correct.env = TRUE,
  kern.smooth = kernel_smoothing,
  nae = "NO",
  thresh.espace.z = density_threshold,
  run.silent.equ = FALSE,
  ncores = n_cores)
equiv.gua_men_NDT

#Background test: guanicoe -> mensalis
gua_to_men_NDT <- humboldt.background.stat(
  g2e = gua_men_NDT,
  rep = n_rep,
  sim.dir = 1,
  env.reso = env_resolution,
  kern.smooth = kernel_smoothing,
  R = grid_resolution,
  correct.env = FALSE,
  thresh.espace.z = density_threshold,
  force.equal.sample = FALSE,
  run.silent.bak = FALSE,
  ncores = n_cores)

#Background test: mensalis -> guanicoe
men_to_gua_NDT <- humboldt.background.stat(
  g2e = gua_men_NDT,
  rep = n_rep,
  sim.dir = 2,
  env.reso = env_resolution,
  kern.smooth = kernel_smoothing,
  R = grid_resolution,
  correct.env = FALSE,
  thresh.espace.z = density_threshold,
  force.equal.sample = FALSE,
  run.silent.bak = FALSE,
  ncores = n_cores)

# Inspect NDT background-test results

gua_to_men_NDT$obs
gua_to_men_NDT$p.D
gua_to_men_NDT$p.I

men_to_gua_NDT$obs
men_to_gua_NDT$p.D
men_to_gua_NDT$p.I


# 12. EXPORT RESULTS -----------------------------------------------------------
dir.create ("output/humboldt", recursive = TRUE, showWarnings = FALSE)

write.csv (pairwise_similarity, "output/humboldt/humboldt_pairwise_similarity.csv",
  row.names = FALSE)
write.csv (background_results, "output/humboldt/humboldt_background_tests.csv",
  row.names = FALSE)

#Save complete R objects so equivalency distributions, background
#randomisations, environmental spaces, and other diagnostics remain available.
saveRDS (humboldt_results,"output/humboldt/humboldt_pairwise_results.rds")

saveRDS(
  list(
    g2e = gua_men_NDT,
    similarity = niche.sim.NDT,
    equivalency = equiv.gua_men_NDT,
    guanicoe_to_mensalis = gua_to_men_NDT,
    mensalis_to_guanicoe = men_to_gua_NDT),
  "output/humboldt/humboldt_gua_men_NDT.rds")


#13. SCHOENER'S D MATRIX -------------------------------------------------------
#Order used throughout the manuscript
groups <- c(
  "cacsilensis",
  "lhyb",
  "guanicoe",
  "mensalis",
  "vhyb",
  "vicugna")

#Place all pairwise comparisons in the lower triangle
D_matrix <- pairwise_similarity %>%
  rowwise() %>%
  mutate(
    pos1 = match(group1, groups),
    pos2 = match(group2, groups),
    
    x = ifelse(pos1 < pos2, group1, group2),
    y = ifelse(pos1 < pos2, group2, group1)
  ) %>%
  ungroup() %>%
  select(x, y, value = D)

#Add diagonal (D = 1)
diagonal <- data.frame(
  x = groups,
  y = groups,
  value = 1)

D_plot_data <- bind_rows(
  D_matrix,
  diagonal)

# Factor order
D_plot_data$x <- factor(
  D_plot_data$x,
  levels = groups)

D_plot_data$y <- factor(
  D_plot_data$y,
  levels = rev(groups))

#Labels used in the manuscript
group_labels <- c(
  cacsilensis = expression(italic(L.~g.~cacsilensis)),
  lhyb = "Guanaco contact-zone",
  guanicoe = expression(italic(L.~g.~guanicoe)),
  mensalis = expression(italic(V.~v.~mensalis)),
  vhyb = "Vicuña contact-zone",
  vicugna = expression(italic(V.~v.~vicugna)))

#Plot
D_plot <- ggplot(
  D_plot_data,
  aes(x = x, y = y, fill = value)) +
  
  geom_tile(
    colour = "white",
    linewidth = 0.7) +
  
  geom_text(
    aes(label = sprintf("%.2f", value)),
    size = 3.8) +
  
  scale_fill_gradient(
    low = "white",
    high = "black",
    limits = c(0, 1),
    name = "Schoener's D") +
  
  scale_x_discrete(
    labels = group_labels,
    drop = FALSE) +
  
  scale_y_discrete(
    labels = group_labels,
    drop = FALSE) +
  
  coord_fixed() +
  
  labs(
    x = NULL,
    y = NULL ) +
  
  theme_classic() +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1),
    
    axis.ticks = element_blank(),
    
    legend.position = "right")

D_plot


#Save figure
dir.create ("figures/manuscript", recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = "figures/manuscript/humboldt_schoeners_D_matrix.pdf",
  plot = D_plot,
  width = 8,
  height = 7,
  units = "in")


#14. SESSION INFORMATION -------------------------------------------------------
sessionInfo()
