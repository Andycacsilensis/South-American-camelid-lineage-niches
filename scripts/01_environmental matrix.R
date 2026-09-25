# ------------------------------------------------------------------------------
# Andrea G. Castillo
# Construction of occurrence and background environmental matrices
#
# Chapter 3:
# "Where is the guanaco irreplaceable? Lineage-level niche structure and
# functional redundancy in wild South American camelids"
#
# Code based on the methodological framework proposed by Jiménez and Soberón (2022).
#
# Date: 14 September 2026

# Load required packages --------------------------------------------------------
library(terra)
library(dplyr)

# Load environmental layers ----------------------------------------------------
# Note: Environmental layers are not hosted in this repository.
# Data sources and download instructions are provided in the README.
bio6 <- rast("data/environmental/bio6.asc")
bio14 <- rast("data/environmental/bio14.asc")
CMI_max <- rast("data/environmental/cmi_max.asc")

stck_var <- c(
  bio6,
  bio14,
  CMI_max)


# Lama guanicoe lineages and contact-zone populations --------------------------

#Read occurrence data
cac.occs <- read.csv("data/occurrences/Lama_g_cacsilensis.csv", header=T)[,2:3] 
ghyb.occs <- read.csv("data/occurrences/Lama_g_hybrid.csv", header=T)[,2:3]    
gua.occs <- read.csv("data/occurrences/Lama_g_guanicoe.csv", header=T)[,2:3]    

#Read accessible-area (M) polygons
M.cac <- vect("data/spatial/rm_lgcacsilensis.shp")
M.ghyb <- vect( "data/spatial/rm_lhybrid.shp")
M.gua <- vect( "data/spatial/rm_lgguanicoe.shp")

#Extract environmental values at occurrence locations
cac.occs <- terra::extract(
  stck_var,
  cac.occs,
  ID = FALSE,
  xy = TRUE
) %>%
  select(x, y, everything()) %>%
  rename(
    decimalLon = x,
    decimalLat = y)

ghyb.occs <- terra::extract(
  stck_var,
  ghyb.occs,
  ID = FALSE,
  xy = TRUE
) %>%
  select(x, y, everything()) %>%
  rename(
    decimalLon = x,
    decimalLat = y)

gua.occs <- terra::extract(
  stck_var,
  gua.occs,
  ID = FALSE,
  xy = TRUE
) %>%
  select(x, y, everything()) %>%
  rename(
    decimalLon = x,
    decimalLat = y)

#Export occurrence environmental matrices
write.csv(cac.occs,  "data/environmental_matrices/cacsilensis_occurrences_lonlat_bio_cmimax_2026.csv", 
          row.names = F)
write.csv(ghyb.occs, "data/environmental_matrices/lhybrid_occurrences_lonlat_bio_cmimax_2026.csv",
          row.names = F)
write.csv(gua.occs, "data/environmental_matrices/guanicoe_occurrences_lonlat_bio_cmimax_2026.csv", 
          row.names = F)

#Generate random background points within each accessible area (M)
#Guanaco groups use 10,000 background points per accessible area.
Mrp.cac <- spatSample(M.cac, size = 10000, method = "random")
Mrp.ghyb <- spatSample(M.ghyb, size = 10000, method = "random")
Mrp.gua <- spatSample(M.gua, size = 10000, method = "random")

#Extract environmental values at background locations
Mrp.cac <- terra::extract(
  stck_var,
  Mrp.cac,
  ID = FALSE,
  xy = TRUE
) %>%
  select(x, y, everything()) %>%
  rename(
    decimalLon = x,
    decimalLat = y)

Mrp.ghyb <- terra::extract(
  stck_var,
  Mrp.ghyb,
  ID = FALSE,
  xy = TRUE
) %>%
  select(x, y, everything()) %>%
  rename(
    decimalLon = x,
    decimalLat = y)

Mrp.gua <- terra::extract(
  stck_var,
  Mrp.gua,
  ID = FALSE,
  xy = TRUE
) %>%
  select(x, y, everything()) %>%
  rename(
    decimalLon = x,
    decimalLat = y)

#Export background environmental matrices
write.csv(Mrp.cac, "data/environmental_matrices/cacsilensis_10K_Mpoints_lonlat_bio_cmimax_2026.csv",
          row.names = F)
write.csv(Mrp.ghyb,  "data/environmental_matrices/lhybrid_10K_Mpoints_lonlat_bio_cmimax_2026.csv",
          row.names = F)
write.csv(Mrp.gua, "data/environmental_matrices/guanicoe_10K_Mpoints_lonlat_bio_cmimax_2026.csv",
          row.names = F)


# Vicugna vicugna lineages and contact-zone populations ------------------------
#Read occurrence data
men.occs <- read.csv("data/occurrences/vv_mensalis_2026.csv", header=T)[,2:3] 
vhyb.occs <- read.csv("data/occurrences/v_hybrid_2026.csv", header=T)[,2:3]   
vic.occs <- read.csv("data/occurrences/vv_vicugna_2026.csv", header=T)[,2:3]

#Read accessible-area (M) polygons
M.men <- vect("data/spatial/rm_vmensalis.shp")
M.vhyb <- vect("data/spatial/rm_vhybrid.shp")
M.vic <- vect( "data/spatial/rm_vvicugna.shp")

#Extract environmental values at occurrence locations
men.occs <- terra::extract(
  stck_var,
  men.occs,
  ID = FALSE,
  xy = TRUE
) %>%
  select(x, y, everything()) %>%
  rename(
    decimalLon = x,
    decimalLat = y)

vic.occs <- terra::extract(
  stck_var,
  vic.occs,
  ID = FALSE,
  xy = TRUE
) %>%
  select(x, y, everything()) %>%
  rename(
    decimalLon = x,
    decimalLat = y)

vhyb.occs <- terra::extract(
  stck_var,
  hyb.occs,
  ID = FALSE,
  xy = TRUE
) %>%
  select(x, y, everything()) %>%
  rename(
    decimalLon = x,
    decimalLat = y)

#Export occurrence environmental matrices
write.csv(men.occs, "data/environmental_matrices/mensalis_occurrences_lonlat_bio_cmimax_2026.csv", 
          row.names = F)
write.csv(hyb.occs,  "data/environmental_matrices/vhybrid_occurrences_lonlat_bio_cmimax_2026.csv",  
          row.names = F)
write.csv(vic.occs,  "data/environmental_matrices/vicugna_occurrences_lonlat_bio_cmimax_2026.csv",
          row.names = F)

#Generate random background points within each accessible area (M)
#Vicuña groups use 5,000 background points per accessible area.
Mrp.men <- spatSample(M.men, size = 5000, method = "random")
Mrp.vhyb <- spatSample(M.vhyb, size = 5000, method = "random")
Mrp.vic <- spatSample(M.vic, size = 5000, method = "random")

#Extract environmental values at background locations
Mrp.men <- terra::extract(
  stck_var,
  Mrp.men,
  ID = FALSE,
  xy = TRUE
) %>%
  select(x, y, everything()) %>%
  rename(
    decimalLon = x,
    decimalLat = y)

Mrp.vhyb <- terra::extract(
  stck_var,
  Mrp.vhyb,
  ID = FALSE,
  xy = TRUE
) %>%
  select(x, y, everything()) %>%
  rename(
    decimalLon = x,
    decimalLat = y)

Mrp.vic <- terra::extract(
  stck_var,
  Mrp.vic,
  ID = FALSE,
  xy = TRUE
) %>%
  select(x, y, everything()) %>%
  rename(
    decimalLon = x,
    decimalLat = y)

#Export background environmental matrices
write.csv(Mrp.men_env,"data/environmental_matrices/mensalis_5K_Mpoints_lonlat_bio_cmimax_2026.csv",
  row.names = FALSE)
write.csv(Mrp.vhyb_env,"data/environmental_matrices/vhybrid_5K_Mpoints_lonlat_bio_cmimax_2026.csv",
  row.names = FALSE)
write.csv(Mrp.vic_env,"data/environmental_matrices/vicugna_5K_Mpoints_lonlat_bio_cmimax_2026.csv",
  row.names = FALSE)




