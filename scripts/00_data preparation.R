#-------------------------------------------------------------------------------
# Andrea G. Castillo
# Downloading GBIF data, cleaning occurrence records, and selecting environmental variables
#
# Chapter 3:
# "Where is the guanaco irreplaceable? Lineage-level niche structure and
# functional redundancy in wild South American camelids"
#
# Occurrence data cleaning follows the framework proposed by Zizka et al. (2020).
#
# Date: 14 September 2026

# Load required packages --------------------------------------------------------
library(terra)
library(dplyr)
library(rgbif)
library(corrplot)


# STEP 1: Download occurrence data from GBIF -----------------------------------

# Lama guanicoe
gbif_lama <- occ_search(
  scientificName = "Lama guanicoe",
  fields = c(
    "scientificName",
    "decimalLatitude",
    "decimalLongitude",
    "year"),
  limit = 20000)

write.csv( gbif_lama$data,"data/occurrences/gbif_lama2026.csv", row.names = FALSE)

# Vicugna vicugna
gbif_vicugna <- occ_search(
  scientificName = "Vicugna vicugna",
  fields = c(
    "scientificName",
    "decimalLatitude",
    "decimalLongitude",
    "year"),
  limit = 20000)

write.csv(gbif_vicugna$data, "data/occurrences/gbif_vicugna2025.csv", row.names = FALSE)


# GBIF references:
# Lama guanicoe: GBIF.org (22 May 2026). GBIF Occurrence Download.
# URl: https://doi.org/10.15468/dl.j3c6qn
#
# Vicugna vicugna: GBIF.org (4 December 2025). GBIF Occurrence Download.
# URL: https://doi.org/10.15468/dl.33q2t7
#
# NOTE: The GBIF DOI references above correspond to archived GBIF occurrence
# downloads used in the study. Queries performed through occ_search() may return 
# different records if the GBIF database is updated.


# STEP 2: Clean occurrence data ------------------------------------------------

# Read compiled occurrence datasets
lama <- read.csv( "data/occurrences/All_Lama_mayo2026.csv",header = TRUE)

vicu <- read.csv( "data/occurrences/Alldata_vicugna2025.csv", header = TRUE)


# i) Remove records with missing longitude, latitude, or year ------------------
lama_1 <- lama[!is.na(lama$decimalLon), ]
lama_2 <- lama_1[!is.na(lama_1$decimalLat), ]
lama_3 <- lama_2[!is.na(lama_2$year), ]

vicu_1 <- vicu[!is.na(vicu$decimalLon), ]
vicu_2 <- vicu_1[!is.na(vicu_1$decimalLat), ]
vicu_3 <- vicu_2[!is.na(vicu_2$year), ]

# ii) Remove duplicate records based on coordinates ----------------------------
lama_4 <- lama_3[!duplicated(lama_3[, 3:4]), ]
vicu_4 <- vicu_3[!duplicated(vicu_3[, 3:4]), ]


# iii) Remove records with zero-valued coordinates -----------------------------
lama_5 <- lama_4[lama_4$decimalLon != 0 & lama_4$decimalLat != 0,]
vicu_5 <- vicu_4[vicu_4$decimalLon != 0 & vicu_4$decimalLat != 0,]

# iv) Exclude records collected before 1980 ------------------------------------
lama_6 <- lama_5[lama_5$year >= 1980, ]
vicu_6 <- vicu_5[vicu_5$year >= 1980, ]

# v) Exclude records without decimal coordinate precision ----------------------
lama_7 <- lama_6[(lama_6$decimalLon - floor(lama_6$decimalLon)) != 0,]
lama_8 <- lama_7[(lama_7$decimalLat - floor(lama_7$decimalLat)) != 0,]

vicu_7 <- vicu_6[(vicu_6$decimalLon - floor(vicu_6$decimalLon)) != 0,]
vicu_8 <- vicu_7[(vicu_7$decimalLat - floor(vicu_7$decimalLat)) != 0,]

# Save cleaned occurrence records ---------------------------------------------
write.csv(lama_8,"data/occurrences/lama_1cleaned_occs.csv", row.names = FALSE)
write.csv(vicu_8, "data/occurrences/vicu_1cleaned_occs.csv", row.names = FALSE)

# vi) Check whether occurrence records fall within the South America shapefile using ArcGIS.
# vii) Spatially thin occurrence records to one record per grid cell using ArcGIS 
# Lama guanicoe
lama_grid <- read.csv("data/occurrences/lama_grid_intersect.csv", header = TRUE)
lama_thinned <- lama_grid[!duplicated(lama_grid$grid_id)]
write.csv(lama_thinned, "data/occurrences/lama_2cleaned_occs.csv",row.names = FALSE)

# Vicugna vicugna
vicu_grid <- read.csv("data/occurrences/vicu_grid_intersect.csv",header = TRUE)
vicu_thinned <- vicu_grid[!duplicated(vicu_grid$grid_id)]
write.csv(vicu_thinned, "data/occurrences/vicu_2cleaned_occs.csv", row.names = FALSE)


# STEP 3: Select environmental variables ---------------------------------------

# Read environmental layers
# NOTE: Environmental layers are stored locally under data/environmental/.
# Because of file size and data-distribution considerations, these layers
# may not be hosted directly in the GitHub repository. Data sources and
# instructions for obtaining the environmental layers will be provided
# in the repository README.

bio1 <- rast("data/environmental/bio1.asc")
bio6 <- rast("data/environmental/bio6.asc")
bio7 <- rast("data/environmental/bio7.asc")
bio11 <- rast("data/environmental/bio11.asc")
bio12 <- rast("data/environmental/bio12.asc")
bio14 <- rast("data/environmental/bio14.asc")
bio17 <- rast("data/environmental/bio17.asc")
CMI_min <- rast("data/environmental/cmi_min.asc")
CMI_max <- rast("data/environmental/cmi_max.asc")
CMI_mean <- rast("data/environmental/cmi_mean.asc")


# Create environmental raster stack
stck_var <- c(
  bio1,
  bio6,
  bio7,
  bio11,
  bio12,
  bio14,
  bio17,
  CMI_min,
  CMI_max,
  CMI_mean)


# Vicugna vicugna --------------------------------------------------------------
vicu.occs<-vicu_thinned 

# Extract environmental values at occurrence locations
vicu_env <- terra::extract(
  stck_var,
  vicu.occs,
  ID = FALSE,
  xy = TRUE) %>%
  select(x, y, everything()) %>%
  rename(
    decimalLon = x,
    decimalLat = y)

# Remove unnecessary columns
vicu_env$year <- NULL
vicu_env$source <- NULL
vicu_env$ID <- NULL
vicu_env$locality <- NULL

# Retain environmental variables only
vicu_env <- vicu_env[, -c(1:3)]

# Replace NoData values with NA and remove missing observations
vicu_env[vicu_env == -9999.0] <- NA
vicu_env <- na.omit(vicu_env)

# Calculate Pearson correlation matrix
corr_vicu <- cor(vicu_env)

# Export correlation matrix
write.table(corr_vicu, "output/correlation/corr_vicu_2026.csv",sep = ";")


# Plot correlation matrix
col_fun <- colorRampPalette(
  c( "#BB4444",
     "#EE9988",
     "#FFFFFF",
     "#77AADD",
     "#4477AA"))

corrplot(
  corr_vicu,
  method = "shade",
  shade.col = NA,
  tl.col = "black",
  tl.srt = 45,
  col = col_fun(200),
  addCoef.col = "black",
  addcolorlabel = "no",
  order = "AOE",
  type = "lower")


# Perform PCA
pca_vicu <- prcomp( vicu_env, scale. = TRUE, center = TRUE)
pca_vicu
summary(pca_vicu)
plot(pca_vicu)
biplot(pca_vicu)


# Lama guanicoe ----------------------------------------------------------------
lama.occs<-lama_thinned 

# Extract environmental values at occurrence locations
lama_env <- terra::extract(
  stck_var,
  lama.occs,
  ID = FALSE,
  xy = TRUE) %>%
  select(x, y, everything()) %>%
  rename(
    decimalLon = x,
    decimalLat = y)

# Remove unnecessary columns
lama_env$year <- NULL
lama_env$source <- NULL
lama_env$ID <- NULL
lama_env$locality <- NULL

# Retain environmental variables only
lama_env <- lama_env[, -c(1:3)]

# Replace NoData values with NA and remove missing observations
lama_env[lama_env == -9999.0] <- NA
lama_env <- na.omit(lama_env)

# Calculate Pearson correlation matrix
corr_lama <- cor(lama_env)

# Export correlation matrix
write.table(corr_lama, "output/correlation/corr_lama_2026.csv",sep = ";")

# Plot correlation matrix
corrplot(
  corr_lama,
  method = "shade",
  shade.col = NA,
  tl.col = "black",
  tl.srt = 45,
  col = col_fun(200),
  addCoef.col = "black",
  addcolorlabel = "no",
  order = "AOE",
  type = "lower")


# Perform PCA
pca_lama <- prcomp(lama_env, scale. = TRUE, center = TRUE)
pca_lama
summary(pca_lama)
plot(pca_lama)
biplot(pca_lama)