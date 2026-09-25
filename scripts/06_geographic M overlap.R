# ------------------------------------------------------------------------------
# Andrea G. Castillo
# Geographic overlap among accessible areas (M)
#
# Chapter 3:
# "Where is the guanaco irreplaceable? Lineage-level niche structure and
# functional redundancy in wild South American camelids"
#
# This script quantifies pairwise geographic overlap between the accessible
# areas (M) of guanaco and vicuña groups. Overlap is calculated as the shared
# area (km2) and as the percentage of each group's accessible area represented
# by the intersection.
#
# Date: 14 September 2026

#Load required packages --------------------------------------------------------
library(sf)
library(dplyr)


#Read accessible-area (M) polygons --------------------------------------------
# Lama guanicoe
M.cac  <- st_read("data/spatial/rm_lgcacsilensis.shp")
M.lhyb <- st_read("data/spatial/rm_lhybrid.shp")
M.gua  <- st_read("data/spatial/rm_lgguanicoe.shp")
# Vicugna vicugna
M.men  <- st_read("data/spatial/rm_vmensalis.shp")
M.vhyb <- st_read("data/spatial/rm_vhybrid.shp")
M.vic  <- st_read("data/spatial/rm_vvicugna.shp")

#Check coordinate reference systems -------------------------------------------
#All input polygons are expected to use WGS84 (EPSG:4326).
st_crs(M.cac)
st_crs(M.lhyb)
st_crs(M.gua)
st_crs(M.men)
st_crs(M.vhyb)
st_crs(M.vic)

#Transform polygons to an equal-area projection -------------------------------
#NOTE: EPSG:4326 uses geographic coordinates in degrees and is therefore not
#appropriate for calculating continental-scale areas in km2.
#All polygons are transformed to an Albers Equal Area projection before
#calculating areas and geographic intersections.

crs_aea <-
  "+proj=aea +lat_1=-5 +lat_2=-42 +lat_0=-22 +lon_0=-60 +datum=WGS84 +units=m +no_defs"

# Lama guanicoe
M.cac_aea  <- st_transform (M.cac, crs_aea)
M.lhyb_aea <- st_transform (M.lhyb, crs_aea)
M.gua_aea  <- st_transform (M.gua, crs_aea)

# Vicugna vicugna
M.men_aea  <- st_transform (M.men, crs_aea)
M.vhyb_aea <- st_transform (M.vhyb, crs_aea)
M.vic_aea  <- st_transform (M.vic, crs_aea)

#Calculate geographic overlap -------------------------------------------------
#This function calculates:
# 1. Total accessible area (M) of the guanaco group.
# 2. Total accessible area (M) of the vicuña group.
# 3. Geographic area shared by both groups.
# 4. Percentage of the guanaco M overlapping the vicuña M.
# 5. Percentage of the vicuña M overlapping the guanaco M.

geo_overlap <- function(
    M_g, M_v) {
  # Dissolve all polygon features within each accessible area.
  M_g <- st_union(M_g)
  M_v <- st_union(M_v)
  
  #Calculate total accessible areas in km2.
  area_g <- as.numeric(st_area(M_g)) / 1e6
  area_v <- as.numeric(st_area(M_v)) / 1e6
  
  #Calculate geographic intersection.
  inter <- suppressWarnings(st_intersection(M_g, M_v))
  
  #Calculate shared area in km2.
  if (
    length(inter) == 0 ||
    all(st_is_empty(inter))) {
    area_shared <- 0
  
  } else {
    
    area_shared <- sum(
      as.numeric(
        st_area(inter))) / 1e6}

  #Calculate overlap relative to the guanaco accessible area.
  overlap_g_pct <- 100 * area_shared / area_g
  
  #Calculate overlap relative to the vicuña accessible area.
  overlap_v_pct <-100 * area_shared / area_v
  
  data.frame(
    guanaco_M_km2 = area_g,
    vicuna_M_km2 = area_v,
    shared_M_km2 = area_shared,
    overlap_from_guanaco_pct = overlap_g_pct,
    overlap_from_vicuna_pct = overlap_v_pct)}


#Define pairwise comparisons ---------------------------------------------------
#All nine combinations between guanaco and vicuña groups are evaluated.
comparisons <- data.frame(
  
  guanaco = c(
    "cacsilensis",
    "cacsilensis",
    "cacsilensis",
    "guanaco_contact",
    "guanaco_contact",
    "guanaco_contact",
    "guanicoe",
    "guanicoe",
    "guanicoe"),
  
  vicuna = c(
    "mensalis",
    "vicuna_contact",
    "vicugna",
    "mensalis",
    "vicuna_contact",
    "vicugna",
    "mensalis",
    "vicuna_contact",
    "vicugna"))

#Store projected M polygons in named lists ------------------------------------
guanaco.M <- list(
  cacsilensis = M.cac_aea,
  guanaco_contact = M.lhyb_aea,
  guanicoe = M.gua_aea)

vicuna.M <- list(
  mensalis = M.men_aea,
  vicuna_contact = M.vhyb_aea,
  vicugna = M.vic_aea)

#Calculate all pairwise geographic overlaps -----------------------------------
results_geo <- lapply(
  seq_len(
    nrow(comparisons)),
  function(i) {
    result <- geo_overlap (
      M_g = guanaco.M[
        [comparisons$guanaco[i]]],
      M_v = vicuna.M[
        [comparisons$vicuna[i]]] )
    
    data.frame(
      guanaco = comparisons$guanaco[i],
      vicuna = comparisons$vicuna[i],
      result
    )}) %>%
  bind_rows()

# Inspect results ---------------------------------------------------------------
results_geo

#Create rounded summary table --------------------------------------------------
results_geo_summary <- results_geo %>%
  select(
    guanaco,
    vicuna,
    shared_M_km2,
    overlap_from_guanaco_pct,
    overlap_from_vicuna_pct
  ) %>%
  mutate(
    shared_M_km2 = round(
      shared_M_km2,
      1),
    overlap_from_guanaco_pct = round(
      overlap_from_guanaco_pct,
      1),
    overlap_from_vicuna_pct = round(
      overlap_from_vicuna_pct,
      1))

results_geo_summary

# Export geographic overlap results --------------------------------------------
write.csv(results_geo, "output/geographic_overlap/geographic_overlap_M.csv",
  row.names = FALSE)