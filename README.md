# Where is the guanaco irreplaceable? Lineage-level niche structure and functional redundancy in wild South American camelids

**Andrea G. Castillo¹˒², Laura Jiménez³˒⁴, and Horacio Samaniego²˒⁵**

¹ Escuela de Graduados, Programa de Doctorado en Ciencias mención Ecología y Evolución, Facultad de Ciencias, Universidad Austral de Chile, Valdivia, Chile  
² Instituto de Conservación, Biodiversidad y Territorio, Laboratorio de Ecoinformática, Universidad Austral de Chile, Valdivia, Chile  
³ Centro de Investigación e Innovación para el Cambio Climático, Universidad Santo Tomás, Santiago, Chile  
⁴ Centro de Modelamiento Matemático, Universidad de Chile, Santiago, Chile  
⁵ Instituto de Sistemas Complejos, Valparaíso, Chile

## What this repository contains

This repository contains the data, processed environmental layers, and R scripts associated with the analysis of environmental niche structure in *Lama guanicoe* and *Vicugna vicugna* at the lineage and contact-zone population levels. The scripts reproduce the niche breadth, niche overlap, niche dynamics, geographic overlap, and niche similarity/divergence analyses described in the associated manuscript.

> **Note on ecological niche modelling:** The processed climatic layers used in the manuscript are provided in this repository; however, the workflow used for environmental-variable preparation and ellipsoidal ecological niche modelling is not reproduced here. These procedures are part of ongoing methodological work led by Laura Jiménez (LJ), including an R package currently in preparation for publication. The modelling component of this repository will be updated, where appropriate, once the associated methodological work becomes publicly available.
>

## Analytical framework

<p align="center">
  <img src="figures/Figure 1_analytical_framework.png" width="900">
</p>

<p align="center">
  <em>Conceptual and analytical framework used to evaluate environmental niche structure and redundancy among guanaco and vicuña lineages and contact-zone populations.</em>
</p>

## Repository layout

```text
.
├── data/
│   ├── environmental/            # Processed environmental variables
│   ├── occurrences/              # Occurrence records and environmental background points
│   └── spatial/                  # Accessible-area (M) polygons by group
│
├── scripts/
│   ├── 00_data_preparation.R     # Data cleaning and preparation
│   ├── 01_environmental_matrix.R # Extraction and organization of environmental data
│   ├── 02_PCA_analysis.R         # Common environmental PCA
│   ├── 03_niche_breadth.R        # Niche breadth and rarefaction analyses
│   ├── 04_niche_overlap.R        # Pairwise niche overlap (ecospat)
│   ├── 05_niche_dynamics.R       # Expansion, stability, and unfilling (ecospat)
│   ├── 06_geographic_M_overlap.R # Geographic overlap among accessible areas (M)
│   └── 07_NOT_NDT.R              # Niche Overlap and Niche Divergence Tests (humboldt)
│
├── figures/                      # Main and supplementary figures
├── supplementary/                # Supplementary tables and analysis outputs
└── README.md                     # Repository documentation
```

## Analyses

The analytical workflow implemented in this repository includes:

- **Environmental PCA** (`02_PCA_analysis.R`) — Characterization of the common environmental space used for lineage- and population-level niche analyses.

- **Niche breadth** (`03_niche_breadth.R`) — Estimation of environmental niche breadth using convex hull area, including rarefaction to account for differences in sample size among groups.

- **Niche overlap** (`04_niche_overlap.R`) — Pairwise estimation of niche overlap among lineages and contact-zone populations using Schoener's D and Warren's I in the `ecospat` framework.

- **Niche dynamics** (`05_niche_dynamics.R`) — Directional assessment of niche expansion, stability, and unfilling between guanaco and vicuña groups using `ecospat`.

- **Geographic overlap** (`06_geographic_M_overlap.R`) — Quantification of the spatial overlap between the accessible areas (M) of guanaco and vicuña groups.

- **Niche Overlap and Niche Divergence Tests** (`07_NOT_NDT.R`) — Niche equivalency and directional background tests implemented with the `humboldt` package, accounting for environmental availability and analogous environmental conditions.

## Data

The `data/` directory contains the processed datasets used by the analytical scripts and is organized into environmental, occurrence/background, and spatial data.

### Environmental data

The `data/environmental/` directory contains the three environmental variables used throughout the analyses:

- `bio6.asc` — Minimum temperature of the coldest month (BIO6)
- `bio14.asc` — Precipitation of the driest month (BIO14)
- `cmi_max.asc` — Maximum Climatic Moisture Index (CMImax)

The climatic variables were derived from [CHELSA v2.1](https://www.chelsa-climate.org/datasets) and processed to the spatial extent and resolution used in the study. Further details on environmental-variable preparation are provided in the associated manuscript and in the modelling note above.

### Occurrence and background data

The `data/occurrences/` directory contains the occurrence records assigned to each lineage or contact-zone population, together with the environmental background points sampled from their corresponding accessible areas (M).

Occurrence records for *Lama guanicoe* and *Vicugna vicugna* were compiled from the Global Biodiversity Information Facility (GBIF), iNaturalist, EutherianCoP (Mondanaro et al., 2025), and the South American Archaeological Isotopic Database (SAAID; Pezo-Lanfranco et al., 2024). For EutherianCoP and SAAID, only modern occurrence records of guanaco and vicuña were incorporated. These data were complemented with records from regional studies spanning the distributions of both species (Marín et al., 2008; Núñez, 2008; Castillo et al., 2018; Mesas et al., 2023; Rojas-Bonzi et al., 2024).

Files containing `occurrences` correspond to observed occurrence records, whereas files containing `Mpoints` represent environmental background points randomly sampled from the corresponding accessible area (M). Background samples comprise 10,000 points for each guanaco group and 5,000 points for each vicuña group.

### Accessible areas (M)

The `data/spatial/` directory contains the spatial files defining the accessible area (M) of each lineage and contact-zone population:

- *L. g. cacsilensis*
- *L. g. guanicoe*
- Guanaco contact-zone populations
- *V. v. mensalis*
- *V. v. vicugna*
- Vicuña contact-zone populations

Each group is stored in a separate directory containing the ESRI Shapefile and its associated component files.

Accessible areas were constructed by integrating spatial information from three sources:

1. current species distribution ranges reported by the [IUCN Red List of Threatened Species](https://www.iucnredlist.org/);
2. historical distribution ranges available from [PHYLACINE 1.2](https://megapast2future.github.io/PHYLACINE_1.2/) (Faurby et al., 2018); and
3. migration routes available through the [Atlas of Ungulate Migration](https://www.cms.int/gium), developed by the Global Initiative on Ungulate Migration (GIUM) under the Convention on the Conservation of Migratory Species of Wild Animals (CMS).

These sources were integrated to construct species-level accessible-area hypotheses. Lineage- and contact-zone-specific M hypotheses were subsequently delimited according to the geographic criteria described in the associated manuscript.
## Requirements

Analyses were conducted in R. The main packages used throughout the workflow include:

- `terra` and `sf` — spatial data processing
- `ade4` — principal component analysis
- `geometry` — convex hull estimation
- `ecospat` — niche overlap and niche dynamics
- `humboldt` — niche overlap and divergence tests
- `ggplot2` — data visualization

Package-specific requirements and additional dependencies are documented within the corresponding scripts.

## Citation

This repository accompanies the manuscript:

> Castillo, A.G., Jiménez, L., & Samaniego, H. *Where is the guanaco irreplaceable? Lineage-level niche structure and functional redundancy in wild South American camelids.* Manuscript in preparation.

Citation information will be updated upon publication.

## License and data reuse

The R scripts developed for this study may be reused with appropriate attribution under the license specified for this repository.

Occurrence, environmental, and spatial data included or referenced in this repository originate partly from third-party sources and remain subject to the terms, licenses, and citation requirements of their original providers. Users should consult the original data sources described in the Data section before redistributing or reusing these materials.

Please cite the associated manuscript and the original data providers when using materials from this repository.
