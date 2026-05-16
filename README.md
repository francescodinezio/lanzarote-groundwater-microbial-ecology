# Lanzarote groundwater microbial ecology

R scripts and supporting workflows used for the ecological, statistical, and visualization analyses presented in:

**Di Nezio et al. — Microbial diversity of groundwater-related environments in Lanzarote**

## Overview

This repository contains the code used to process and analyze bacterial 16S rRNA gene amplicon data from groundwater-dependent ecosystems of Lanzarote (Canary Islands, Spain). The analyses focus on:

* Whole-community bacterial diversity and composition
* Functional subsets of the microbiome
* Potentially pathogenic taxa
* Predicted human-associated taxa
* Environment-associated taxa
* Beta-diversity partitioning
* Indicator taxa analyses
* Random forest classification
* Figure generation for the manuscript and supplementary material

The repository includes scripts for downstream ecological analyses and visualization of ASV-based datasets generated after DADA2 processing.

## Data description

The analyses are based on ASV tables generated from 16S rRNA gene amplicon sequencing targeting the V3–V4 region using primers 341F/805R.

Habitats included in the study:

* Anchialine caves
* Anchialine pools
* Enclosed marine bays
* Saltwork pans
* Wells and water galleries
* Spring-fed ponds

## Main analyses

The scripts include:

* Community matrix cleaning and filtering
* Taxonomic aggregation
* Alpha diversity analyses
* NMDS ordination
* PERMANOVA analyses
* Beta-diversity decomposition
* Indicator species analyses
* Functional group assignment
* Heatmap generation
* Random forest classification analyses
* Figure assembly and export

## Requirements

Main R packages used include:

* tidyverse
* vegan
* indicspecies
* BAT
* patchwork
* glmmTMB
* emmeans
* performance
* randomForest
* caret
* ape
* phangorn

Additional custom helper scripts are stored in the `functions/` directory.

## Reproducibility

The repository is intended to provide full reproducibility of the ecological and statistical analyses presented in the manuscript. Raw sequencing data are available in the corresponding public sequence repository described in the manuscript.

## Citation

If using this repository, please cite:

Di Nezio F., Martinez A., et al. *Microbial diversity of groundwater-related environments in Lanzarote*.

## License

This repository is released under the MIT License.

