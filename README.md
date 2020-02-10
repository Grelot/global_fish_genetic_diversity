Codes for the paper : "Global determinants of freshwater and marine fish genetic diversity"
================================================


[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/3545)

#### Stephanie Manel, Pierre-Edouard Guerin, David Mouillot, Simon Blanchet, Laure Velez, Camille Albouy, Loic Pellissier


Montpellier, 2017-2019  

Published in Nature Communications, 2020  
full-text acces: https://rdcu.be/b1sXy

---------------------------------------------------------------------

<img src="https://github.com/Grelot/global_fish_genetic_diversity/blob/master/99-utils/logo_shiny.png" href="https://shiny.cefe.cnrs.fr/wfgd/" width="32" height="32"> A web application is available to display Figure 1 with more details: https://shiny.cefe.cnrs.fr/wfgd/  


<img src="https://github.com/Grelot/global_fish_genetic_diversity/blob/master/99-utils/gitlab_logo.png" href="https://gitlab.mbb.univ-montp2.fr/reservebenefit/worldmap_fish_genetic_diversity" width="32" height="32"> Codes also availables on gitlab: https://gitlab.mbb.univ-montp2.fr/reservebenefit/worldmap_fish_genetic_diversity   



---------------------------------------------------------------------


# Table of contents

1. [Introduction](#1-introduction)
2. [Installation](#2-installation)
    1. [Prerequisites](#21-prerequisites)
    2. [Singularity container](#22-singularity-container)
    3. [Data Files](#23-data-files)
    4. [Set up](#24-set-up)
3. [Scripts Code Source](#3-scripts-code-source)  
4. [Running the pipeline](#4-running-the-pipeline)
    1. [Filter raw data](#41-filter-raw-data)
    2. [Georeferenced sequences alignments by species](#42-data-files)
    3. [Species sequence pairwise comparison](#43-species-sequence-pairwise-comparison)
    4. [Genetic diversity calculation](#44-genetic-diversity-calculation)
    5. [Statistical analysis](#45-statistical-analysis)
        1. [Merge genetic data with environmental data by cell](#451-merge-genetic-data-with-environmental-data-by-cell)        
        2. [Figures](#452-figures)
        3. [Analysis](#453-analysis)
        4. [Supplementary figures](#454-supplementary-figures)
    6. [Taxonomy and habitat attributed to each individual sequences](#46-taxonomy-and-habitat-attributed-to-each-individual-sequences)


# 1. Introduction

This repository contains all the scripts to reproduce the results of the paper Manel et al. (2019) from the georeferenced barcode sequences of the supergroup "actinopterygii" downloaded from [BOLD](http://www.boldsystems.org/index.php/Public_SearchTerms?taxon=&searchMenu=records&query=actinopterygii) on 17th september 2018. 

The pipeline is composed of 6 steps : 

1. [Filter raw data](#41-filter-raw-data)
2. [Georeferenced sequences alignments by species](#42-data-files)
3. [Species sequence pairwise comparison](#43-species-sequence-pairwise-comparison)
4. [Genetic Diversity calculation](#44-genetic-diversity-calculation)
5. [Statistical analysis](#45-statistical-analysis)
6. [Taxonomy and habitat attributed to each individual sequences](#46-taxonomy-and-habitat-attributed-to-each-individual-sequences)

Figures and statistical analysis can be reproduced directly (see [Figures](#452-figures) section) without running the whole pipeline.

Only datafiles necessary to initiate the whole pipeline as well as to produce figures and statisticial analysis are provided. 

# 2. Installation

## 2.1 Prerequisites
You must install the following softwares and packages to run all steps:
For [Figures and statiscal analysis](#45-statistical-analysis), only R packages are needed. 

- [JULIA Version 1.1.0](https://julialang.org/)
    - `julia-module` DataFrames
    - `julia-module` DelimitedFiles
    - `julia-module` DataFramesMeta
    - `julia-module` StatsBase
    - `julia-module` Statistics
    - `julia-module` CSV
- [R Version 3.2.3](https://cran.r-project.org/)
    - `R-package` raster
    - `R-package` plotrix
    - `R-package` sp
    - `R-package` maptools
    - `R-package` parallel
    - `R-package` png
    - `R-package` plyr
    - `R-package` shape
    - `R-package` MASS
    - `R-package` hier.part
    - `R-package` countrycode
    - `R-package` sjPlot
    - `R-package` gridExtra
    - `R-package` ggplot2
    - `R-package` lme4
    - `R-package` SpatialPack
    - `R-package` rgeos | if install.packages("rgeos") failed, then try: install.packages("https://cran.r-project.org/src/contrib/Archive/rgeos/rgeos_0.3-26.tar.gz", type="source")
    - `R-package` rgdal | it may require to install "libgdal-dev"
    - `R-package` [rfishbase](https://cran.r-project.org/web/packages/rfishbase/rfishbase.pdf)
    - `R-package` pgirmess
    - `R-package` car
- [Python Version 3.6.8](https://www.python.org/downloads/release/python-368/)
    - `python3-module` argparse
    - `python3-module` re
    - `python3-module` [ete3](http://etetoolkit.org/)
    - `python3-module` numpy
    - `python3-module` csv
    - `python3-module` re
    - `python3-module` csv
    - `python3-module` difflib
- [MUSCLE Version 3.8.31](https://www.drive5.com/muscle/)

## 2.2 Singularity container

Alternatively, you can download and use a singularity container with all prerequisites (R, Julia, Python, Muscle).

### Install Singularity

See https://www.sylabs.io/docs/ for instructions to install Singularity.

### Download the container

```
singularity pull --name global_fish_genetic_diversity.simg shub://Grelot/global_fish_genetic_diversity:global_fish_genetic_diversity
```

### Use the container

This command will spawn a shell environment with all prerequisites.
```
singularity shell global_fish_genetic_diversity.simg
```



## 2.3 Data files
The included data files are :

* [seqbold_data.tsv](02-raw_data/seqbold_data.tsv) : georeferenced barcode sequences of the supergroup "actinopterygii" downloaded from [BOLD](http://www.boldsystems.org/index.php/Public_SearchTerms?taxon=&searchMenu=records&query=actinopterygii) on 17th september 2018
* [grid_equalarea200km](01-infos/grid_equalarea200km) : shapefile of worldmap equal area projection epsg:4326 with nested equal area grids (cell sizes of 200km)
* [equalarea_id_coords.tsv](01-infos/equalarea_id_coords.tsv) : ID and left/right/top/bottom coordinates of each equal area into the shapefile grid_equalarea200km.
* [marine_actinopterygii_species.txt](01-infos/marine_actinopterygii_species.txt) : list of "actinopterygii" saltwater species according to [fishbase](http://www.fishbase.org/)
* [marine_bo_o2dis.asc](01-infos/spatial_layers/marine_bo_o2dis.asc) : spatial layer of marine oxygen concentration [mol/l] from [gmed](http://gmed.auckland.ac.nz/)
* [marine_bo_sst_mean.asc](01-infos/spatial_layers/marine_bo_sst_mean.asc) : spatial layer of sea surface temperature from [gmed](http://gmed.auckland.ac.nz/)
* [marine_velocity_mean.asc](01-infos/spatial_layers/marine_velocity_mean.asc) : spatial layer of velocity of velocity (marine) from [gmed](http://gmed.auckland.ac.nz/)
* [freshwater_wc2.0_bio_10m_01.tif](01-infos/spatial_layers/freshwater_wc2.0_bio_10m_01.tif) : spatial layer of global mean temperature from [worldclim](http://worldclim.org/version2)
* [freshwater_velocity_mean.tif](01-infos/spatial_layers/freshwater_velocity_mean.tif) : spatial layer of velocity (freshwater) from [gmed](http://gmed.auckland.ac.nz/)
* [datatoFigshare](01-infos/datatoFigshare) : shapefile of drainage basins from [a global database on freshwater fish species occurrence in drainage basins](https://www.nature.com/articles/sdata2017141)
* [datacell_grid_descriteurs.csv](01-infos/datacell_grid_descriteurs.csv) : table of bathymetry, chlorophyll and other information attributed to each cell extracted from [gebco](https://www.gebco.net/) and [bio-oracle](http://www.bio-oracle.org/)
* [equalarea_id_coordsCA_FWRS_MR_RS.csv](01-infos/equalarea_id_coordsCA_FWRS_MR_RS.csv) : Table of species diversity by 200km square gridcell from [OBIS](http://www.iobis.org)
* [ne_50m_land](01-infos/ne_50m_land) : shapefile of worldcoast from [naturalearthdata](http://www.naturalearthdata.com)
* [ne_50m_rivers_lake_centerlines_scale_rank](01-infos/ne_50m_rivers_lake_centerlines_scale_rank) : shapefile of riverlines from [naturalearthdata](http://www.naturalearthdata.com)
* [GSHHS_h_L2.shp](01-infos/big_lakes/GSHHS_h_L2.shp) : shape polygon file of big lakes from [naturalearthdata](http://www.naturalearthdata.com)
* [EnvFreshwater.csv](01-infos/EnvFreshwater.csv) : slope and flow information for each geographical cell with a river from [gebco](https://www.gebco.net/)
* [distanceCote](01-infos/distanceCote.csv) : distance from shore for each cell calculated from [gmed](http://gmed.auckland.ac.nz/)
* _ornament_ [fishes images](01-infos/images) :free silhouette images of fishes from [phylopic](http://phylopic.org/)

## 2.4 Set Up
Clone the project and switch to the main folder, it's your working directory

```
git clone http://gitlab.mbb.univ-montp2.fr/reservebenefit/worldmap_fish_genetic_diversity.git
cd worldmap_fish_genetic_diversity
```

You're ready to run the analysis. Now follow the instructions at [Running the pipeline](#4-running-the-pipeline)

# 3. Scripts code source

For more information, we provide a detailed description of each scripts in [SCRIPTS.md](SCRIPTS.md)

# 4. Running the pipeline

## 4.1 Filter raw data
1. Keep only the CO1 sequences with lat/lon information
    * input : [seqbold_data.tsv](02-raw_data/seqbold_data.tsv)
    * output : [co1_ssll_seqbold_data.tsv](03-filtered_data/co1_ssll_seqbold_data.tsv)
```
bash 00-scripts/step1/filter_raw_data.sh
```

## 4.2 Georeferenced sequences alignments by species
1. Align sequences from the same species with MUSCLE and create coordinates _.coord_ file for each sequence
    * input : [co1_ssll_seqbold_data.tsv](03-filtered_data/co1_ssll_seqbold_data.tsv)
    * output : [05-species_alnt](05-species_alnt) .fasta &.coords files
```
bash 00-scripts/step2/seq_alnt_filtered_data.sh
```

2. According to a list of marine species, move fasta and coords files into _marine_ or _freshwater_ folder
    * input : [marine_actinopterygii_species.txt](01-infos/marine_actinopterygii_species.txt), [05-species_alnt](05-species_alnt)
    * output : [/06-species_alnt_cluster/marine](/06-species_alnt_cluster/marine), [/06-species_alnt_cluster/freshwater](/06-species_alnt_cluster/freshwater)
```
mkdir 06-species_alnt_cluster/total
mkdir 06-species_alnt_cluster/freshwater
mkdir 06-species_alnt_cluster/marine
bash 00-scripts/step2/cluster_freshwater_vs_marine.sh
```

3. Attribute at each individual sequences an ID of cell of the shapefile of worldmap equal area projection from its coordinates
    * input : [grid_equalarea200km](01-infos/grid_equalarea200km), [/06-species_alnt_cluster/marine](/06-species_alnt_cluster/marine), [/06-species_alnt_cluster/freshwater](/06-species_alnt_cluster/freshwater) .coords files
    * output : [/06-species_alnt_cluster/marine](/06-species_alnt_cluster/marine), [/06-species_alnt_cluster/freshwater](/06-species_alnt_cluster/freshwater) .equalareacoords files
```
Rscript 00-scripts/step2/equalareacoords.R
```

## 4.3 Species sequence pairwise comparison
1. Generate individual sequences pairwise comparison data matrices for each species for both each cell and each latitudinal band from species sequences alignments and cell locations.
    * input : [06-species_alnt_cluster](06-species_alnt_cluster)
    * output : [07-master_matrices](07-master_matrices) 
```
julia 00-scripts/step3/master_matrices.jl
```

## 4.4 Genetic diversity calculation
1. Attribute mean genetic diversity value at both each cell and each latitudinal band. In latitudinal band case, we filter out species with no genetic diversity or/and less than 3 individuals.Standard deviation are estimated from 1000 bootstrapped replications.
    * input : [07-master_matrices](07-master_matrices)
    * output : [equalarea_numbers.csv](08-genetic_diversity/)
```
julia 00-scripts/step4/equalarea_numbers.jl
julia 00-scripts/step4/latband_numbers.jl
```

2. Generate a table of cell_coordinates, mean of Genetic diversity, number of species, mean/sd number of individuals by species into each cell.
    * inputs : [equalarea_numbers.csv](08-genetic_diversity/), [gdval_by_area.csv](08-genetic_diversity/) 
    * outputs : [metrics_by_area_freshwater.csv](08-genetic_diversity/metrics_by_area_freshwater.csv), [metrics_by_area_marine.csv](08-genetic_diversity/metrics_by_area_marine.csv)
```
bash 00-scripts/step4/gdval_by_cell.sh
julia 00-scripts/step4/metrics_by_area_and_species.jl
```

## 4.5 Statistical analysis
To generate figures (or analysis), you can simply type the proposed commands on your terminal or alternatively open the R script and run the commands in R 

### 4.5.1 Merge genetic data with environmental data by cell
1. Generate a table of center of cell (xy) coordinates, ID of cell, mean of Genetic diversity, number of species, mean/sd number of individuals by species into each cell, bathymetry, chlorophyll concentration, oxygen concentration, temperature, drainage basin surface area into each cell
    * inputs : [metrics_by_area_freshwater.csv](08-genetic_diversity/metrics_by_area_freshwater.csv),[datatoFigshare](01-infos/datatoFigshare),[datacell_grid_descriteurs.csv](01-infos/datacell_grid_descriteurs.csv)
    * output : [total_data_genetic_diversity_with_all_descripteurs.tsv](09-all_descripteurs/total_data_genetic_diversity_with_all_descripteurs.tsv)
```
Rscript 00-scripts/step5/descripteurs.R
```

### 4.5.2 Figures
1.  Map of the global distribution of genetic diversity for marine species
    * inputs : [total_data_genetic_diversity_with_all_descripteurs.tsv](09-all_descripteurs/total_data_genetic_diversity_with_all_descripteurs.tsv), [grid_equalarea200km](01-infos/grid_equalarea200km), [ne_50m_rivers_lake_centerlines_scale_rank](01-infos/ne_50m_rivers_lake_centerlines_scale_rank), [GSHHS_h_L2.shp](01-infos/big_lakes/GSHHS_h_L2.shp)
    * output : [figure1.eps](10-figures/figure1.eps)
```
Rscript 00-scripts/step5/figures/figure1.R
```

2. Congruence between fish genetic and species diversity
    * inputs : [total_data_genetic_diversity_with_all_descripteurs.tsv](09-all_descripteurs/total_data_genetic_diversity_with_all_descripteurs.tsv), [equalarea_id_coordsCA_FWRS_MR_RS.csv](01-infos/equalarea_id_coordsCA_FWRS_MR_RS.csv), [grid_equalarea200km](01-infos/grid_equalarea200km), [ne_50m_rivers_lake_centerlines_scale_rank](01-infos/ne_50m_rivers_lake_centerlines_scale_rank), [GSHHS_h_L2.shp](01-infos/big_lakes/GSHHS_h_L2.shp)
    * output : [figure2.tiff](10-figures/figure2.tiff)
```
Rscript 00-scripts/step5/figures/figure2.R
```

3. Determinant of the patterns of fish genetic diversity
    * inputs : [total_data_genetic_diversity_with_all_descripteurs.tsv](09-all_descripteurs/total_data_genetic_diversity_with_all_descripteurs.tsv),[EnvFreshwater.csv](01-infos/EnvFreshwater.csv),[distanceCote](01-infos/distanceCote.csv),[equalarea_id_coordsCA_FWRS_MR_RS.csv](01-infos/equalarea_id_coordsCA_FWRS_MR_RS.csv)
    * output : [figure3.pdf](10-figures/figure3.pdf)
```
Rscript 00-scripts/step5/figures/figure3.R
```

### 4.5.3 Analysis

1. Wilcoxon test to assess whether genetic diversity means differ between marine and freshwater species.
```
Rscript 00-scripts/step5/analysis/wilcoxon_tests.R
```

2. Sensitivity analysis
    * inputs : [total_data_genetic_diversity_with_all_descripteurs.tsv](09-all_descripteurs/total_data_genetic_diversity_with_all_descripteurs.tsv), [equalarea_id_coordsCA_FWRS_MR_RS.csv](01-infos/equalarea_id_coordsCA_FWRS_MR_RS.csv), [grid_equalarea200km](01-infos/grid_equalarea200km)
    * outputs : [marine_sensitivity.csv](13-analysis/marine_sensitivity_model.csv), [freshwater_sensitivity_model.csv](13-analysis/freshwater_sensitivity_model.csv)
```
Rscript 00-scripts/step5/analysis/sensitive_analysis_model.R
```

3. Sensitivity analysis based on taxonomic coverage

    * inputs : [total_data_genetic_diversity_with_all_descripteurs.tsv](09-all_descripteurs/total_data_genetic_diversity_with_all_descripteurs.tsv), [equalarea_id_coordsCA_FWRS_MR_RS.csv](01-infos/equalarea_id_coordsCA_FWRS_MR_RS.csv), [grid_equalarea200km](01-infos/grid_equalarea200km)
    * outputs : [marine_sensitivity_covtax.csv](13-analysis/marine_sensitivity_covtax.csv), [freshwater_sensitivity_covtax.csv](13-analysis/freshwater_sensitivity_covtax.csv)
```
Rscript 00-scripts/step5/analysis/sensitive_analysis_covtax.R
```

### 4.5.3 Supplementary figures
1. Spatial autocorrelogramme based on the I-Moran coefficient
    * inputs : [total_data_genetic_diversity_with_all_descripteurs.tsv](09-all_descripteurs/total_data_genetic_diversity_with_all_descripteurs.tsv), [grid_equalarea200km](01-infos/grid_equalarea200km)
    * output : [figureS1.pdf](10-figures/figureS1.pdf)
```
Rscript 00-scripts/step5/supplementary_figures/figureS1.R
```
2. Global distribution of higher and lower percentiles of genetic diversity
    * inputs : [total_data_genetic_diversity_with_all_descripteurs.tsv](09-all_descripteurs/total_data_genetic_diversity_with_all_descripteurs.tsv), [grid_equalarea200km](01-infos/grid_equalarea200km), [ne_50m_rivers_lake_centerlines_scale_rank](01-infos/ne_50m_rivers_lake_centerlines_scale_rank), [GSHHS_h_L2.shp](01-infos/big_lakes/GSHHS_h_L2.shp)
    * output : [figureS2.tiff](10-figures/figureS2.tiff)
```
Rscript 00-scripts/step5/supplementary_figures/figureS2.R
```
3. Latitudinal distribution of species diversity
    * inputs : [total_data_genetic_diversity_with_all_descripteurs.tsv](09-all_descripteurs/total_data_genetic_diversity_with_all_descripteurs.tsv), [equalarea_id_coordsCA_FWRS_MR_RS.csv](01-infos/equalarea_id_coordsCA_FWRS_MR_RS.csv), [grid_equalarea200km](01-infos/grid_equalarea200km)
    * output : [figureS3.pdf](10-figures/figureS3.pdf)
```
Rscript 00-scripts/step5/supplementary_figures/figureS3.R
```
4. Regional effect on the global genetic diversity pattern
    * input : [total_data_genetic_diversity_with_all_descripteurs.tsv](09-all_descripteurs/total_data_genetic_diversity_with_all_descripteurs.tsv)
    * output : [figureS4.pdf](10-figures/figureS4.pdf)
```
Rscript 00-scripts/step5/supplementary_figures/figureS4.R
```

5. Sampling effect
    * inputs : [total_data_genetic_diversity_with_all_descripteurs.tsv](09-all_descripteurs/total_data_genetic_diversity_with_all_descripteurs.tsv), [grid_equalarea200km](01-infos/grid_equalarea200km), [ne_50m_rivers_lake_centerlines_scale_rank](01-infos/ne_50m_rivers_lake_centerlines_scale_rank), [GSHHS_h_L2.shp](01-infos/big_lakes/GSHHS_h_L2.shp)
    * output : [figureS5.tiff](10-figures/figureS5.tiff)
```
Rscript 00-scripts/step5/supplementary_figures/figureS5.R
```

6. Taxonomic coverage of the sequences used by the model
    * input : [watertype_all_modeles_effectives_family.csv](11-sequences_taxonomy_habitat/watertype_all_modeles_effectives_family.csv)
    * output : [figureS6.pdf](10-figures/figureS6.pdf)
```
Rscript 00-scripts/step5/supplementary_figures/figureS6.R
```

7. Spatial distribution of taxonomic coverage
    * inputs : [total_data_genetic_diversity_with_all_descripteurs.tsv](09-all_descripteurs/total_data_genetic_diversity_with_all_descripteurs.tsv), [grid_equalarea200km](01-infos/grid_equalarea200km), [ne_50m_rivers_lake_centerlines_scale_rank](01-infos/ne_50m_rivers_lake_centerlines_scale_rank), [GSHHS_h_L2.shp](01-infos/big_lakes/GSHHS_h_L2.shp), [equalarea_id_coordsCA_FWRS_MR_RS.csv](01-infos/equalarea_id_coordsCA_FWRS_MR_RS.csv)
    * ouputs : [species_taxon_coverage_number_cells.csv](12-taxonomic_coverage/species_taxon_coverage_number_cells.csv), [figureS7.tiff](10-figures/figureS7.tiff)
```
Rscript 00-scripts/step5/supplementary_figures/figureS7.R
```

8. Intraspecific genetic diversity mean in each 10° latitudinal bands (not in the paper)
    * inputs : [freshwater_latbands_bootstraps.csv](08-genetic_diversity/freshwater_latbands_bootstraps.csv), [marine_latbands_bootstraps.csv](08-genetic_diversity/marine_latbands_bootstraps.csv)
    * output : [figureS8.pdf](10-figures/figureS8.pdf)
```
Rscript 00-scripts/step5/supplementary_figures/figureS8.R
```


## 4.6 Taxonomy and habitat attributed to each individual sequences

1. Write a table of individual sequences with geographical cell localisation
    * inputs : [total_data_genetic_diversity_with_all_descripteurs.tsv](09-all_descripteurs/total_data_genetic_diversity_with_all_descripteurs.tsv), [06-species_alnt_cluster](06-species_alnt_cluster), [co1_ssll_seqbold_data.tsv](03-filtered_data/co1_ssll_seqbold_data.tsv)
    * output : [map_marine_sequences.csv](11-sequences_taxonomy_habitat/map_marine_sequences.csv)
```
python3 00-scripts/step6/sequences_table.py
```
2. Check if the watertype _marine_|_freshwater_ for each species by cell is correct according to the model _marine_|_freshwater_
    * inputs : [map_marine_sequences.csv](11-sequences_taxonomy_habitat/map_marine_sequences.csv), [metrics_by_area_marine.csv](08-genetic_diversity/metrics_by_area_marine.csv), [spatial_layers/](01-infos/spatial_layers/),
    * output : [wrong_freshwater_sequences.csv](11-sequences_taxonomy_habitat/wrong_freshwater_sequences.csv)
```
Rscript 00-scripts/step6/check_freshwater_assignation.R
```
3. Assign habitat (demersal, pelagic...) information to each individual sequences according to their attributed species name
    * input : [map_marine_sequences.csv](11-sequences_taxonomy_habitat/map_marine_sequences.csv)
    * output : [sequences_withdemerpelag.csv](11-sequences_taxonomy_habitat/sequences_withdemerpelag.csv)
```
Rscript 00-scripts/step6/sequences_demerpelag.R
```
4. Cure habitat assignation and species name for each individual sequences
    * input : [sequences_withdemerpelag.csv](11-sequences_taxonomy_habitat/sequences_withdemerpelag.csv)
    * output : [cured_sequences_withdemerpelag.csv](11-sequences_taxonomy_habitat/cured_sequences_withdemerpelag.csv)
```
Rscript 00-scripts/step6/sequences_cure_species_name.R
```
5. cure family column by renaming [BOLD](http://www.boldsystems.org/) family by its equivalent into [NCBI taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy)
    * input : [cured_sequences_withdemerpelag.csv](11-sequences_taxonomy_habitat/cured_sequences_withdemerpelag.csv)
    * output : [cured_family_sequences_withdemerpelag.csv](11-sequences_taxonomy_habitat/cured_family_sequences_withdemerpelag.csv)
```
bash 00-scripts/step6/rename_family_bold_to_ncbi.sh
```
6. Write a table of number of species/number of sequences by taxonomic order/family used for each model
    * inputs : [cured_family_sequences_withdemerpelag.csv](11-sequences_taxonomy_habitat/cured_family_sequences_withdemerpelag.csv)
    * output : [watertype_all_modeles_effectives_family.csv](11-sequences_taxonomy_habitat/watertype_all_modeles_effectives_family.csv)
```
python3 00-scripts/step6/sequences_taxonomy.py
```
