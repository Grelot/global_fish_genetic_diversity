# Detailed description of scripts

summary, inputs, outputs of scritps

Scripts are classified by step (see [Table of contents](#-table-of-contents)) and programmatic language(BASH  , PYTHON, JULIA, R).

# Table of contents

1. [Filter raw data](#1-filter-raw-data)
2. [Georeferenced sequences alignments by species](#2-data-files)
3. [Species sequence pairwise comparison](#3-species-sequence-pairwise-comparison)
4. [Genetic Diversity calculation](#4-genetic-diversity-calculation)
5. [Statistical analysis](#5-statistical-analysis)
6. [Taxonomy and habitat attributed to each individual sequences](#6-taxonomy-and-habitat-attributed-to-each-individual-sequences)

# 1 Filter raw data
#### BASH scripts
* [filter_raw_data.sh](00-scripts/step1/filter_raw_data.sh) : it keeps only the CO1 sequences with lat/lon information.
    - input :
        * [seqbold_data.tsv](02-raw_data/seqbold_data.tsv) : georeferenced barcode sequences from the supergroup "actinopterygii"  from [BOLD](http://www.boldsystems.org/)
    - output :
        * [co1_ssll_seqbold_data.tsv](03-filtered_data/co1_ssll_seqbold_data.tsv) : table of fitlered CO1 sequences with lat/lon and [BOLD](http://www.boldsystems.org/)'s taxonomy information
* [get_geonames_coordinates.sh](00-scripts/step1/get_geonames_coordinates.sh) : it uses [geonames.org](http://www.geonames.org/) to attribute coordinates lat/lon of individual sequences from their textual information of location when lat/lon is missing.

#### PYTHON2 scripts
* [lat_long_DMS_DD_converter.py](00-scripts/step1/lat_long_DMS_DD_converter.py) : it converts from DMS format to DD format the given coordinates.

# 2 Georeferenced sequences alignments by species
#### BASH scripts
* [seq_alnt_filtered_data.sh](00-scripts/step2/seq_alnt_filtered_data.sh) : aligns sequences from the same species with MUSCLE and creates coordinates _.coord_ file for each sequence.
    - inputs :
        * [co1_ssll_seqbold_data.tsv](03-filtered_data/co1_ssll_seqbold_data.tsv) : table of fitlered CO1 sequences with lat/lon
    - outputs :
        * [{species}.fasta](05-species_alnt) : alignment files of each {species} species
        * [{species}.coords](05-species_alnt) : coordinates lat/lon of each individual sequences of each {species} species
* [cluster_freshwater_vs_marine.sh](00-scripts/step2/cluster_freshwater_vs_marine.sh) : according to a list of marine species, moves the fasta and coords files into marine, freshwater repertories.
    - inputs : 
        * [05-species_alnt](05-species_alnt) : species fasta and coordinates files
        * [marine_actinopterygii_species.txt](01-infos/marine_actinopterygii_species.txt) : list of "actinopterygii" saltwater species according to [fishbase](http://www.fishbase.org/)
    - outputs :
        * [06-species_alnt_cluster/freshwater](/06-species_alnt_cluster/freshwater) : _freshwater_ species fasta & coordinates files
        * [06-species_alnt_cluster/marine](/06-species_alnt_cluster/marine) : _marine_ species fasta & coordinates files

#### PYTHON2 scripts
* [fasta_coords_files_species_generator.py](00-scripts/step2/fasta_coords_files_species_generator.py) : extracts sequences and associated coordinates from the filtered data.
    - input :
        * [co1_ssll_seqbold_data.tsv](03-filtered_data/co1_ssll_seqbold_data.tsv) : table of fitlered CO1 sequences with lat/lon
    - outputs :
        * [{species}.fasta](04-species_seq) : alignment files of each {species} species
        * [{species}.coords](04-species_seq) : coordinates lat/lon of each individual sequences of each {species} species

#### R scripts
* [equalareacoords.R](00-scripts/step2/equalareacoords.R) : attributes at each individual sequences an ID of cell of the shapefile of worldmap equal area projection from its coordinates.
    - inputs :
        * [06-species_alnt_cluster/freshwater/{species}.coords](/06-species_alnt_cluster/freshwater) : coordinates files by _freshwater_ {species} species
        * [06-species_alnt_cluster/marine/{species}.coords](/06-species_alnt_cluster/marine) : coordinates files by _marine_ {species} species
        * [grid_equalarea200km](01-infos/grid_equalarea200km) : shapefile of worldmap equal area projection epsg:4326 with nested equal area grids (cell sizes of 200km)
    - outputs :
        * [06-species_alnt_cluster/freshwater/{species}.equalareacoords](/06-species_alnt_cluster/freshwater) : cells of _freshwater_ {species} species
        * [06-species_alnt_cluster/marine/{species}.equalareacoords](/06-species_alnt_cluster/marine) : cells of _marine_ {species} species

# 3 Species sequence pairwise comparison
#### JULIA scripts
* [Lib_Compare_Pairwise.jl](00-scripts/step3/Lib_Compare_Pairwise.jl) : functions to compute the Genetic Diversity value from a set of sequences.
* [Lib_Create_Master_Matrices.jl](00-scripts/step3/Lib_Create_Master_Matrices.jl) : functions to create master data matrices that are used to compute genetic diversity.
* [master_matrices.jl](00-scripts/step3/master_matrices.jl) : generates master data matrices from species sequences alignments.
    - input :
        * [06-species_alnt_cluster](06-species_alnt_cluster) : .fasta, .equalareacoords files for each species
    - output :
        * [07-master_matrices](07-master_matrices) : individual sequences pairwise comparison data matrices for each species for each cell 

# 4 Genetic Diversity calculation
#### BASH scripts
* [gdval_by_cell.sh](00-scripts/step4/gdval_by_cell.sh) : generates CSV files with 2 columns : cell ID and mean genetic diversity per species into the cell.

#### JULIA scripts   
* [Lib_GD_summary_functions.jl](00-scripts/step4/Lib_GD_summary_functions.jl) : functions to calculate genetic diversity at species level and cell level
* [equalarea_numbers.jl](00-scripts/step4/equalarea_numbers.jl) : attributes mean genetic diversity at each equal area grid cell. Genetic diversity is calculated from master data matrices.
    - input :
        * [07-master_matrices](07-master_matrices) : individual sequences pairwise comparison data matrices for each species for each cell 
    - output :
        * [equalarea_numbers.csv](08-genetic_diversity/) : genetic diversity by cell
* [metrics_by_area_and_species.jl](00-scripts/step4/metrics_by_area_and_species.jl): it generates files for statistical analysis at next step : mean genetic diversity per cell, genetic diversity per species per cell, number of individuals per species, number of species per cell, cell coordinates, cell ID...
    - input : 
        * [equalarea_numbers.csv](08-genetic_diversity/) : genetic diversity by cell
        * [gdval_by_area.csv](08-genetic_diversity/) : CSV files with 2 columns : cell ID and mean genetic diversity per species into the cell
    - outputs :
        * [metrics_by_area_marine.csv](08-genetic_diversity/metrics_by_area_marine.csv) : table of ID_cell,ISO3,is_sea,cloMeanVal,cloMinVal,cloMaxVal,bathyVal,AP,HDI_2015,fshD information attributed to each cell with _marine_ species
        * [metrics_by_area_freshwater.csv](08-genetic_diversity/metrics_by_area_freshwater.csv) : table of ID_cell,ISO3,is_sea,cloMeanVal,cloMinVal,cloMaxVal,bathyVal,AP,HDI_2015,fshD information attributed to each cell with _freshwater_ species

* [latband_numbers.jl](00-scripts/step4/latband_numbers.jl) : attributes mean genetic diversity at each latitudinal band.
    - input : 
        * [pairwise_latbands.csv](07-master_matrices/) : genetic diversity per species by 10° latitudinal band
    - outputs :
        * [08-genetic_diversity/marine_latbands_numbers.csv](08-genetic_diversity/marine_latbands_numbers.csv) : table of ID_latband,ISO3,is_sea,cloMeanVal,cloMinVal,cloMaxVal,bathyVal,AP,HDI_2015,fshD information attributed to each latitudinal latband with _marine_ species 
        * [08-genetic_diversity/freshwater_latbands_numbers.csv](08-genetic_diversity/freshwater_latbands_numbers.csv) : table of ID_latband,ISO3,is_sea,cloMeanVal,cloMinVal,cloMaxVal,bathyVal,AP,HDI_2015,fshD information attributed to each latitudinal latband with _freshwater_ species
        * [marine_latbands_bootstraps.csv](08-genetic_diversity/marine_latbands_bootstraps.csv) : table of ID_latband,ISO3,is_sea,cloMeanVal,cloMinVal,cloMaxVal,bathyVal,AP,HDI_2015,fshD information attributed to each latitudinal latband with _marine_ species 
        * [freshwater_latbands_bootstraps.csv](08-genetic_diversity/freshwater_latbands_bootstraps.csv) : table of ID_latband,ISO3,is_sea,cloMeanVal,cloMinVal,cloMaxVal,bathyVal,AP,HDI_2015,fshD information attributed to each latitudinal latband with _freshwater_ species
  
* [Lib_bootstrap.jl](00-scripts/step4/Lib_bootstrap.jl) : functions to bootstrap by species latitudinal band genetic diversity


# 5 Statistical analysis
#### R scripts
* [functions.R](00-scripts/step5/figures/functions.R) : library of R functions required to run figures scripts.
* [descripteurs.R](00-scripts/step5/descripteurs.R) : from genetic data and geographic,environmental data, this script generates a table with cell as row and genetic,environmental,geographic variables as column.
    - inputs :
        * [metrics_by_area_freshwater.csv](08-genetic_diversity/metrics_by_area_freshwater.csv) : table of ID_cell,ISO3,is_sea,cloMeanVal,cloMinVal,cloMaxVal,bathyVal,AP,HDI_2015,fshD information attributed to each cell with _freshwater_ species
        * [metrics_by_area_marine.csv](08-genetic_diversity/metrics_by_area_marine.csv) : table of ID_cell,ISO3,is_sea,cloMeanVal,cloMinVal,cloMaxVal,bathyVal,AP,HDI_2015,fshD information attributed to each cell with _marine_ species
        * [marine_bo_o2dis.asc](01-infos/spatial_layers/marine_bo_o2dis.asc) : spatial layer of marine oxygen concentration [mol/l]
        * [marine_bo_sst_mean.asc](01-infos/spatial_layers/marine_bo_sst_mean.asc) : spatial layer of sea surface temperature
        * [marine_velocity_mean.asc](01-infos/spatial_layers/marine_velocity_mean.asc) : spatial layer of velocity of velocity (marine)
        * [freshwater_wc2.0_bio_10m_01.tif](01-infos/spatial_layers/freshwater_wc2.0_bio_10m_01.tif) : spatial layer of global mean temperature
        * [freshwater_velocity_mean.tif](01-infos/spatial_layers/freshwater_velocity_mean.tif) : spatial layer of velocity (freshwater)
        * [datatoFigshare](01-infos/datatoFigshare) : shapefile of drainage basins
        * [datacell_grid_descriteurs.csv](01-infos/datacell_grid_descriteurs.csv) : table of ID_cell,ISO3,is_sea,cloMeanVal,cloMinVal,cloMaxVal,bathyVal,AP,HDI_2015,fshD information attributed to each cell
    - output :
        * [total_data_genetic_diversity_with_all_descripteurs.tsv](09-all_descripteurs/total_data_genetic_diversity_with_all_descripteurs.tsv) : table of center of cell (xy) coordinates, ID of cell, mean of Genetic diversity, number of species, mean/sd number of individuals by species into each cell, bathymetry, chlorophyll concentration, oxygen concentration, temperature, drainage basin surface area into each cell
* [figure1.R](00-scripts/step5/figures/figure1.R) : from table of genetic,environmental,geographic variables by cell and shapefiles, it generates maps of the global distribution of genetic diversity as a _tiff_ file.
    - inputs :
        * [total_data_genetic_diversity_with_all_descripteurs.tsv](09-all_descripteurs/total_data_genetic_diversity_with_all_descripteurs.tsv) : table of center of cell (xy) coordinates, ID of cell, mean of Genetic diversity, number of species, mean/sd number of individuals by species into each cell, bathymetry, chlorophyll concentration, oxygen concentration, temperature, drainage basin surface area into each cell
        * [grid_equalarea200km](01-infos/grid_equalarea200km) : shapefile of worldmap equal area projection epsg:4326 with nested equal area grids (cell sizes of 200km)
        * [ne_50m_land](01-infos/ne_50m_land) : shapefile of worldcoast from (http://www.naturalearthdata.com)
        * [ne_50m_rivers_lake_centerlines_scale_rank](01-infos/ne_50m_rivers_lake_centerlines_scale_rank) : shapefile of riverlines from (http://www.naturalearthdata.com)
        * [GSHHS_h_L2.shp](01-infos/big_lakes/GSHHS_h_L2.shp) : shape polygon file of big lakes from (http://www.naturalearthdata.com)
    - output : 
        * [figure1.tiff](10-figures/figure1.tiff)
* [figure2.R](00-scripts/step5/figures/figure2.R) : from table of genetic,environmental,geographic variables by cell and species diversity , this script generates figures of the congruence between fish genetic and species diversity.
    - inputs :
        * [total_data_genetic_diversity_with_all_descripteurs.tsv](09-all_descripteurs/total_data_genetic_diversity_with_all_descripteurs.tsv)
        * [grid_equalarea200km](01-infos/grid_equalarea200km) : shapefile of worldmap equal area projection epsg:4326 with nested equal area grids (cell sizes of 200km)
        * [ne_50m_land](01-infos/ne_50m_land) : shapefile of worldcoast 
        * [ne_50m_rivers_lake_centerlines_scale_rank](01-infos/ne_50m_rivers_lake_centerlines_scale_rank) : shapefile of riverlines
        * [GSHHS_h_L2.shp](01-infos/big_lakes/GSHHS_h_L2.shp) : shape polygon file of big lakes
        * [equalarea_id_coordsCA_FWRS_MR_RS.csv](01-infos/equalarea_id_coordsCA_FWRS_MR_RS.csv) : Table of species diversity by 200km square cell
    - output :
        * [figure2.tiff](10-figures/figure3.tiff)
* [figure3.R](00-scripts/step5/figures/figure3.R) : from table of genetic,environmental,geographic variables by cell, this script generates figures of determinant of the patterns of fish genetic diversity.
    - inputs :
        * [total_data_genetic_diversity_with_all_descripteurs.tsv](09-all_descripteurs/total_data_genetic_diversity_with_all_descripteurs.tsv)
        * [grid_equalarea200km](01-infos/grid_equalarea200km) : shapefile of worldmap equal area projection epsg:4326 with nested equal area grids (cell sizes of 200km)
        * [EnvFreshwater.csv](01-infos/EnvFreshwater.csv) : slope and flow information for each geographical cell with a river 
        * [distanceCote](01-infos/distanceCote.csv) : distance from shore for each cell
    - output :
        * [figure3.pdf](10-figures/figure2.pdf)

* [figureS1.R](00-scripts/step5/supplementary_figures/figureS1.R) : from table of genetic,environmental,geographic variables by cell and shapefiles, it generates Spatial autocorrelogramme based on the I-Moran coefficient figure.
    - inputs :
        * [total_data_genetic_diversity_with_all_descripteurs.tsv](09-all_descripteurs/total_data_genetic_diversity_with_all_descripteurs.tsv) : table of center of cell (xy) coordinates, ID of cell, mean of Genetic diversity, number of species, mean/sd number of individuals by species into each cell, bathymetry, chlorophyll concentration, oxygen concentration, temperature, drainage basin surface area into each cell
        * [grid_equalarea200km](01-infos/grid_equalarea200km) : shapefile of worldmap equal area projection epsg:4326 with nested equal area grids (cell sizes of 200km)
    output :
        * [figureS1.pdf](10-figures/figureS1.pdf)
* [figureS2.R](00-scripts/step5/supplementary_figures/figureS2.R) : from table of genetic,environmental,geographic variables by cell, this script generates figure of global distribution of higher and lower percentiles.
    - inputs :
         * [total_data_genetic_diversity_with_all_descripteurs.tsv](09-all_descripteurs/total_data_genetic_diversity_with_all_descripteurs.tsv) : table of center of cell (xy) coordinates, ID of cell, mean of Genetic diversity, number of species, mean/sd number of individuals by species into each cell, bathymetry, chlorophyll concentration, oxygen concentration, temperature, drainage basin surface area into each cell
        * [grid_equalarea200km](01-infos/grid_equalarea200km) : shapefile of worldmap equal area projection epsg:4326 with nested equal area grids (cell sizes of 200km)
        * [ne_50m_land](01-infos/ne_50m_land) : shapefile of worldcoast
        * [ne_50m_rivers_lake_centerlines_scale_rank](01-infos/ne_50m_rivers_lake_centerlines_scale_rank) : shapefile of riverlines
        * [GSHHS_h_L2.shp](01-infos/big_lakes/GSHHS_h_L2.shp) : shape polygon file of big lakes
        * [figureS2.tiff](10-figures/figureS2.tiff)

* [figureS3.R](00-scripts/step5/supplementary_figures/figureS3.R) : it generates barplot of species diversity distribution accross latitudinal band (10°) and a boxplot of _marine_, _freshwater_ species diversity by cells.
    - inputs :
        * [total_data_genetic_diversity_with_all_descripteurs.tsv](09-all_descripteurs/total_data_genetic_diversity_with_all_descripteurs.tsv)
        * [grid_equalarea200km](01-infos/grid_equalarea200km) : shapefile of worldmap equal area projection epsg:4326 with nested equal area grids (cell sizes of 200km)
        * [equalarea_id_coordsCA_FWRS_MR_RS.csv](01-infos/equalarea_id_coordsCA_FWRS_MR_RS.csv) : Table of species diversity by 200km square cell
    - output :
        * [figureS4.pdf](10-figures/figureS3.pdf)

* [figureS4.R](00-scripts/step5/supplementary_figures/figureS4.R) : from table of genetic,environmental,geographic variables by cell, this script generates figure of regional effect on the global genetic diversity pattern.
    - input :
        * [total_data_genetic_diversity_with_all_descripteurs.tsv](09-all_descripteurs/total_data_genetic_diversity_with_all_descripteurs.tsv)
    - output :
        * [figureS4.pdf](10-figures/figureS3.pdf)

* [figureS5.R](00-scripts/step5/supplementary_figures/figureS5.R) : it generates maps of the number of species by cell, number of sequences by cell and number of sequences by species by cell.
    - inputs :
        * [total_data_genetic_diversity_with_all_descripteurs.tsv](09-all_descripteurs/total_data_genetic_diversity_with_all_descripteurs.tsv)
        * [grid_equalarea200km](01-infos/grid_equalarea200km)
        * [ne_50m_land](01-infos/ne_50m_land)
        * [ne_50m_rivers_lake_centerlines_scale_rank](01-infos/ne_50m_rivers_lake_centerlines_scale_rank)
        * [GSHHS_h_L2.shp](01-infos/big_lakes/GSHHS_h_L2.shp)
    - output : 
        * [figureS5.tiff](10-figures/figureS5.tiff)

* [figureS6.R](00-scripts/step5/supplementary_figures/figureS6.R) : barplot of sequences number and species number by taxonomic family/order
    - input :
        * [watertype_all_modeles_effectives_family.csv](11-sequences_taxonomy_habitat/watertype_all_modeles_effectives_family.csv)
    - output :
        * [figureS6.pdf](10-figures/figureS6.pdf)

* [figureS7.R](00-scripts/step5/supplementary_figures/figureS7.R) :  it generates maps of the taxonomic coverage by cells for _marine_ and _freshwater_ species
    - inputs :
        * [total_data_genetic_diversity_with_all_descripteurs.tsv](09-all_descripteurs/total_data_genetic_diversity_with_all_descripteurs.tsv)
        * [grid_equalarea200km](01-infos/grid_equalarea200km) : shapefile of worldmap equal area projection epsg:4326 with nested equal area grids (cell sizes of 200km)
        * [ne_50m_land](01-infos/ne_50m_land) : shapefile of worldcoast 
        * [ne_50m_rivers_lake_centerlines_scale_rank](01-infos/ne_50m_rivers_lake_centerlines_scale_rank) : shapefile of riverlines
        * [GSHHS_h_L2.shp](01-infos/big_lakes/GSHHS_h_L2.shp) : shape polygon file of big lakes
        * [equalarea_id_coordsCA_FWRS_MR_RS.csv](01-infos/equalarea_id_coordsCA_FWRS_MR_RS.csv) : Table of species diversity by 200km square cell
    - output :
        * [figureS7.tiff](10-figures/figureS7.tiff)

* [figureS8.R](00-scripts/step5/supplementary_figures/figureS8.R) :  it generates barplot of species diversity distribution accross latitudinal band (10°) at species level
    - inputs :
        * [marine_latbands_bootstraps.csv](08-genetic_diversity/marine_latbands_bootstraps.csv) : table of ID_latband,ISO3,is_sea,cloMeanVal,cloMinVal,cloMaxVal,bathyVal,AP,HDI_2015,fshD information attributed to each latitudinal latband with _marine_ species
        * [freshwater_latbands_bootstraps.csv](08-genetic_diversity/freshwater_latbands_bootstraps.csv) : table of ID_latband,ISO3,is_sea,cloMeanVal,cloMinVal,cloMaxVal,bathyVal,AP,HDI_2015,fshD information attributed to each latitudinal latband with _freshwater_ species
    - output :
        * [figureS8.pdf](10-figures/figureS8.pdf)


# 6 Taxonomy and habitat attributed to each individual sequences
#### BASH scripts
* [rename_family_bold_to_ncbi.sh](00-scripts/step6/rename_family_bold_to_ncbi.sh) : from cured table of individual sequences, cure family column by renaming [BOLD](http://www.boldsystems.org/) family by its equivalent into [NCBI taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy).
    - input :
        * [cured_sequences_withdemerpelag.csv](11-sequences_taxonomy_habitat/cured_sequences_withdemerpelag.csv) : cured table of individual sequences with habitat column
    - output :
        * [cured_family_sequences_withdemerpelag.csv](11-sequences_taxonomy_habitat/cured_family_sequences_withdemerpelag.csv) : cured taxonomy/family table of individual sequences with habitat column

#### PYTHON3 scripts
* [sequences_table.py](00-scripts/step6/sequences_table.py) : from table of genetic,environmental,geographic variables by cell, name of fasta files into [06-species_alnt_cluster](06-species_alnt_cluster) and CO1 sequences with lat/lon information, it writes a table of sequences with geographical cell localisation .
    - inputs :
        * [total_data_genetic_diversity_with_all_descripteurs.tsv](09-all_descripteurs/total_data_genetic_diversity_with_all_descripteurs.tsv) : table with cell as row and genetic,environmental,geographic variables as column
        * [06-species_alnt_cluster](06-species_alnt_cluster) : folder containing {species}.fasta files with {species} as [BOLD](http://www.boldsystems.org/)'s name of the species
        * [co1_ssll_seqbold_data.tsv](03-filtered_data/co1_ssll_seqbold_data.tsv) : table of fitlered CO1 sequences with lat/lon and [BOLD](http://www.boldsystems.org/)'s taxonomy information
    - output : 
        * [map_marine_sequences.csv](11-sequences_taxonomy_habitat/map_marine_sequences.csv) : table of individual sequences with geographical cell localisation
* [sequences_taxonomy.py](00-scripts/step6/sequences_taxonomy.py) : from cured table of individuals sequences and data table used for the different models, write a table of number of species/number of sequences by taxonomic order/family used for each model.
    - inputs :
        * [cured_family_sequences_withdemerpelag.csv](11-sequences_taxonomy_habitat/cured_family_sequences_withdemerpelag.csv) : cured taxonomy/family table of individual sequences with habitat column
        * [models](09-all_descripteurs/models) : folder which contains model's table of cells
    - output :
        * [watertype_all_modeles_effectives_family.csv](11-sequences_taxonomy_habitat/watertype_all_modeles_effectives_family.csv) : table of number of species/number of sequences by taxonomic order/family used for each model

#### R scripts
* [check_freshwater_assignation.R](00-scripts/step6/check_freshwater_assignation.R) : from model's table of cells and table of sequences, it writes a list of species with wrong watertype assignment according to the model using [rfishbase](https://cran.r-project.org/web/packages/rfishbase/rfishbase.pdf).
    - inputs :
        * [map_marine_sequences.csv](11-sequences_taxonomy_habitat/map_marine_sequences.csv) : table of individual sequences with geographical cell localisation
        * [freshwaterDG.txt](09-all_descripteurs/models/freshwaterDG.txt) : freshwater model's table of cells (model's table of cells are stored into folder [models](09-all_descripteurs/models))
    - output :
        * [wrong_freshwater_sequences.csv](11-sequences_taxonomy_habitat/wrong_freshwater_sequences.csv) : table of individual sequences with wrong watertype assignment
* [sequences_demerpelag.R](00-scripts/step6/sequences_demerpelag.R) : from table of sequences with geographical cell localisation, it assigns habitat (demersal, pelagic...) information based on species name attributed to the sequence and write a new table with habitat column.
    - input :
        * [map_marine_sequences.csv](11-sequences_taxonomy_habitat/map_marine_sequences.csv) : table of individual sequences with geographical cell localisation
    - output :
        * [sequences_withdemerpelag.csv](11-sequences_taxonomy_habitat/sequences_withdemerpelag.csv) :  table of individual sequences with habitat column
* [sequences_cure_species_name.R](00-scripts/step6/sequences_cure_species_name.R) : cure failed habitat assignment and wrong species name which are not recognized by [fishbase](https://www.fishbase.de/) database.
    - input :
        * [sequences_withdemerpelag.csv](11-sequences_taxonomy_habitat/sequences_withdemerpelag.csv) : table of individual sequences with habitat column
    - output :
        * [cured_sequences_withdemerpelag.csv](11-sequences_taxonomy_habitat/cured_sequences_withdemerpelag.csv) : cured table of individual sequences with habitat column
