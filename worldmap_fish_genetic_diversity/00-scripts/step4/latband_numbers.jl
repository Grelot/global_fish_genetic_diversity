##########################################################################
## Codes for the paper:
## Global determinants of freshwater and marine fish genetic diversity
## Authors :
## Stephanie Manel, Pierre-Edouard Guerin, David Mouillot,
## Simon Blanchet, Laure Velez, Camille Albouy, Loic Pellissier
## 
## Montpellier, 2017-2019
## Submited to Nature communications, 2019
##
##########################################################################
##
## attributes mean genetic diversity at each latitudinal band.
## Genetic diversity is calculated from master data matrices
## Also performs bootstrap on species for each latitudinal band
##
##########################################################################
include("Lib_GD_summary_functions.jl")
include("Lib_bootstrap.jl")

using DataFrames, CSV

for cluster in ["freshwater", "marine"]
    mm =ready_data("./07-master_matrices/$(cluster)_pairwise_latbands.csv")
    ## mean GD by latband
    latband, latband_full = calc_cellvalues(mm)
    ## GD by species by latband and number of available sequences
    latband_gd_species=calcGDspeciespersite(mm)
    ## filter species with GD=0
    latband_gd_species1=filter(row -> row[:GDspecies] !=0, latband_gd_species)   
    ## filter species with less than 3 sequences
    latband_gd_species2=filter(row -> row[:nb_indv] > 2, latband_gd_species1) 
    ## bootstraping with 1000 replicats
    latband_boot = bootstrap_species(latband_gd_species2, 1000)
    ## write CSV
    CSV.write("./08-genetic_diversity/$(cluster)_latbands_numbers.csv", latband)
    CSV.write("./08-genetic_diversity/$(cluster)_latbands_bootstraps.csv", latband_boot)
end

