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
## This code is based on :
## An Anthropocene Map of Genetic Diversity
##
## Andreia Miraldo, Sen Li, Michael K. Borregaard,
## Alexander Floréz-Rodriguéz, Shyam Gopalakrishnan, Mirneza Risvanovic,
## Zhiheng Wang, Carsten Rahbek, Katharine A. Marske & David Nogués-Bravo
##
## Submitted to Science, 2016
## Code in this file by Michael K. Borregaard
## modified by P.E. GUERIN
## Montpellier 2017
#########################################################################
##
## generates files for statistical analysis at next step : 
## mean genetic diversity per cell, genetic diversity per species per cell,
## number of individuals per species, number of species per cell, 
## cell coordinates, cell ID
##
##
#########################################################################

include("Lib_GD_summary_functions.jl")
using DataFrames, CSV
for cluster in ["freshwater", "marine"]
    ## metrics by area (coords, number of individuals, number of species,GDval, mean, sd, ...)
    mm =ready_data("07-master_matrices/$(cluster)_pairwise_equalarea.csv")
    coord = CSV.read("01-infos/equalarea_id_coords.tsv",delim="\t")
    names!(coord, [:site,:left,:right,:bottom,:top])
    total = calcGDmeanmednumberind(mm)
    totcoo = join(coord, total, on = :site)
    #writetable("site_coords_GD_nb_species_indv_marine.csv", totcoo)
    ## keep only areas with IUCN referenced species presence
    iucn_site = CSV.read("08-genetic_diversity/gdval_by_area_$(cluster).csv",header=1,delim=';')
    names!(iucn_site,[:site,:GD])
    test=join(iucn_site,totcoo, on = :site)
    tes = DataFrame(site=test[:site], left=test[:left], top=test[:top], right=test[:right],bottom=test[:bottom], GD_mean=test[:GD_mean],GD_median=test[:GD_median], nb_species=test[:nb_species],nb_indv_mean=test[:nb_indv_mean],nb_indv_median=test[:nb_indv_median],nb_indv_sd=test[:nb_indv_sd])
    CSV.write("08-genetic_diversity/metrics_by_area_$(cluster).csv", tes)
    ## metrics by species on area (coords, number of individuals, GDvaln mean, sed, ...)
    total2 = calcGDspeciespersite(mm)
    totcoo2 = join(coord, total2, on = :site)
    #writetable("site_coords_GD_species_per_site_marine.csv", totcoo2)
    ## keep only areas with IUCN referenced species presence    
    test=join(iucn_site,totcoo2, on = :site)
    tes2 = DataFrame(species=test[:species],
    site=test[:site],
    left=test[:left],
    top=test[:top],
    right=test[:right],
    bottom=test[:bottom],
    nb_indv=test[:nb_indv],
    GD_species=test[:GDspecies],
    GD_area=test[:GD])
    CSV.write("08-genetic_diversity/metrics_by_area-species_$(cluster).csv", tes2)
end
