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
## attributes mean genetic diversity at each equal area grid cell.
## Genetic diversity is calculated from master data matrices
##
##########################################################################
include("Lib_GD_summary_functions.jl")
using DataFrames, CSV

for cluster in ["freshwater", "marine"]
    mm =ready_data("./07-master_matrices/$(cluster)_pairwise_equalarea.csv")
    equalarea, equalarea_full = calc_cellvalues(mm)
    CSV.write("./08-genetic_diversity/$(cluster)_equalarea_numbers.csv", equalarea)
end


