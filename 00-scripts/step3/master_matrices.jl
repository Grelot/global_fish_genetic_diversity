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
## generates master data matrices from species sequences alignments.
##
##########################################################################

include("Lib_Create_Master_Matrices.jl")
for cluster in ["total", "freshwater", "marine"]
    create_master_matrices("./06-species_alnt_cluster/$(cluster)")
end
