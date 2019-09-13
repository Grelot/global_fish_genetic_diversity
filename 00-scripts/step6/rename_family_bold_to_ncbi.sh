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
## rename BOLD family by their equivalent into NCBI taxonomy
##
##########################################################################

## input table of indivual sequences with taxonomy, geographical cell and habitat
FICHIER="11-sequences_taxonomy_habitat/cured_sequences_withdemerpelag.csv"
## output table of individual sequences with cured family taxonomy name
OUTPUTF="11-sequences_taxonomy_habitat/cured_family_sequences_withdemerpelag.csv"

## replace BOLD family name by NCBI family name
## /!\ we group subfamily Distichodontidae with family-level Citharinidae
sed -e 's/Anoplogastridae/Anoplogasteridae/g' \
-e 's/"Scaridae/Labridae/g' \
-e 's/Scomberesocidae/Belonidae/g' \
-e 's/Achiropsettidae/Rhombosoleidae/g' \
-e 's/Caesionidae/Lutjanidae/g' \
-e 's/Latidae/Centropomidae/g' \
-e 's/Rivulidae/Aplocheilidae/g' \
-e 's/Eleginopsidae/Eleginopidae/g' \
-e 's/Zeniontidae/Zenionidae/g' \
-e 's/Horabagridae/Bagridae/g' \
-e 's/Distichodontidae/Citharinidae/g' \
-e 's/Cynolebiidae/Aplocheilidae/g' \
$FICHIER > $OUTPUTF

