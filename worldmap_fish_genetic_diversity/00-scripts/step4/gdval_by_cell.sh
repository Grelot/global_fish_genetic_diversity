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
## generates CSV files with 2 columns :
## cell ID and mean genetic diversity per species into the cell
##
##########################################################################
for CATEGORIE in `echo "freshwater marine"`;
do
echo "site;GD" > 08-genetic_diversity/gdval_by_area_$CATEGORIE".csv"
awk 'NR>1' 08-genetic_diversity/$CATEGORIE"_equalarea_numbers.csv" | while read ligne;
do
echo $ligne | awk -F ',' '{ print $1";"$2}'
done >> 08-genetic_diversity/gdval_by_area_$CATEGORIE".csv"
done


