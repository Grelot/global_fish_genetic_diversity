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
## according to a list of marine species, 
## moves the fasta and coords files into marine, freshwater repertories.
## cluster .coords and .fasta files into freshwater or marine 
## fisbase attribute of the species
##
##########################################################################
for ls_species in `ls ./05-species_alnt/*fasta`;
do
 species_name=`basename $ls_species .fasta`
 troncated_species=`basename $ls_species .fasta | sed 's/_/ /g' | awk '{ print $1}'`; 
 if grep "$troncated_species" ./01-infos/marine_actinopterygii_species.txt;
 #marine
 then
 absolutepathspeciesfasta=`readlink -f ./05-species_alnt/$species_name".fasta"`;
 ln -s $absolutepathspeciesfasta ./06-species_alnt_cluster/marine/$species_name".fasta";
 absolutepathspeciescoords=`readlink -f ./05-species_alnt/$species_name".coords"`;
 ln -s $absolutepathspeciescoords ./06-species_alnt_cluster/marine/$species_name".coords";
 #freshwater
 else
 absolutepathspeciesfasta=`readlink -f ./05-species_alnt/$species_name".fasta"`;
 ln -s $absolutepathspeciesfasta ./06-species_alnt_cluster/freshwater/$species_name".fasta";
 absolutepathspeciescoords=`readlink -f ./05-species_alnt/$species_name".coords"`;
 ln -s $absolutepathspeciescoords ./06-species_alnt_cluster/freshwater/$species_name".coords";
 fi
 ln -s $absolutepathspeciesfasta ./06-species_alnt_cluster/total/$species_name".fasta";
 ln -s $absolutepathspeciescoords ./06-species_alnt_cluster/total/$species_name".coords";
done

