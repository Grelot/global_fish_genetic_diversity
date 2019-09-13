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
## Keep only the CO1 sequences with lat/lon information
##
##########################################################################
#convert raw data into tsv format
#sed 's/ /_/g' ./02-raw_data/seqbold_data.txt > ./02-raw_data/seqbold_data.tsv
#remove samples with no lat/lon information
awk '{ if( $34 != "_" && $35 != "_") print $0 }' ./02-raw_data/seqbold_data.tsv > ./03-filtered_data/ll_seqbold_data.tsv
#remove samples with no species name
awk '{ if($21 != "_") print $0 }' ./03-filtered_data/ll_seqbold_data.tsv > ./03-filtered_data/sll_seqbold_data.tsv
#remove samples with no sequences or sequences with IUAPC ambiguity i.e. polybase characters (e.g. "N")
awk '{ test= match( $45,"[RWSNYMKHBVD_]"); if( test ==0) print $0}' ./03-filtered_data/sll_seqbold_data.tsv > ./03-filtered_data/ssll_seqbold_data.tsv
#get information lat/lon with geonames for data with missing lat/lon information but region
./00-scripts/step1/get_geonames_coordinates.sh
#merge geonames-georeferenced sequences with already georeferenced sequences
cat ./03-filtered_data/new_geo_actinopterygii_sl ./03-filtered_data/ssll_seqbold_data.tsv > ./03-filtered_data/total_ssll_seqbold.tsv

#get locus name for each sample
awk '{if(NR>1)print}' ./03-filtered_data/total_ssll_seqbold.tsv | while read ligne;
do
process_id=`echo $ligne | awk '{print $1}'`
requete="http://www.boldsystems.org/index.php/Public_RecordView?processid="$process_id
#echo $process_id
locus_id=`lynx -dump $requete | grep Locus:`
echo $process_id" "$locus_id
done > ./03-filtered_data/locushead

#get sample's ID with a locus which is not "Cytochrome Oxidase Subunit 1 5' Region"
cat ./03-filtered_data/locushead | while read ligne;
do
  locus_name=`echo $ligne | cut -f 3- -d " "`
  locus_ID=`echo $ligne | cut -f 1 -d " "`
  #echo $locus_name
  if [ "$locus_name" != "Cytochrome Oxidase Subunit 1 5' Region" ]
  then
    echo $locus_ID
  fi
done > ./03-filtered_data/wronglocusID

#remove samples with locus sequence that is not "Cytochrome Oxidase Subunit 1 5' Region"
grep -avf ./03-filtered_data/wronglocusID ./03-filtered_data/total_ssll_seqbold.tsv > ./03-filtered_data/co1_ssll_seqbold_data.tsv


