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
## aligns sequences from the same species with MUSCLE
## and creates coordinates file for each sequence.
##
##########################################################################
#directory for .coords and .fasta files
#04-species_seq
#write .coords and .fasta files from CO1 georeferenced sequences file into 04-species_seq
python3 ./00-scripts/step2/fasta_coords_files_species_generator.py -o ./04-species_seq/ -f ./03-filtered_data/co1_ssll_seqbold_data.tsv

#directory for muscle3 intra-species sequences alignment results
#05-species_alnt
#align raw fasta sequences of the same species with MUSCLE3 then convert it into ODD-format fasta
for species_path in `ls ./04-species_seq/*fasta`;
do
  species_name=`basename $species_path .fasta | tr -d "."`;
  echo $species_name;
  muscle3 -in $species_path -out ./05-species_alnt/$species_name.tmp;
  sed -e 's/A/1/g' -e 's/C/2/'g -e 's/G/3/g' -e 's/T/4/g' -e 's/-/0/g' -e 's/>[0-9]\+/J/g' ./05-species_alnt/$species_name.tmp | tr -d '\n' | tr 'J' '\n' | tail -n +2 | sed -e 's/./& /g' -e 's/[ \t]$//g' > ./05-species_alnt/$species_name.fasta;
  printf "\n" >> ./05-species_alnt/$species_name.fasta;
done
rm ./05-species_alnt/*.tmp

#symbolic links from coords files to 05-species_alnt
for filecoords in `ls ./04-species_seq/*.coords`;
do
  linknamefilecoords=`basename $filecoords .coords | tr -d "."`;
  absolutepathfilecoords=`readlink -f $filecoords`;
  ln -s $absolutepathfilecoords ./05-species_alnt/$linknamefilecoords".coords";
done  

