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
## using rfishbase, assign habitat (demersal, pelagic...)
## information to each individual sequences
## according to their attributed species name
##
##
##########################################################################
## Libraries
library(rfishbase)

##########################################################################
## load data
## input table of sequences
sequencesf="11-sequences_taxonomy_habitat/map_marine_sequences.csv"
## output table of sequences with habitat column
outputf="11-sequences_taxonomy_habitat/sequences_withdemerpelag.csv"

### total raw sequence
sequences=read.table(sequencesf, header=F, sep=",")
colnames(sequences)=c("id","watertype","species","curedspecies","genus","family","order","lat","lon","cell")
## list of species names
messpecies=as.character(gsub('_', ' ', sequences$curedspecies))
## habitat assignation according to the specie name
test=species(messpecies,fields = species_fields$habitat)
## missing data
#length(which(is.na(test$DemersPelag)))

##########################################################################
## write sequences table with taxonomy information and HABITAT ASSIGNATION
sequences_withDemerpelag=cbind(sequences, test$DemersPelag)
colnames(sequences_withDemerpelag)=c("id","watertype","boldspecies","curedspecies","genus","family","order","lat","lon","cell","demerpelag")
write.table(sequences_withDemerpelag,outputf,row.names=F,sep=",",col.names=T)

##########################################################################
## see the species without habitat assignation which require curation
species_withoutDemerpelag=unique(sort(sequences_withDemerpelag[which(is.na(sequences_withDemerpelag$demerpelag)),]$curedspecies))
