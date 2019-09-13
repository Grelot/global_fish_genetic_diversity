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
## check if the watertype [marine|freshwater] for each species by cell
## is correct according to the model [marine|freshwater]
##
##########################################################################
## Libraries
library(rfishbase)

##########################################################################
## load data
## input table of sequences
sequencesf="11-sequences_taxonomy_habitat/map_marine_sequences.csv"
## input cells of the freshwater model
freshwaterModelCellFile="09-all_descripteurs/models/freshwaterDG.txt"
## output csv file with list of species with wrong watertype assignment
wrongfreshwaterSeqFile="11-sequences_taxonomy_habitat/wrong_freshwater_sequences.csv"

# liste des especes d'actinopterygien
actino = species_list(Class="Actinopterygii")
cactino= country(actino)
## exctract "saltwater"=="marine" species
saltactino=data.frame(species=cactino$Species,saltwater=as.integer(cactino$Saltwater))
saltactino.uniq=unique(saltactino[order(saltactino$species),])
actino.saltwater=saltactino.uniq[which(saltactino.uniq$saltwater ==1),]

## data used for the model freshwater
#dcm=read.table("dcM.ind.339.txt",header=T)
dcf=read.table(freshwaterModelCellFile,header=T)

### total raw sequence
sequences=read.table(sequencesf,header=F,sep=",")
colnames(sequences)=c("id","watertype","species","curedspecies","genus","family","order","lat","lon","cell")

##########################################################################
## raw sequence kept for the model
seqmod=sequences[which(sequences$cell %in% dcf$cell),]
seqmodF=seqmod[which(seqmod$watertype=="freshwater"),]

## presumed marine species from the sample sequence
uniqspecies=unique(sort(gsub('_', ' ', seqmodF$species)))
## unknow marine species
unknow=uniqspecies[which(!uniqspecies %in% actino)]

wrongFreshwaterSpec=c(
"Paradiplogrammus bairdi",
"Pseudoicichthys australis",
"Trachinocephalus myops"
)

wrongFreshwaterSpec=c(
"Paradiplogrammus bairdi",
"Pseudoicichthys australis",
"Pseudopentaceros richardsoni",
"Pseudotrematomus bernacchii",
"Pseudotrematomus eulepidotus",
"Pseudotrematomus hansoni",
"Pseudotrematomus lepidorhinus",
"Pseudotrematomus nicolai",
"Pseudotrematomus pennellii",
"Pseudotrematomus scotti",
"Pseudotrematomus tokarevi",
"Trachinocephalus myops",
"Verulux cypselurus")

## write output csv file with list of species with wrong watertype assignment
wrongFreshwaterSeq=seqmod[which(seqmod$species %in% gsub(' ', '_', wrongFreshwaterSpec)),]
write.table(wrongFreshwaterSeq,wrongfreshwaterSeqFile,row.names=F,sep=",",col.names=T)

## (optionnal) visualize geographical cells with wrong-assigned watertype individuals
library(maps)
map("world", fill=TRUE, col="white", bg="lightblue", ylim=c(-90, 90), mar=c(0,0,0,0))
points(wrongFreshwaterSeq$lat,wrongFreshwaterSeq$lon, col="red", pch=16)


