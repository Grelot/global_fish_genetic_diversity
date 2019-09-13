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
## functions used to generate figures
## Supplementary Figure S6.
## Taxonomic coverage of the sequences used by the model
## (a) Number of species per order 
## (b) Number of sequences per order 
##
##
##########################################################################
## Libraries
## no library required 

##########################################################################
## functions

## sum of number of element by 2 columns moF moM
## according to a category defined by the third column moCategory
aggregateBy2Colons <- function(moF,moM,moCategory){
	numberSpecModT=aggregate(cbind(moF,moM), by=list(Category=moCategory), FUN=sum)
 	numbSpecMoSumColons=data.frame(category=numberSpecModT$Category, numberOfElements=numberSpecModT$moM+numberSpecModT$moF)
 	return(numbSpecMoSumColons)
}


##########################################################################
## load data
tc=read.csv("11-sequences_taxonomy_habitat/watertype_all_modeles_effectives_family.csv",header=T,sep=",")
## correct wrong size of family in NCBI taxonomy
tc[which(tc$taxonFamily == "Citharinidae"),]$ncbiNumberOfSpecies = tc[which(tc$taxonFamily == "Citharinidae"),]$ncbiNumberOfSpecies + 91
tc[which(tc$taxonFamily == "Gadidae"),]$ncbiNumberOfSpecies = 53

## number of species by family in two model(marine+freshwater)
monModelMARINENumbSpec=tc$dcMNumberOfSpecies
monModelFRESHWATERNumbSpec=tc$dcF346FreshwaterNumberOfSpecies
## number of sequences by family in two model(marine+freshwater)
monModelMARINENumbSeq=tc$dcMNumberOfSequences
monModelFRESHWATERNumbSeq=tc$dcF346FreshwaterNumberOfSequences


## models with at least 4 sequences by species and at least 2 species by cell
#monModelMARINENumbSpec=tc$S2I4C29MNumberOfSpecies
#monModelFRESHWATERNumbSpec=tc$S2I5C32FNumberOfSpecies
## models with at least 3 sequences by species and at least 2 species by cell
#monModelMARINENumbSpec=tc$S2I3C66MNumberOfSpecies
#monModelFRESHWATERNumbSpec=tc$S2I3C75FNumberOfSpecies
## models with at least 2 sequences by species and at least 8 species by cell
#monModelMARINENumbSpec=tc$S8I2C36MNumberOfSpecies
#monModelFRESHWATERNumbSpec=tc$S8I2C34FNumberOfSpecies
## models with at least 2 sequences by species and at least 3 species by cell
#monModelMARINENumbSpec=tc$S2I3C66MNumberOfSpecies
#monModelFRESHWATERNumbSpec=tc$S3I2C75FNumberOfSpecies

## number of species by family in two model(marine+freshwater)
NSpeF<-aggregateBy2Colons(monModelMARINENumbSpec,monModelFRESHWATERNumbSpec,tc$taxonFamily)
TotSpeFamilyNumbSeq<-aggregate(tc$ncbiNumberOfSpecies, by=list(Category=tc$taxonFamily), FUN=sum)
propFam=data.frame(proportion=(NSpeF$numberOfElements/TotSpeFamilyNumbSeq$x)*100, family=as.vector(TotSpeFamilyNumbSeq$Category))

## number of species by ORDER in two model(marine+freshwater)
NSpeO<-aggregateBy2Colons(monModelMARINENumbSpec,monModelFRESHWATERNumbSpec,tc$taxonOrder)
## Number of species order total
TotSpeOrderNumbSeq<-aggregate(tc$ncbiNumberOfSpecies, by=list(Category=tc$taxonOrder), FUN=sum)
propOrder=data.frame(proportion=(NSpeO$numberOfElements/TotSpeOrderNumbSeq$x)*100,
	order=as.vector(TotSpeOrderNumbSeq$Category))

## number of sequences by family in two model(marine+freshwater)
NSeqF<-aggregateBy2Colons(monModelMARINENumbSeq,monModelFRESHWATERNumbSeq,tc$taxonFamily)

## number of sequences by order in two model(marine+freshwater)
NSeqO<-aggregateBy2Colons(monModelMARINENumbSeq,monModelFRESHWATERNumbSeq,tc$taxonOrder)


##########################################################################
## Quartile of the proportion of species and number of sequences
## for the 63 orders and 480 families
summary(propFam)
summary(propOrder)


##########################################################################
## barplots
## Number of species per family
#barplot(propFam$proportion, las=2,names.arg=propFam$family,ylab="Prop species per order",main ="full models",cex.names = 0.5 )
## Number of species per order 

## Number of sequences by family
#barplot(NSeqF$numberOfElements, las=2,names.arg=NSeqF$category,ylab="Number sequences per order",cex.names = 0.5 )
## Number of sequences by order


##########################################################################
## write pdf files
pdf("10-figures/figureS6.pdf",width=12,height=12,paper='special')

layout(matrix(c(1,2),nrow=2))
par(mar=c(6,4,4,4)) # c(bottom, left, top, right)
barplot(propOrder$proportion,las=2,names.arg=propOrder$order,ylab="Proportion species per order",cex.names = 0.6 )
title("(a)", adj  = 0)
barplot(NSeqO$numberOfElements, las=2,names.arg=NSeqO$category,ylab="Number sequences per order",cex.names = 0.6 )
title("(b)", adj  = 0)

dev.off()
