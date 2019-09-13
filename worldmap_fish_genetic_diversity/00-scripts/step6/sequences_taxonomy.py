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
## 
#==============================================================================
# NOTICE
#==============================================================================
# write a table of number of species/number of sequences
# by taxonomic order/family used for each model
#==============================================================================
# MODULES
#==============================================================================
from ete3 import Tree, TreeStyle, NodeStyle, TextFace, NCBITaxa
import re
import numpy as np
#import matplotlib.pyplot as plt
import csv


#==============================================================================
# CLASSES
#==============================================================================
## object "sequence"
class Sequence:
    def __init__(self, id, watertype, species_name,cured_name, genus,family, order, lat, lon, lieu):
        self.id = id
        self.watertype = watertype
        self.species_name = species_name
        self.cured_name = cured_name
        self.genus = genus
        self.family = family
        self.order = order
        self.lat= lat
        self.lon = lon
        self.lieu = lieu
"""
object "cells used by a model"
this object is defined by :
 - the file which stores information used to generate the model
 - a list of cells used by the model
 - a list of individual sequences into the cells used by the model
 - a list of taxonomic families with corresponding number of species used by the model
 - a list of taxonomic families with corresponding number of sequences used by the model
 - the number of individual sequences by family
""" 
class CellModele:
    def __init__(self,id,nom_fichier,watertype):
        self.watertype= watertype
        self.id=id
        self.fichier=nom_fichier
        self.cellList=[]
        self.seqList=[]
        self.numberOfSpecFamDic={}
        self.numberOfSeqFamDic={}


#==============================================================================
#FUNCTIONS
#==============================================================================
## extract the column "colons" of a csv file "fs" and return it as a list
def colon_fichier(fs,colons,seps):
    with open(fs,'r') as f:
        reader = csv.DictReader(f, delimiter=seps)
        listeCell=[]
        for row in reader:
            listeCell.append(row[colons])
    return listeCell

## calculate de number of BOLD species into a family
def numbSpecFam(lseqList,lnomFam):
    speciesFam=[]
    for seq in lseqList:
            if seq.family == lnomFam:
                if seq.cured_name not in speciesFam:
                    speciesFam.append(seq.cured_name)
    numberofspecFam=len(speciesFam)
    return numberofspecFam


#==============================================================================
#ARGUMENTS
#==============================================================================
## output table of number of species/number of sequences by taxonomic order/family used for each model
outputf="11-sequences_taxonomy_habitat/watertype_all_modeles_effectives_family.csv"

## input table of individual sequences with cured family taxonomy name
#fichier="map_marine_sequences.csv"
#fichier="../worldmap_fish_genetic_diversity/00-scripts/step5/review_pe/cured_sequences_withdemerpelag.csv"
fichier="11-sequences_taxonomy_habitat/cured_family_sequences_withdemerpelag.csv"

## extraire les lieux conserves pour les differents modeles

## list of name of files of each modele
fmodeles=["09-all_descripteurs/models/dcM.514Marine.txt",
"09-all_descripteurs/models/dcF346Freshwater.txt",
"09-all_descripteurs/models/S2I3C66M.txt",
"09-all_descripteurs/models/S2I3C75F.txt",
"09-all_descripteurs/models/S2I4C29M.txt",
"09-all_descripteurs/models/S2I5C32F.txt",
"09-all_descripteurs/models/S3I2C74M.txt",
"09-all_descripteurs/models/S3I2C75F.txt",
"09-all_descripteurs/models/S8I2C36M.txt",
"09-all_descripteurs/models/S8I2C34F.txt"]


#==============================================================================
#MAIN
#==============================================================================
## create a list of modeles with attribute "id" and "fichier" (name of file) and "watertype"
mod_watertype=["marine","freshwater"]*5
listCellMod=[]
idcount=0
for nomf in fmodeles:
    idcount+=1
    listCellMod.append(CellModele(idcount,nomf,mod_watertype[idcount-1]))

## attribute a list of cells for each model
for cemo in listCellMod:
    cellMo=colon_fichier(cemo.fichier,'cell',' ')
    cemo.cellList=cellMo

## sequences information
fs=open(fichier,'r')
sequenceList=[]
for ligne in fs.readlines()[1:]:
    ligneSplit=ligne.split(',')
    lid=ligneSplit[0]
    lwatertype=ligneSplit[1].replace('"', '')
    lspecname=ligneSplit[2].replace('"', '')
    lcuredname=ligneSplit[3].replace('"', '')
    lgenus=ligneSplit[4].replace('"', '')
    lfamily=ligneSplit[5].replace('"', '')
    lorder=ligneSplit[6].replace('"', '')
    llat=ligneSplit[7]
    llon=ligneSplit[8]
    llieu=ligneSplit[9]
    sequenceList.append(Sequence(lid,lwatertype,lspecname,lcuredname,lgenus,lfamily,lorder,llat,llon,llieu))

uSpec=np.unique([i.cured_name for i in sequenceList])
uFam=np.unique([i.family for i in sequenceList])

## select sequences of models
for seq in sequenceList:
    for cemo in listCellMod:
        if seq.lieu in cemo.cellList and seq.watertype == cemo.watertype:
            cemo.seqList.append(seq)


## all the ncbi tree taxonomy
ncbi = NCBITaxa()

## get rank information of nodes of the phylogenetic tree "actinopterygii"
actinodesc=ncbi.get_descendant_taxa('Actinopterygii',intermediate_nodes=True,rank_limit='species',collapse_subspecies=True)
actino=ncbi.translate_to_names(actinodesc)
actinorank=ncbi.get_rank(actinodesc)

### FAMILY
## extract family names into the class "actinopterygii"
actinofamilies=[]
for key, value in actinorank.items():
    if value == "family":
    	actinofamilies.append(key)

## calculate number of NCBI species into a family
ncbiActinoFamDic={}
for idFam in actinofamilies:
    nomFam=ncbi.translate_to_names([idFam])[0]
    nbSpecIntoFam=len(ncbi.get_descendant_taxa(idFam))
    ncbiActinoFamDic[nomFam]=nbSpecIntoFam

## calculate number of BOLD species into a family
seqActinoFamDic={}
for idFam in actinofamilies:
    nomFam=ncbi.translate_to_names([idFam])[0]
    print(nomFam)
    seqActinoFamDic[nomFam]=numbSpecFam(sequenceList,nomFam)
    for cemo in listCellMod:
        cemo.numberOfSpecFamDic[nomFam]=numbSpecFam(cemo.seqList,nomFam)

## find family BOLD name that are not in NCBI taxonomy
wrongBOLDFamily=[]
for seq in sequenceList:
    cc=0
    for idFam in actinofamilies:
        nomFam=ncbi.translate_to_names([idFam])[0]
        if seq.family == nomFam:
            cc=1
            break
    if cc !=1:
        if seq.family not in wrongBOLDFamily:
            wrongBOLDFamily.append(seq.family)

## calculate the number of sequences into a family
def numberofSequencesFam(lseqList, lnomFam):
    countSequences=0
    for seq in lseqList:
        if seq.family == lnomFam:
            countSequences+=1
    return countSequences

## number of sequences by family
NumberOfSeqActinoFamDic={}
for idFam in actinofamilies:
    nomFam=ncbi.translate_to_names([idFam])[0]
    print(nomFam)
    NumberOfSeqActinoFamDic[nomFam]=numberofSequencesFam(sequenceList,nomFam)
    for cemo in listCellMod:
        cemo.numberOfSeqFamDic[nomFam]=numberofSequencesFam(cemo.seqList,nomFam)

## concatenate results into a string
my_str="taxonOrder,taxonFamily,ncbiNumberOfSpecies,boldNumberOfSpecies,boldNumberOfSequences"
my_str+=","
my_str+=",".join([cemo.fichier.split(".")[0]+"NumberOfSpecies" for cemo in listCellMod])
my_str+=","
my_str+=",".join([cemo.fichier.split(".")[0]+"NumberOfSequences" for cemo in listCellMod])
my_str+="\n"
for idFam in actinofamilies:
    nomFam=ncbi.translate_to_names([idFam])[0]
    for seq in sequenceList:
        if nomFam == seq.family:
            nomOrder = str(seq.order)
            break
    my_str+=",".join((nomOrder,nomFam,str(ncbiActinoFamDic[nomFam]),str(seqActinoFamDic[nomFam]),str(NumberOfSeqActinoFamDic[nomFam])))
    my_str+=","
    my_str+=",".join([str(cemo.numberOfSpecFamDic[nomFam]) for cemo in listCellMod])
    my_str+=","
    my_str+=",".join([str(cemo.numberOfSeqFamDic[nomFam]) for cemo in listCellMod])
    my_str+="\n"

## write output files table of number of species/number of sequences
## by taxonomic order/family for each model
with open(outputf, 'w') as csvfile:
    writer = csv.writer(csvfile)
    csvfile.write(my_str)


#==============================================================================
# END OF SCRIPTS
#==============================================================================
## this part of the code is not necessary
"""
N = len(ncbiActinoFamDic)
ind = np.arange(N)    # the x locations for the groups
width = 0.05       # the width of the bars: can also be len(x) sequence

p1 = plt.bar(ind,ncbiActinoFamDic.values())
p2 = plt.bar(ind,seqActinoFamDic.values())

plt.ylabel('Number of species')
plt.title('')
plt.xticks(ind,seqActinoFamDic.keys() )
plt.yticks(np.arange(0, 81, 10))
plt.show()

### plot Number of sequences by family
N = len(ncbiActinoFamDic)
ind = np.arange(N)    # the x locations for the groups
width = 0.05       # the width of the bars: can also be len(x) sequence
plt.figure(1)
ax = plt.subplot(111)
ax.bar(ind,np.repeat(0,len(ind)), color="black")
ax.bar(ind,NumberOfSeqActinoFamDic.values(), color="black")
plt.ylabel('Number of sequences')
ax.title('')
ax.set_xticks(ind)
ax.set_xticklabels(tuple(seqActinoFamDic.keys()), rotation=90)
for xs, lab in zip(ax.get_xticklabels(),seqActinoFamDic.keys()):
    xs.set_fontsize(5)
    if seqActinoFamDic[lab] != 0:
        xs.set_color('black')        
    else:
        xs.set_color('lightgrey')

plt.yticks(np.arange(0,max(NumberOfSeqActinoFamDic.values()) , 100))
plt.subplots_adjust(bottom=0.2)
plt.show()



### ORDER
actinoOrders=[]
for key, value in actinorank.items():
    if value == "order":
    	actinoOrders.append(key)

## number of species by order according to NCBI
ncbiActinoOrdDic={}
for idOrd in actinoOrders:
    nomOrd=ncbi.translate_to_names([idOrd])[0]
    nbSpecIntoOrd=len(ncbi.get_descendant_taxa(idOrd))
    ncbiActinoOrdDic[nomOrd]=nbSpecIntoOrd



## number of species by order
seqActinoOrdDic={}
for idOrd in actinoOrders:
    nomOrd=ncbi.translate_to_names([idOrd])[0]
    print(nomOrd)
    speciesOrd=[]
    for seq in sequenceList:
        if seq.order == nomOrd:
            if seq.cured_name not in speciesOrd:
                speciesOrd.append(seq.cured_name)
    numberofspecOrd=len(speciesOrd)
    print(numberofspecOrd)
    seqActinoOrdDic[nomOrd]=numberofspecOrd

## number of sequences by order
NumberOfSeqActinoOrdDic={}
for idOrd in actinoOrders:
    nomOrd=ncbi.translate_to_names([idOrd])[0]
    print(nomOrd)    
    countSeq=0
    for seq in sequenceList:
        if seq.order == nomOrd:
            countSeq+=1
    NumberOfSeqActinoOrdDic[nomOrd]=countSeq

### write table order number species number seq number seq ncbi

my_str="taxonOrder,ncbiNumberOfSpecies,boldNumberOfSpecies,boldNumberOfSequences\n"
for idOrd in actinoOrders:
    nomOrd=ncbi.translate_to_names([idOrd])[0]
    my_str+=",".join((nomOrd,str(ncbiActinoOrdDic[nomOrd]),str(seqActinoOrdDic[nomOrd]),str(NumberOfSeqActinoOrdDic[nomOrd])))
    my_str+="\n"

with open('effectives_order.csv', 'w') as csvfile:
    writer = csv.writer(csvfile)
    csvfile.write(my_str)




my_str="taxonOrder,taxonFamily,ncbiNumberOfSpecies,boldNumberOfSpecies,boldNumberOfSequences\n"
for idFam in actinofamilies:
    nomFam=ncbi.translate_to_names([idFam])[0]
    for seq in sequenceList:
        if nomFam == seq.family:
            nomOrder = str(seq.order)
            break
    my_str+=",".join((nomOrder,nomFam,str(ncbiActinoFamDic[nomFam]),str(seqActinoFamDic[nomFam]),str(NumberOfSeqActinoFamDic[nomFam])))
    my_str+="\n"


with open('effectives_family.csv', 'w') as csvfile:
    writer = csv.writer(csvfile)
    csvfile.write(my_str)





### plot Number of species by order
N = len(ncbiActinoOrdDic)
ind = np.arange(N)    # the x locations for the groups
width = 0.05       # the width of the bars: can also be len(x) sequence
plt.figure(1)
ax = plt.subplot(111)
ax.bar(ind,ncbiActinoOrdDic.values(), color="black")
ax.bar(ind,seqActinoOrdDic.values(), color="red")
plt.ylabel('Number of species')
ax.title('')
ax.set_xticks(ind)
ax.set_xticklabels(tuple(seqActinoOrdDic.keys()), rotation=90)
for xs, lab in zip(ax.get_xticklabels(),seqActinoOrdDic.keys()):
    if seqActinoOrdDic[lab] != 0:
        xs.set_color('black')
    else:
        xs.set_color('lightgrey')

plt.yticks(np.arange(0,max(ncbiActinoOrdDic.values()) , 100))
plt.subplots_adjust(bottom=0.2)
plt.show()
plt.savefig("taxon_order_coverage.png")

### plot Number of sequences by order
N = len(ncbiActinoOrdDic)
ind = np.arange(N)    # the x locations for the groups
width = 0.05       # the width of the bars: can also be len(x) sequence
plt.figure(1)
ax = plt.subplot(111)
ax.bar(ind,np.repeat(0,len(ind)), color="black")
ax.bar(ind,NumberOfSeqActinoOrdDic.values(), color="black")
plt.ylabel('Number of sequences')
ax.title('')
ax.set_xticks(ind)
ax.set_xticklabels(tuple(seqActinoOrdDic.keys()), rotation=90)
for xs, lab in zip(ax.get_xticklabels(),seqActinoOrdDic.keys()):
    if seqActinoOrdDic[lab] != 0:
        xs.set_color('black')
    else:
        xs.set_color('lightgrey')

plt.yticks(np.arange(0,max(NumberOfSeqActinoOrdDic.values()) , 200))
plt.subplots_adjust(bottom=0.2)
plt.show()
"""