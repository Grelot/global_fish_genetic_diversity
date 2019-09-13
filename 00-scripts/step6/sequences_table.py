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
# write a table of sequences with geographical cell localisation
#==============================================================================
# MODULES
#==============================================================================
import difflib
import os
import re
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
## object "phylum"
class Phylum:
    def __init__(self,spec_name,genus,family,order):
        self.spec_name = spec_name
        self.genus = genus
        self.family = family
        self.order = order

#==============================================================================
#ARGUMENTS
#==============================================================================

#fichier="09-all_descripteurs/total_data_genetic_diversity_with_all_descripteurs.tsv"
#folder="06-species_alnt_cluster"
#rawf="../pipeline_propre/03-filtered_data/co1_ssll_seqbold_data.tsv"

## cell descripteurs table
fichier="09-all_descripteurs/total_data_genetic_diversity_with_all_descripteurs.tsv"
## folder containing fasta and coords files of each cell
folder="06-species_alnt_cluster"
## table of filtered sequences from BOLD
rawf="03-filtered_data/co1_ssll_seqbold_data.tsv"
## output csv file
outputf="11-sequences_taxonomy_habitat/map_marine_sequences.csv"

#==============================================================================
#MAIN
#==============================================================================

## list of the cells used for model
csv.field_size_limit(100000000)
fd=open(fichier,'r')
listLieu=[]
for ligne in fd.readlines()[1:]:
    listLieu.append(ligne.split("\t")[2])

## dictionnary [species]={genus,family,order} from filtered BOLD table
dicPhylum= dict()
with open(rawf,'r',encoding="ISO-8859-1") as fr:
    for ligne in fr.readlines():
        print(ligne)
        ligneSplit=re.sub(r' +','\t',ligne).split("\t")
        #ligneSplit=re.split(r'[ \t]', ligne[0])
        l_spec=ligneSplit[20]
        l_genus=ligneSplit[18]
        l_family=ligneSplit[14]
        l_order=ligneSplit[12]
        #print(l_spec,l_genus,l_family,l_order)
        #print(ligne)
        dicPhylum[l_spec]=Phylum(l_spec,l_genus,l_family,l_order)

## read intermediates files to extract species and location of each sequences
listSequence= []
listFichiers = dict()
listFichiers["marine"]=[]
listFichiers["freshwater"]=[]
listEquf = dict()
counter_seq=1
for watertype in ["marine","freshwater"]:
    foldert=folder+"/"+watertype
    for root, dirs, files in os.walk(foldert):
        for filename in files:
            listFichiers[watertype].append(filename)
    listEquf[watertype]=[name for name in listFichiers[watertype] if 'equalareacoords' in name]
    for equf in listEquf[watertype]:
        l_species=equf.split(".")[0]
        l_spec_split=l_species.split('_')
        keys=difflib.get_close_matches(l_species,list(dicPhylum.keys()))
        key=keys[0]       
        l_cured=dicPhylum[key].spec_name
        l_genus=dicPhylum[key].genus
        l_family=dicPhylum[key].family
        l_order=dicPhylum[key].order
        l_open=open(foldert+"/"+equf,'r')
        for ligne in l_open.readlines():
            ligneSplit=ligne.split("\t")
            if len(ligneSplit) > 2:
                l_lat=ligneSplit[0]
                l_lon=ligneSplit[1]
                l_lieu=ligneSplit[2]
                #print(counter_seq,watertype,l_species,l_cured,l_genus,l_family,l_order,l_lat,l_lon,l_lieu)
                listSequence.append(Sequence(counter_seq,watertype,l_species,l_cured,l_genus,l_family,l_order,l_lat,l_lon,l_lieu))
                counter_seq+=1
            else:
                print(l_open)

## check species name (optionnal)
listSpecies=[]
for i in listSequence:
    if "," in i.cured_name:
        print(i.species_name)
    if i.cured_name not in listSpecies:
        listSpecies.append(i.cured_name)

## print tables
file_strings=""
for seq in listSequence:
    local_string=",".join([str(seq.id),seq.watertype,seq.species_name,seq.cured_name,seq.genus,seq.family,seq.order,seq.lat,seq.lon,seq.lieu])
    file_strings+=local_string

# write it
with open(outputf, 'w') as csvfile:
    writer = csv.writer(csvfile)
    csvfile.write(file_strings)

