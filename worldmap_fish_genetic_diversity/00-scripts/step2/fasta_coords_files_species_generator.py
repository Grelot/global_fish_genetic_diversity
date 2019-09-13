#==============================================================================
# Codes for the paper:
# Global determinants of freshwater and marine fish genetic diversity
# Authors :
# Stephanie Manel, Pierre-Edouard Guerin, David Mouillot,
# Simon Blanchet, Laure Velez, Camille Albouy, Loic Pellissier
# 
# Montpellier, 2017-2019
# Submited to Nature communications, 2019
#
#==============================================================================
#NOTICE
#==============================================================================

'''
extracts sequences and associated coordinates from the filtered data

create 2 files from the database file co1_ssll.tsv
for each specie i :
    specie_i.coords
    specie_i.fasta

specie_i.coords description :
(indv 1) lat lon
(indv 2) lat lon
...

specie_i.fasta description :
(indv 1) aligned_seq1
( indv 2) aligned_seq2
...


aligned_seq description :
alignment produced by muscle3
FORMAT: 0 1 2 3 4 where 0=gap 1=A 2=C 3=G 4=T

'''

#==============================================================================
#MODULES
#==============================================================================
import re
import sys
import subprocess
import argparse
import os.path

#==============================================================================
#CLASS
#==============================================================================

class Individu:
    def __init__(self,lat,lon,seq):
        self.lat = lat
        self.lon = lon
        self.seq = seq

class Species:
    def __init__(self,name):
        self.name = name
        self.listOfIndv = []
        self.counter = 0
    def indv_append(self,lat,lon,seq):
        self.listOfIndv.append(Individu(lat,lon,seq))


#==============================================================================
#FUNCTIONS
#==============================================================================




#==============================================================================
#ARGUMENTS
#==============================================================================

parser = argparse.ArgumentParser(description='specie seq geo')
parser.add_argument("-o","--output", type=str)
parser.add_argument("-f","--inputFile",type=str)


#==============================================================================
#MAIN
#==============================================================================
args = parser.parse_args()
resultsPath = args.output
inputFile = args.inputFile

dicOfSpecies={}

#read input file and create species and indv objects
with open(inputFile,'r') as readFile:
    for ligne in readFile.readlines()[1:]:
        ligneSplit=ligne.split()
        speciesName=ligneSplit[20]
        indvLat=ligneSplit[33]
        indvLon=ligneSplit[34]
        indvSeq=ligneSplit[44]
        if speciesName not in dicOfSpecies:
            dicOfSpecies[speciesName] = Species(speciesName)
        dicOfSpecies[speciesName].indv_append(indvLat,indvLon,indvSeq)
        dicOfSpecies[speciesName].counter+=1
readFile.close()

#write species.coords and species.fasta files into resultsPath
for key, val in dicOfSpecies.items():
    if val.counter > 1:
        with open("{0}/{1}.coords".format(resultsPath, key),'w') as coordsFile:
            coordsText=""
            for indv in val.listOfIndv:
                coordsText+="{0}\t{1}\n".format(indv.lat,indv.lon)
            coordsFile.write(coordsText)
        coordsFile.close()
        with open("{0}/{1}.fasta".format(resultsPath, key),'w') as fastaFile:
            fastaText=""
            id_indv=1
            for indv in val.listOfIndv:
                fastaText+=">{0}\n".format(id_indv)
                fastaText+="{0}\n".format(indv.seq)
                id_indv+=1
            fastaFile.write(fastaText)
        fastaFile.close()
    else:
        print(val.name, val.listOfIndv, "only one indiv", sep=" ")
#fin du script
