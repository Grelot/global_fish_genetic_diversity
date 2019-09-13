# -*- coding: utf-8 -*-
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
# NOTICE
#==============================================================================
# converts from DMS format to DD format the given coordinates.
#==============================================================================
# MODULES
#==============================================================================
import re
import argparse
#==============================================================================
# FUNCTIONS
#==============================================================================
def dms2dd(degrees, minutes, seconds, direction):
    dd = round(float(degrees) + float(minutes)/60 + float(seconds)/(60*60),3);
    if direction == 'S' or direction == 'W':
        dd *= -1
    return dd;

def dd2dms(deg):
    d = int(deg)
    md = abs(deg - d) * 60
    m = int(md)
    sd = (md - m) * 60
    return [d, m, sd]

def parse_dms(dms):
    #format : N 45° 24' 5''	W 73° 48' 50''
    #dms = dms.replace("''","' ")
    parts = re.split('[^\d\w]+', dms)
    #print parts
    lat = dms2dd(parts[1], parts[2], parts[3], parts[0])
    lng = dms2dd(parts[5], parts[6], parts[7], parts[4])
    return (lat, lng)

'''
dd = parse_dms("36°57'9' N 110°4'21' W")
print(dd)
print(dd2dms(dd[0]))
'''

#==============================================================================
#ARGUMENTS
#==============================================================================

parser = argparse.ArgumentParser(description='remove AAA* suffix from fasta file')
parser.add_argument("-dms","--dms_coo",type=str)
parser.add_argument("-dd","--dd_coo", type=str)

args = parser.parse_args()
dms_coo = args.dms_coo
dd_coo = args.dd_coo


#==============================================================================
#MAIN
#==============================================================================

#dans le cas ou on donne les coos en DMS

dd_converted_coo = parse_dms(dms_coo)
print(dd_converted_coo[0],dd_converted_coo[1],sep=" ")
#print(dms2dd(dd_converted_coo[0]))
