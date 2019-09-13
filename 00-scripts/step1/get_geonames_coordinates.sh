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
## Uses [http://www.geonames.org/] to find 
## missing coordinates of individual sequences 
## from their textual information of location.
##
##########################################################################
#get samples with no lat/lon information but region name
awk '{if ( ($34=="_" || $35=="_") && ($40!="_" ))  print $0}' ./02-raw_data/seqbold_data.tsv > ./03-filtered_data/no_geo_region

#sed -e 's/[\,.:;$~)(=*%?!-]//g'
cat ./03-filtered_data/no_geo_region | while read ligne;
do
deb_ligne=`echo $ligne | cut -f 1-33 -d " "`
fin_ligne=`echo $ligne | cut -f 36- -d " "`
recherche=`echo $ligne | awk '{print $38_$39_$40}' | sed 's/_/+/g' | sed -e 's/[\,.:;$~)(=*%?!-]/+/g' | sed -e 's/\([^[:blank:]]\)\([[:upper:]]\)/\1+\2/g' | sed -e 's/+$//g'`
requete="http://www.geonames.org/search.html?q="$recherche
#echo $requete
result=`lynx -dump $requete | grep "Lat/Lng" | head -1 | awk '{ print $3"\t"$5}' | sed 's/\[10\]//g'`
#echo $result
if [[ -z $result ]]
then
  #echo "zut"
  lat=`lynx -dump $requetelynx -dump $requete | egrep -o "[SN].[0-9]+°.[0-9]+'.[0-9]+''" | head -1 `
  #echo $lat
  lng=`lynx -dump $requetelynx -dump $requete | egrep -o "[WE].[0-9]+°.[0-9]+'.[0-9]+''" | head -1 `
  #echo $lng
  if [[ ( -z $lat ) || ( -z $lng ) ]]
  then
    :
  else
    result=`python3 00-scripts/step1/lat_long_DMS_DD_converter.py -dms "$lat $lng"`
    echo -e "$deb_ligne $result $fin_ligne"
  fi
else
  result=`echo $result | sed -e 's/\t/ /g' `
  echo -e "$deb_ligne $result $fin_ligne"
fi
done > ./03-filtered_data/new_geo_actinopterygii
#remove entries with no species name
awk '{ if($21 != "_") print $0 }' ./03-filtered_data/new_geo_actinopterygii > ./03-filtered_data/new_geo_actinopterygii_s
#remove samples with no sequences or remove sequences with IUAPC ambiguity i.e. polybase characters (e.g. "N")
awk '{ test= match( $45,"[RWSNYMKHBVD_]"); if( test ==0) print $0}' ./03-filtered_data/new_geo_actinopterygii_s > ./03-filtered_data/new_geo_actinopterygii_sl

