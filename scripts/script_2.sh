#!/bin/bash

##This will merge the overlapping TE intervals for each method (REPET, TIDAL, TEMP) preserving the info about each TE and each strain

#pathinout: path to the input/output files (i.e., /Users/home)
#pathbedtools: path to sortBed and merge (i.e., /Users/home/bedtools2-master/bin/)

pathinout=$1
pathbedtools=$2

##REPET
#Remove header
sed -e '1d' $pathinout/11_strains_repet_concatenated.txt > $pathinout/11_strains_repet_concatenated_indexed.1.txt

#Sort file 
cd $pathbedtools
./sortBed -i $pathinout/11_strains_repet_concatenated_indexed.1.txt > $pathinout/11_strains_repet_concatenated_indexed_sorted.txt

#Merge all overlapping intervals preserving the info with the ID
./bedtools merge -i $pathinout/11_strains_repet_concatenated_indexed_sorted.txt -d 1 -c 24,25 -o collapse > /$pathinout/11_strains_repet_concatenated_indexed_sorted_bedtools.txt

##TIDAL
#Remove header				
sed -e '1d' $pathinout/11_strains_tidal_concatenated.txt > $pathinout/11_strains_tidal_concatenated_indexed.1.txt

#Sort file 
cd $pathbedtools
./sortBed -i $pathinout/11_strains_tidal_concatenated_indexed.1.txt > $pathinout/11_strains_tidal_concatenated_indexed_sorted.txt

#Merge all overlapping intervals preserving the info with the ID
./bedtools merge -i $pathinout/11_strains_tidal_concatenated_indexed_sorted.txt -d 1 -c 18,19 -o collapse > $pathinout/11_strains_tidal_concatenated_indexed_sorted_bedtools.txt

##TEMP
#Remove header				
sed -e '1d' $pathinout/11_strains_temp_concatenated.txt > $pathinout/11_strains_temp_concatenated_indexed.1.txt

#Sort file 
cd $pathbedtools
./sortBed -i $pathinout/11_strains_temp_concatenated_indexed.1.txt > $pathinout/11_strains_temp_concatenated_indexed_sorted.txt

#Merge all overlapping intervals preserving the info with the ID
./bedtools merge -i $pathinout/11_strains_temp_concatenated_indexed_sorted.txt -d 1 -c 14,15 -o collapse > $pathinout/11_strains_temp_concatenated_indexed_sorted_bedtools.txt
