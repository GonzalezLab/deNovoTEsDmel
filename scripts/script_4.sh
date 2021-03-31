#!/bin/bash

## This will merge overlapping TE intervals for the three methods (REPET, TIDAL, TEMP) preserving the info about each TE and each strain

#pathinout: path to the input/output files (i.e., /Users/home)
#pathbedtools: path to sortBed and merge (i.e., /Users/home/bedtools2-master/bin/)

pathinout=$1
pathbedtools=$2

#Remove header 
sed -e '1d' $pathinout/all_tes_non_redundant_indexed.txt >  $pathinout/all_tes_non_redundant_indexed.1.txt

#Sort all file 
cd $pathbedtools
./sortBed -i $pathinout/all_tes_non_redundant_indexed.1.txt > $pathinout/all_tes_non_redundant_indexed_sorted.txt

#Merge all overlapping intervals preserving the info with the ID
./bedtools merge -i $pathinout/all_tes_non_redundant_indexed_sorted.txt -d 1 -c 7 -o collapse > $pathinout/all_tes_non_redundant_indexed_sorted_bedtools.txt


				