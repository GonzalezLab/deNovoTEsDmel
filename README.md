# deNovoTEsDmel

## DESCRIPTION:
There are five scripts to determine the TEs detected by one, two (three different combinations) or three methods (REPET, TIDAL and TEMP)

***script_1.R***: for each method, this script concatenates TEs intervals detected in all strains and filter out those TEs in heterochromatin regions.  

***script_2.sh***: this script merges the overlapping TE intervals for each method (REPET, TIDAL, TEMP) while preserving the information about each TE and each strain.  

***script_3.R***: this script removes redundancy and assigns an ID for each TE and each method.  

***script_4.sh***: this script merges the overlapping TE intervals detected with the three methods while preserving the information about each TE, method and strain.  

***script_5.R***: this script finds, classifies and counts the number of TEs found by one, two and three methods. It needs the file Reference_tes.txt. This file should be placed in the same folder as the input/output files.  

## HOW TO RUN:
`Rscript --vanilla script_1.R </path_1/> </path_2/> </path_3/> </path_4/>`    
>path_1: path to the REPET input files   
>path_2: path to the TIDAL input files   
>path_3: path to the TEMP input files   
>path_4: path to the output files   

`./script_2.sh </path_1> </path_2/>`   
>path_1: path to the input/output files    
>path_2: path to sortBed and merge in bedtools utilities  

`Rscript --vanilla script_3.R </path_1/>`     
>path_1: path to the input/output files   

`./script_4.sh </path_1> </path_2/>`   
>path_1: path to the input/output files   
>path_2: path to sortBed and merge in bedtools utilities   

`Rscript --vanilla script_5.R </path_1/>`     
>path_1: path to the input/output files   
