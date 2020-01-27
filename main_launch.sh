#! /usr/bin/env bash
current_dir=`pwd`
PATH=$current_dir/custom_scripts:$PATH
export PATH
#Automatization process for the first dataset used.
#Prepare the standard table with table_preprepare bash script. Parameters: genotipo.txt and fenotipo.txt tables as
#input and combination.txt as standard table filename output. The numbers 1 and 2 are the headers of the genotipo.txt
#and fenotipo.txt tables respectively
table_preprepare.sh genotipo.txt fenotipo.txt combination.txt 1 2
#Then analyse the standard table with analyse bash script. Parameters: combination.txt as the table filename input,
# 1 as the header of the table, and report.html as the html filename output. 2 as the hypergeometric value and
#10 as the threshold of patients needed to take account of a GVR.
analyse.sh combination.txt 1 report.html 2 10