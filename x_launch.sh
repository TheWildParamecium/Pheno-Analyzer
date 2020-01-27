#! /usr/bin/env bash
grep "chrX" ./input_files/genotipo.txt > ./input_files/xgenotipo.txt
current_dir=`pwd`
PATH=$current_dir/custom_scripts:$PATH
export PATH
#Automatization process for the first dataset used.
#Prepare the standard table with table_preprepare bash script. Parameters: genotipo.txt and fenotipo.txt tables as
#input and combination.txt as standard table filename output. The numbers 1 and 2 are the headers of the genotipo.txt
#and fenotipo.txt tables respectively
table_preprepare.sh xgenotipo.txt fenotipo.txt xcombination.txt 0 2
#Then analyse the standard table with analyse bash script. Parameters: combination.txt as the table filename input,
# 1 as the header of the table, and report.html as the html filename output. This will be changed later, because
#currently the ruby script is generating a json file with the results to be loaded in javascript.
analyse.sh xcombination.txt 1 xreport.html