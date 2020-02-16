#! /usr/bin/env bash
current_dir=`pwd`
PATH=$current_dir/custom_scripts:$PATH
export PATH

#folder variables
initfilesfolder=./input_files
interfilesfolder=./intermediate_files

#Initializes the folders (if needed) where intermediate files will be stored
mkdir -p intermediate_files

#Getting table columns of interest of the genotype's table. It also orders the new table by id column
cut -f 1,3,4,5 $initfilesfolder/$1 | sort -n -k 1 > $interfilesfolder/genotabla.txt

#Getting table columns of interest of the phenotype's table
cut -f 1,2,6-25 $initfilesfolder/$2 > $interfilesfolder/fenotabla.txt

#Excecutes ruby script in order to get the standard table
get_standard_table.rb -g $interfilesfolder/genotabla.txt -G $4 -p $interfilesfolder/fenotabla.txt -P $5 -o $interfilesfolder/$3 -s $interfilesfolder/supl_pats_data.text
