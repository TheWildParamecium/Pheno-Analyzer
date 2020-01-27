This is a repository for keeping different versions of the master's degree final proyect.

Contains:

1) table_preprepare.sh : bash script for automatization of the tasks of preparing the input tables to be procesed by 
      get_standard_table.rb. It first remove unnecessary columns of genotype and fenotype talbes and then calls get_standard_table.rb with       genotype and phenotype tables, number of header rows for each table and the name of the output table
      as parameters. It returns a combination of both tables in a standard format required later for analyse.sh. This script and   
      get_standard_table.rb are specific of the input tables used this time, and different tables may require another preparation scripts. 

2)get_standard_table.rb : ruby script that deals with the tasks of combining input tables (genotype and fenotype tables), and store it
     in a suitable format required later for working with analyse_patients.rb. 

3)analyse.sh :  bash script for automatization of the tasks of analysing the standard table created by table_prepepare.sh. It calls
     analyse_patients.rb with the standard table filename, the number of header rows and the name of the html_file as parameters, and 
     initially it returned a html_file with patients' dataset information displayed with html tables, but will be deprecated, as now it
     return a json file that will give the data to a javascript that will put the information in the html file. This script and
     analyse_patients.rb can be used with whatever table that has be saved in the standard format (see combination.txt for more details).

4)analyse_patients.rb : ruby script that getting a standard table as input (that is, in a [id,[[chromosome, mutation_start,                    mutation_stop], [another_chromosome, mutation_start, mutation_stop],hp_terms(phenotypes),hp_codes], required for a good data parsing) 
     it returns a bunch of information (like number of patient, number of abnormal phenotypes, number of phenotypes per patient, number of
     chromosomes affected per patient and frequencies of which chromosomes are affected) initally in a html_file, now in a json file, that
     later is moved to reporte_html, which shows the information in a enjoyable view.
     
5)reporte_html/ : directory where the patient's dataset information is showed in report.html, that load the data saved in results.json by
      using the script.js contained. Css stylesheets has been used to show the data in a enjoyable format. 
      
6)Other data that has been used by the different script, as the input data and the different files created by the scripts.

7)tablas_ejemplo1.sh initializes the tasks flow by calling table_preprepare.sh and analyse.sh with the required parameters.
