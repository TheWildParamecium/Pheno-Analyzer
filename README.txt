This repository was made for keeping the code delevoped to build an automated workflow analysis of patient's data.
It was developed for a final project of a master's degree in Advanced Biotechnology.

It was titled:
Quality analysis of clinical information in patient cohorts for its use in clinical diagnostic support tools

And the main idea can be seen in the abstract:
The purpose of this work is to deﬁne translations of patients medical records to ontologies of computational   analysis such as the Human Phenotype Ontology (HPO). It will be look at how informative are the proﬁles obtained,
the existing correlation with genomic data and its applicability to support tools for clinical diagnosis.
For this purpose, an automatic worﬂow has been developed, in order to analyze the data of a cohort of patients
with Copy Number Variations (CNVs), a type of mutation that is believed to be responsible for cognitive
impairment, behavioral alterations, abnormalities in the cardiovascular system, etc. The workﬂow has two
distinct paths, one that analyzes the quality of the descriptions of the patient cohort based on three main
parameters: size of the proﬁles, speciﬁcity of the terms and phenotypic space used. With these parameters,
the quality of phenotypic proﬁles of the patients can be characterized. The other way consists in the
execution of a phenotype-genotype association analysis. For this purpose, signiﬁcant relationships between 
genomic variant region among several patients and anomalous phenotypes observed were searched. 
When signiﬁcant mutated region-phenotype relationships have been found, the genes involved are analyzed to 
determine their functions and their impact on the phenotype. It will also be assessed how the quality of the
phenotypic proﬁle information inﬂuences the results obtained.

The tool has two main parts.

1)Dataset-specific part: this part deals with the work of cleaning up the data, removing empty fields and generate
a file with a standard format, in order to make the rest of the tool automated. This part then is specific of 
the source of data collected. In this work the dataset used can be seen at input_files folder, in the files 
fenotipo.txt and genotipo.txt, which contains phenotypic and genotypic information of the patients dataset. 
This part will be expanded in the future in order to include other format's files or sources. The standard table
produced has the following format: patient ID, chromosome affected, start and stop position (as chromosome 
coordinates) of the mutation, and codes of the Human Phenotype Ontology defining the phenotypic profile of the
patient. The columns are divided by tabs (TSV format) and HPO codes are listed as comma-separated values with 
no space between them. If a patient has more than a mutation, then there will be more data with the same
 patient ID in separate rows (as much as number of mutations defined).

2)Automated analysis part: this part works always in the same way, despite of the dataset used to analyse.
 It uses combination.txt as the source for the analysis, needing to be in the standard format defined in the
 paragraph above. Then it downloads some needed files (HPO ontology and Human Reference Genome flat files) and
 performs the phenotype-genotype association analysis. Finally, it does the quality analysis of phenotypic 
 descriptions of the patient's dataset and puts all the data in a html report, making easy the visualization
 of the knowledge obtained through the analysis. 

The whole workflow (dataset-specific and automated analysis) is launch with main_launch.sh. It also can be found a relaunch used in the Master's project to analyse only patients with mutations at chromosome X (x_launch.sh). 
table_preprepare.sh script deals with dataset-specific part of the workflow, and analyse.sh with the automated 
analysis. It is needed to mention that in order to work with this tools, some R libraries (clusterProfiler, org.Hs.eg.db, optparse, stringr, ReactomePA, optparse, biomaRt) and ruby gems (NetAnalyzer, optparse, erb) are needed to be installed, and ruby gem NetAnalyzer needs ruby version 2.4.4 or below in order to work. Some of the 
R libraries can be found at general repository of CRAN and others at bioconductor. 

The whole workflow was developed and executed using an Ubuntu 18.04 console version, operated in Windows. This application can be found at Microsoft Store.

Soon a comprehensive description of the scripts of the tool will be made. 