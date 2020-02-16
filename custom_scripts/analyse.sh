#! /usr/bin/env bash
current_dir=`pwd`
PATH=$current_dir/custom_scripts:$PATH
export PATH

#folder variables
reportfolder=./reporte_html
interfilesfolder=./intermediate_files

#Downloads files needed for analysis (HPO ontology and gff Genome flat files)
if [ -f $interfilesfolder/hp.obo ]; then
    echo "Requested files(hp.obo) already exist, so skipping proccess of download"
else
    wget https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/master/hp.obo -O $interfilesfolder/hp.obo
fi
if [ -f $interfilesfolder/ref_GRCh37.p13_top_level.gff3 ]; then
    echo "Requested files(ref_GRCh37.p13_top_level.gff3) already exist, so skipping proccess of download"
else
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/ARCHIVE/ANNOTATION_RELEASE.105/GFF/ref_GRCh37.p13_top_level.gff3.gz 
    unzip ref_GRCh37.p13_top_level.gff3.gz -d $interfilesfolder/
fi

#Making the report folder (it will contain the report of patients in html format)
mkdir -p $reportfolder

#Getting GVRs. It searches common CNVs in different patients for creating GVRs. Then outputs a text file 
#with pairwise relationships of the type GVRiD - phenotype_term_1, needed for NetAnalyzer to work.
get_gvr_stats_score.rb -s $interfilesfolder/$1 -S $2 -T $interfilesfolder/Tripartite_network.txt -G $interfilesfolder/GVR_annots.txt

#Using NetAnalyzer(a ruby gem) to obtain statistically significant relationships between GVRs and phenotypes
LAYERS='fenotipo,HP:;genotipo,chr;paciente,ID'
NetAnalyzer.rb -i $interfilesfolder/Tripartite_network.txt -l $LAYERS -m 'hypergeometric' -u 'fenotipo,genotipo;paciente' -a $interfilesfolder/hypergeometric_values.txt -N

#Finding/Mapping genes inside each GVR and producing a text file with that genes for a gene enrichment process
gene_finder.rb -v $interfilesfolder/hypergeometric_values.txt -G $interfilesfolder/GVR_annots.txt \
    -d $4 -p $5 -f $interfilesfolder/ref_GRCh37.p13_top_level.gff3 \
    -c $interfilesfolder/candidate_genes.txt -H $interfilesfolder/hp.obo

#Gen enrichment with R. It gets some gene data fields like GO terms, description, etc, using a query in BioMart.
#It returns gene data enriched as output. It is not needed currently the workflow, as the later analysis 
#using R libraries can perform without that data. But it could be useful for future approaches so it has
#been kept.
gene_enrichment.r --candidategenes $interfilesfolder/candidate_genes.txt \
    --enrichedgenes $interfilesfolder/enriched_genes.txt \

#Merge all gene enriched data (this means, fusion GVR data with enriched gen data using GVRid and GeneID 
#as a join point). It outputs this enriched data to a text file used later in analyse_patients script. 
#It also outputs a file used in GO, KEGG and Reactome over-representation analysis (also called functional
#analysis)
fusion_gene_data.rb -e $interfilesfolder/enriched_genes.txt \
    -f $interfilesfolder/ref_GRCh37.p13_top_level.gff3 \
    -v $interfilesfolder/hypergeometric_values.txt \
    -G $interfilesfolder/GVR_annots.txt \
    -d $4 -p $5 -t $interfilesfolder/pre_geneset_data.txt \
    -m $interfilesfolder/gvr_complement_data_for_geneset.txt

#Performing functional analysis. It outputs the results produced, with statistically
#significant functional terms attached to a certain GVR, if found.
#Gene Ontology analysis of the Biological Proccess subontology 
gene_set_analysis.r --pregeneset $interfilesfolder/pre_geneset_data.txt \
    --results $interfilesfolder/geneset_results_go_bp.txt \
    --analysis GO --subontology BP
#Gene Ontology analysis of the Molecular Function subontology
gene_set_analysis.r --pregeneset $interfilesfolder/pre_geneset_data.txt \
    --results $interfilesfolder/geneset_results_go_mf.txt \
    --analysis GO --subontology MF
#Functional analysis of KEGG database 
gene_set_analysis.r --pregeneset $interfilesfolder/pre_geneset_data.txt \
    --results $interfilesfolder/geneset_results_kegg.txt \
    --analysis KEGG
#Functional analysis of Reactome database
gene_set_analysis.r --pregeneset $interfilesfolder/pre_geneset_data.txt \
    --results $interfilesfolder/geneset_results_reactome.txt \
    --analysis reactome

#Assessing the quality of phenotypic profile of patient's dataset. It also integrates the results
#produced by functional analysis. Finally it produces a html file as a report with the patient's dataset results.
analyse_patients.rb -s $interfilesfolder/$1 -S $2 \
    -H $interfilesfolder/hp.obo \
    -o $reportfolder/results.json \
    -R $reportfolder/$3 -t templates/patient_report.erb \
    -v $interfilesfolder/hypergeometric_values.txt \
    -G $interfilesfolder/GVR_annots.txt -d $4 -p $5 \
    -z $interfilesfolder/geneset_results_go_bp.txt \
    -Z $interfilesfolder/geneset_results_go_mf.txt \
    -x $interfilesfolder/geneset_results_kegg.txt \
    -X $interfilesfolder/geneset_results_reactome.txt \
    -m $interfilesfolder/gvr_complement_data_for_geneset.txt \
    -w $interfilesfolder/supl_pats_data.text