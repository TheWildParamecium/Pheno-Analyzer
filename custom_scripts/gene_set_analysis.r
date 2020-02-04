#!/usr/bin/env Rscript
#LOADING LIBRARIES
library("clusterProfiler")
library("org.Hs.eg.db")
library("optparse")
library("stringr")
library("ReactomePA")
#FUNCTION DEFINITIONS
#MAIN
#CREATING OPTPARSE OPTIONS
option_list = list( 
    make_option(c("-r", "--results"), action="store_true",  
        type="character", default="../intermediate_files/gene_set_results.txt",
        help="Save the results produced by gene set analysis. By default [default]"),
    make_option(c("-t", "--pregeneset"), action="store_true",  
        type="character", default="../intermediate_files/GVR_annots.txt",
        help="Load the enriched genes' file for gene set analysis. By default [default]"),
    make_option(c("-s", "--subontology"), action="store_true",  
        type="character", default="BP",
        help="It sets the the subontology of GO to be used. Supports BP, MF and CC. By default [default]"),
    make_option(c("-a", "--analysis"), action="store_true",  
        type="character", default="GO",
        help="It sets the ontology to be used. Supports KEGG, reactome and GO. By default [default]")
        )
opt = parse_args(OptionParser(option_list=option_list))

#Loading genes file and stratifying by GVRs
genes = read.delim(opt$pregeneset)
names(genes) = c("GVR", "ENTREZ", "SYMBOL")
gvrs = lapply(split.data.frame(genes, genes$GVR), function(x) x$ENTREZ)

if(opt$analysis == "GO"){
#Applying gene set analysis for every GVR if the case is GO
results = list()
for (name in names(gvrs)) {
    results[[name]] = enrichGO(gene = as.character(gvrs[[name]]),
                OrgDb         = org.Hs.eg.db,
                ont           = opt$subontology,
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01)
}
}else if(opt$analysis == "KEGG"){
#Applying gene set analysis for every GVR if the case is KEGG
results = list()
for (name in names(gvrs)) {
    results[[name]] = enrichKEGG(gene = as.character(gvrs[[name]]),
                      organism         = 'human',
                      pvalueCutoff  = 0.01)
}
}else if(opt$analysis == "reactome"){
#Applying gene set analysis for every GVR if the case is Reactome
results = list()
for (name in names(gvrs)) {
    results[[name]] = enrichPathway(gene = as.character(gvrs[[name]]),
                      organism         = 'human',
                      pvalueCutoff  = 0.01)
}
}

#Saving data in a txt file
txt = c()
for (name in names(results)) {
    if (length(results[[name]]$Description) != 0){
        #txt <- c(txt, stringr::str_interp("${name}\t${results[[name]]$ID}\t${results[[name]]$GeneRatio}\t${results[[name]]$BgRatio}\t${results[[name]]$pvalue}\t${results[[name]]$Description}\t${results[[name]]$geneID}\n"))
        txt <- c(txt, stringr::str_interp("${name}\t${results[[name]]$Description}\n")) 
    } 
}
writeLines(txt, opt$results)