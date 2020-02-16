#!/usr/bin/env Rscript

#LOADING REQUIRED LIBRARIES
library("optparse")
library("biomaRt")

#MAIN
#CREATING OPTPARSE OPTIONS (parameters that will be received from bash script) 
option_list = list( 
    make_option(c("-c", "--candidategenes"), action="store_true",  
        type="character", default="../intermediate_files/candidate_genes.txt",
        help="Load the genes data file used for enrichment. By default [default]"),
    make_option(c("-e", "--enrichedgenes"), action="store_true",  
        type="character", default="../intermediate_files/enriched_genes.txt",
        help="Save the enriched genes in a file. By default [default]") )
opt = parse_args(OptionParser(option_list=option_list))

#Loading genes that will be given more information (of GO terms, description, etc)
candidate_genes_list = read.delim(opt$candidategenes)
candidate_genes = as.vector(unlist(candidate_genes_list))

#Loading BiomaRt data
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
#Performing a query with our genes
enriched_genes = getBM(attributes = c('hgnc_symbol', 'entrezgene_id', 'description',
                    'go_id', 'name_1006', 'kegg_enzyme'),
                        filters = 'hgnc_symbol', 
                        values = candidate_genes, 
                        mart = ensembl)
#Writting results to a text file                        
write.table(enriched_genes, file=opt$enrichedgenes, quote=FALSE, sep='\t', row.names=F)