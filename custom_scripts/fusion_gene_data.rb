#! /usr/bin/env ruby

#LOADING LIBRARIES
require 'optparse'

#FUNCTION DEFINITIONS
#Loading GVR with statistically significant phenotypes associated, produced by NetAnalyzer.
#It also loads his hypergeometric value (it measures the strengh of the bond) 
def load_gvr_and_hyper_values(gvr_path, hyper_path, threshold)
    gvr_annots = Hash.new { |hash, key| hash[key] = {'start' => nil, 'stop' => nil, 'chr' => nil, 'n_pats' => nil, 'phens' => [] } }

    File.open(gvr_path).each do |line|
        fields = line.chomp.split("\t")
        gvr_annots[fields[0]] = {'start' => fields[1], 'stop' => fields[2], 'chr' => fields[3], 'n_pats' => fields[4], 'phens' => [] }
    end

    File.open(hyper_path).each do |line|
        fields = line.chomp.split("\t")
        if fields[2].to_f >= threshold
            gvr_annots[fields[1]]['phens'].push(fields[0])
        end
    end
    return gvr_annots
end

#Loading the file with genomic annotations of the Human Genome, needed for mapping genes inside GVR coordinates
def load_gff_file(gff_path)
    gen_annots = {}
    chr = ""
    File.open(gff_path).each do |line|
        if !( line.start_with?("#") )
            fields = line.chomp.split("\t")
            if fields[2] == "region"
                metadata = fields[8].split(";")
                n = (( metadata[3] ).split("=")[1])
                chr = "chr#{n}"
            elsif fields[2] == "gene"
                metadata = fields[8].split(";")
                id = (metadata.shift.split("="))[1]
                gen_annots[id] = {}
                gen_annots[id]["pseudo"] = "false"
                gen_annots[id]["chr"] = chr
                gen_annots[id]["start"] = fields[3]
                gen_annots[id]["stop"] = fields[4]

                metadata.each do |item|
                    key = item.split("=")[0]
                    value = item.split("=")[1]
                    gen_annots[id][key] = value
                end
            end
        end
    end
    return gen_annots
end

#Loading the enriched genes' file produced by gene_enrichment.r script.
def load_enriched_genes(enriched_path)
    enriched_genes_data = {}
    counter = 0
    headers = []
    names = []
    File.open(enriched_path).each do |line|
        if counter == 0
            headers = line.chomp.split("\t")
            headers.shift
            counter += 1
        else
            fields = line.chomp.split("\t")
            name = fields.shift

            if !(names.include?(name)) 
                names.push(name)
                enriched_genes_data[name] = {}
                enriched_genes_data[name][headers[0]] = fields[0]
                enriched_genes_data[name][headers[1]] = [fields[1]]
                enriched_genes_data[name][headers[2]] = [fields[2]]
                enriched_genes_data[name][headers[3]] = [fields[3]]
                enriched_genes_data[name][headers[4]] = []
                if (fields.length == 5)
                    enriched_genes_data[name][headers[4]] = [fields[4]]
                end
            else
                if !(enriched_genes_data[name][headers[1]].include?(fields[1]))
                    enriched_genes_data[name][headers[1]].push(fields[1])
                end
                
                if !(enriched_genes_data[name][headers[2]].include?(fields[2]))
                    enriched_genes_data[name][headers[2]].push(fields[2])
                end
            
                if !(enriched_genes_data[name][headers[3]].include?(fields[3]))
                    enriched_genes_data[name][headers[3]].push(fields[3])
                end

                if (fields.length == 5)
                    if !(enriched_genes_data[name][headers[4]].include?(fields[4]))
                        enriched_genes_data[name][headers[4]].push(fields[4])
                    end
                end
            end
        end
    end
    return enriched_genes_data
end

#Merging all data together
def fusion_gvr_and_genes(gvrs_data, full_genes_data)
    gvrs_data_with_genes = gvrs_data.dup
    gvrs_data.each do |gvr_id, gvr_fields|
        gvrs_data_with_genes[gvr_id]["genes"] = {}
        full_genes_data.each do |gen_name, gen_fields|
            if gen_fields["chr"] == gvr_fields["chr"]
                if (gen_fields["start"] > gvr_fields["start"]) && (gen_fields["start"] < gvr_fields["stop"])
                    gvrs_data_with_genes[gvr_id]["genes"][gen_name] = gen_fields 
                elsif (gen_fields["stop"] > gvr_fields["start"]) && (gen_fields["stop"] < gvr_fields["stop"])
                    gvrs_data_with_genes[gvr_id]["genes"][gen_name] = gen_fields
                end
            end
        end
    end
    return gvrs_data_with_genes
end

#Setting enriched genes info to the genes inside each GVR
def enrich_genes(filtered_genes, enriched_genes)
    full_info_genes = {}
    enriched_genes.each do |gen_name, values|
        filtered_genes.each do |gen_id, metadata|
            if metadata["Name"] == gen_name
                full_info_genes[gen_name] = values
                full_info_genes[gen_name]["id"] = gen_id
                metadata.each do |key, value|
                    full_info_genes[gen_name][key] = value
                end
            end
        end
    end
    return full_info_genes
end

#MAIN
#Defining the script parser options (parameters that will be received from bash script) 
options = {}
OptionParser.new do |opts|

    options[:gvr_annots] = nil
    opts.on( '-G', '--gvr_annots PATH', 'GVR annotations data' ) do |gvr_annots|
        options[:gvr_annots] = gvr_annots
    end

    options[:hyper_values] = nil
    opts.on( '-v', '--hyper_values PATH', 'hypergeometric values file, with associations between GVR and phenotypes' ) do |hyper_values|
        options[:hyper_values] = hyper_values
    end

    options[:h_threshold] = 2
    opts.on( '-d', '--hyper_threshold INTEGER', 'hypergeometric values threshold of the associations' ) do |h_threshold|
        options[:h_threshold] = h_threshold.to_i
    end

    options[:patients_threshold] = 10
    opts.on( '-p', '--npatients_threshold INTEGER', 'patients number threshold used for relationships' ) do |patients_threshold|
        options[:patients_threshold] = patients_threshold.to_i
    end

    options[:gff_path] = nil
    opts.on( '-f', '--gff_file PATH', 'pathway indicating the gff file used for finding genes inside the GVR' ) do |gff_path|
        options[:gff_path] = gff_path
    end

    options[:enriched_genes] = nil
    opts.on( '-e', '--enriched_genes PATH', 'pathway indicating the enriched genes file path' ) do |enriched_genes|
        options[:enriched_genes] = enriched_genes
    end

    options[:pregene_set] = nil
    opts.on( '-t', '--pre_geneset_file PATH', 'pathway indicating the output filename with genes for gen set analysis' ) do |pregene_set|
        options[:pregene_set] = pregene_set
    end

    options[:gvr_complement] = nil
    opts.on( '-m', '--gvr_complement PATH', 'pathway indicating the output filename with gvr complement data for gene set' ) do |gvr_complement|
        options[:gvr_complement] = gvr_complement
    end
end.parse!

#Getting GVR and their statistically significant phenotype traits 
gvr_annots_and_hyper_values = load_gvr_and_hyper_values(options[:gvr_annots], options[:hyper_values], options[:h_threshold])

#Filter GVR with the desired thresholds (Hipergeometric value of x or more and y or more patients,
#defined and taken from the bash script)
filtered_GVRs = gvr_annots_and_hyper_values.select { |gvr, fields| (fields["phens"].length > 0) && (fields["n_pats"].to_i >= options[:patients_threshold]) }

#Loading GFF file with genes and annotations
genes_gff3_data = load_gff_file(options[:gff_path])

#Filtering GFF data by chromosome
filter_chrs = []
filtered_GVRs.each do |gvr, fields|
    filter_chrs.push(fields["chr"])
end
filter_chrs = filter_chrs.uniq
genes_data_filtered = genes_gff3_data.select{|id, fields| filter_chrs.include?(fields["chr"])}

#Setting additional information from the previous enrichment process to genes inside each GVR
partial_enriched_genes = load_enriched_genes(options[:enriched_genes])
full_enriched_genes_data_filtered = enrich_genes(genes_data_filtered, partial_enriched_genes)
gvrs_with_enriched_genes = fusion_gvr_and_genes(filtered_GVRs, full_enriched_genes_data_filtered)

#Writting data to a file

#For testing purposes
gvr_total_counter = 0
gvr_with_no_genes_counter = 0
puts("-"*70)
#For testing purposes end

File.open(options[:pregene_set], "w") do |file|
    gvrs_with_enriched_genes.each do |id, fields|
        gvr_total_counter += 1 #For testing purposes
        if (fields["n_pats"]).to_i > (options[:patients_threshold]).to_i
            if fields["genes"].length > 0
                fields["genes"].each do |gen_name, content|
                    file.puts("#{id}\t#{content["entrezgene_id"]}\t#{gen_name}")
                end
            else
                gvr_with_no_genes_counter += 1 #For testing purposes
                puts("GVR with id #{id} was discarded because no gene was found inside his coordinates")
            end
        end
    end
end

#For testing purposes
puts("-"*70)
print("From a total of #{gvr_total_counter} GVRs, #{gvr_with_no_genes_counter} were discarded because no genes were found inside his coordinates, and finally from the remaining #{gvr_total_counter - gvr_with_no_genes_counter} GVRs, ")
#For testing purposes

File.open(options[:gvr_complement], "w") do |file|
    gvrs_with_enriched_genes.each do |id, fields|
        file.puts("#{id}\t#{fields["chr"]}\t#{fields["n_pats"]}\t#{fields["phens"]}\t#{fields["start"]}\t#{fields["stop"]}\n")
    end
end
