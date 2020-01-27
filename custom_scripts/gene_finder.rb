#! /usr/bin/env ruby

#LOADING LIBRARIES
require 'optparse'

#FUNCTION DEFINITIONS

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
                gen_annots[id]["pseudo"] = metadata.include?("pseudo=true") ? "true" : "false"
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

#For loading the hp.obo ontology file within ruby's script
def get_hpo_data(hpofile="hp.obo")
    id_ph = nil
    alt_array = []
    hpo_data = Hash.new { |hash, key| hash[key] = {'is_obsolete' => false, 'parents' => [], "hpterm" => "Patients without phenotype described" } }
    hpo_id = ""
    File.readlines(hpofile).each do |line|
        if line.start_with? "id:" #if number is an id
            if !id_ph.nil?  #Setting alternatives ids gathered to the previous hp term
                alt_array.each do |term|
                    hpo_data[term] = hpo_data[id_ph]
                end
                id_ph = nil
                alt_array = []
            end
            hpo_id = line[4..13].chomp.dup  #Get the id term of the hpo
            #By default the term will not be obsolet if it is not specified
     
        elsif line.start_with? "is_obsolete:" #If the hp term if obsolete
            hpo_data[hpo_id]["is_obsolete"] = true 
    
        elsif line.start_with? "replaced_by:" #If the hp term was obsolete
            hpo_data[hpo_id]["replaced_by"] = line[13..22].chomp #Keeps the term for which it was replaced by
        
        elsif line.start_with? "is_a:" #Condition for saving the parent(s) term(s) id(s) of that hp term
            hpo_data[hpo_id]["parents"].push(line[6..15].chomp)
        
        elsif line.start_with? "name:" 
            hpo_data[hpo_id]["hpterm"] = line[6...line.length].chomp

        elsif line.start_with? "alt_id" 
            alt_array.push(line[8...line.length].chomp)
            id_ph = hpo_id
        end 
        
    end

    return hpo_data
end

def give_genes_to_each_GVR(gvrs_data, genes_data)
    candidate_genes = {}
    genes_data.each do |gen_id, gen_fields|
        gvrs_data.each do |gvr_id, gvr_fields|
            if gen_fields["chr"] == gvr_fields["chr"]
                if (gen_fields["start"] > gvr_fields["start"]) && (gen_fields["start"] < gvr_fields["stop"])
                    candidate_genes[gen_id] = gen_fields
                elsif (gen_fields["stop"] > gvr_fields["start"]) && (gen_fields["stop"] < gvr_fields["stop"])
                    candidate_genes[gen_id] = gen_fields
                end
            end
        end
    end
    return candidate_genes
end

#MAIN
#Defining the script parser options 
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

    options[:patients_threshold] = 2
    opts.on( '-p', '--npatients_threshold INTEGER', 'patients number threshold used for relationships' ) do |patients_threshold|
        options[:patients_threshold] = patients_threshold.to_i
    end

    options[:gff_path] = nil
    opts.on( '-f', '--gff_file PATH', 'pathway indicating the gff file used for finding genes inside the GVR' ) do |gff_path|
        options[:gff_path] = gff_path
    end

    options[:candidate_genes] = nil
    opts.on( '-c', '--candidate_genes PATH', 'pathway indicating the candidate genes output file' ) do |candidate_genes|
        options[:candidate_genes] = candidate_genes
    end

    options[:hp_file] = nil
    opts.on( '-H', '--hp_file PATH', 'Load phenotable file to the script' ) do |hp_file|
        options[:hp_file] = hp_file
    end
end.parse!

#Loading HPO ontology 
hpo_data = get_hpo_data(options[:hp_file])

#Getting GVR and their statistically significant phenotype traits 
gvr_annots_and_hyper_values = load_gvr_and_hyper_values(options[:gvr_annots], options[:hyper_values], options[:h_threshold])

#Filter GVR with the desired thresholds (Hipergeometric value of 2 or more and 10 or more patients)
filtered_GVRs = gvr_annots_and_hyper_values.select { |gvr, fields| (fields["phens"].length > 0) && (fields["n_pats"].to_i >= options[:patients_threshold]) }

#Loading GFF file with genes and annotaciones
genes_data = load_gff_file(options[:gff_path])

#Filtering GFF data by chromosome
filter_chrs = []
filtered_GVRs.each do |gvr, fields|
    filter_chrs.push(fields["chr"])
end
filter_chrs = filter_chrs.uniq
genes_data_filtered = genes_data.select{|id, fields| filter_chrs.include?(fields["chr"])}

#Getting candidate genes of each GVR
candidate_genes = give_genes_to_each_GVR(filtered_GVRs, genes_data_filtered)

#Writing candidate genes on a file

counter = 0
File.open(options[:candidate_genes], "w") do |file|
    candidate_genes.each do |gen_id, fields|
        if fields["pseudo"] == "false"
            file.puts(fields["Name"])
        else
            counter += 1
        end
    end
end

puts("-"*70)
puts("Se descartaron #{counter} pseudogenes en el análisis")
puts("-"*70)

