#! /usr/bin/env ruby

#LOADING REQUIRED LIBRARIES
require 'optparse'
require "erb"
require 'json'

#FUNCTION DEFINITIONS

#MERGING FUNCTIONAL ANALYSIS RESULTS WITH OTHER GVR DATA
def get_common_get_set_results(go_bp, go_mf, kegg, reactome)
    results = {}
    keys = []
    keys = keys.concat(go_bp.keys)
    keys = keys.concat(go_mf.keys)
    keys = keys.concat(kegg.keys)
    keys = keys.concat(reactome.keys)

    keys.each do |gvr|
        results[gvr] = { "descr_go_mf" => ["No results were found"], 
            "descr_go_bp" => ["No results were found"],
            "descr_kegg" => ["No results were found"], 
            "descr_reactome" => ["No results were found"] }
    end

    keys.each do |gvr|
        if go_bp.keys.include?(gvr)
            results[gvr] = results[gvr].merge(go_bp[gvr])
        end
        if go_mf.keys.include?(gvr)
            results[gvr] = results[gvr].merge(go_mf[gvr])
        end
        if kegg.keys.include?(gvr)
            results[gvr] = results[gvr].merge(kegg[gvr])
        end
        if reactome.keys.include?(gvr)
            results[gvr]= results[gvr].merge(reactome[gvr])
        end
    end
    return results
end

def fusion_hashmaps(primary, secondary)
    results = primary.dup
    iter_keys = primary.keys
    iter_keys.each do |key|
        secondary[key].each do |field, value|
            if !(field.empty?) && !(value.empty?)
                results[key][field] = value
            end
        end
    end
    return results
end

#LOADING FILES
def load_patients_sup_data(sup_path)
    quotes = []
    File.open(sup_path).each do |line|
        quotes.push(line.chomp)
    end
    return quotes
end

def load_enriched_complement_gvr_for_geneset(gvr_compl_path)
    gvr_compl = {}
    File.open(gvr_compl_path).each do |line|
        fields = line.chomp.split("\t")
        gvr_id = fields.shift
        gvr_compl[gvr_id] = {}
        gvr_compl[gvr_id]["chr"] = fields.shift
        gvr_compl[gvr_id]["n_pats"] = fields.shift

        hpcodes = fields.shift.split('[')[1]
        hpcodes = hpcodes.split("]")[0]
        hpcodes = hpcodes.split(",").map{|dataa| dataa.strip.gsub!('"', '')}
        gvr_compl[gvr_id]["hpcodes"] = hpcodes

        gvr_compl[gvr_id]["start"] = fields.shift
        gvr_compl[gvr_id]["stop"] = fields.shift
    end
    return gvr_compl
end

def load_geneset_results(geneset_path, ont)
    gvrs = {}
    gvr_id = ""
    descr = []
    File.open(geneset_path).each do |line|
        if (line.length > 2)         
            if line.start_with?("chr")
                fields = line.chomp.split("\t")
                gvr_id = fields.shift
                gvrs[gvr_id] = {"descr_#{ont}" => []}

                descr = fields.shift
                if descr.start_with?("c(")
                    descr = descr.split("c(")[1]
                    if descr.end_with?(")")
                        descr = descr.split(")")[0]
                    end
                    descr = descr.split(",").map{|dataa| dataa.gsub!('"', '')}
                else
                    descr = [descr]
                end
            elsif line.end_with?(")")
                fields = line.chomp.split("\t")
                descr = fields.shift.split(")")[0]
                descr = (descr.split(",")).map{|dataa| dataa.gsub!('"', '')}
            else
                fields = line.chomp.split("\t")
                descr = fields.shift.split(",")
                descr = descr.map{|dataa| dataa.gsub!('"', '')}
            end
            gvrs[gvr_id]["descr_#{ont}"] = gvrs[gvr_id]["descr_#{ont}"].concat(descr)
            descr = []
            
        end
    end
    return gvrs
end

def load_file(file_path, initial_line = 0)
    data = {}
    headers = []
    count = 0
    File.open(file_path).each do |line|
        if count >= initial_line
            fields = line.chomp.split("\t")
            pat_id = fields.shift
            if data[pat_id].nil?
                data[pat_id] = {}
                data[pat_id]["genetics"] = [fields[0..2]]
                data[pat_id]["hpcodes"] = fields[3].nil? ? nil : fields[3].split(",").map{|term| term.strip} 
            else
                data[pat_id]["genetics"] << fields[0..2]
            end
        else
            headers << line.chomp.split("\t")
        end
        count += 1
    end
    return data, headers
end

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


#Support function of get_statistics, to get frequencies
def get_frequencies(data, get_total_quantity = false, total_pats = 0)
    total_occs = 0
    data.each do |key,value| 
        total_occs += value
    end

    total_occs = total_occs.to_f 
    if total_pats == 0
        total_pats = total_occs
    end

    frequencies = {}
    data.each do |key,value| 
        frequencies[key] = ((value/total_pats)*100).round(1)
    end

    if get_total_quantity
        return total_occs, frequencies
    else
        return frequencies
    end
end


#Support function of get_statistics, to get the standard deviation of a sample
def get_std_deviation(data)
    sum = 0
    total = data.length.to_f
    data.each do |num|
        if(!num.nil?)
            sum += num
        end
    end 
    mean = sum / total
    squared_sub = 0
    data.each do |num|
        if(!num.nil?)
            squared_sub += (num - mean)**2
        end
    end
    result = ((squared_sub / total)**(0.5)).round(2)
    return result
end

#Counting the total number of phenotypes, his frequencies, number of 
#chromosomes affected and which chromosomes are affected
def get_statistics(pat_dataset)
    
    n_patients = pat_dataset.keys.length.to_f
    hpcodes_occurrences = Hash.new(0)
    chromosomes_affected = Hash.new(0)
    n_phens_per_pat = Hash.new(0)
    n_chroms = Hash.new(0) 

    pat_dataset.each do |patient_id, patient_data| 
        n_chroms[pat_dataset[patient_id]["genetics"].length.to_s] += 1 

        pat_dataset[patient_id]["genetics"].each do |chromosome_data| 
            chromosomes_affected[chromosome_data[0]] += 1
        end

        if pat_dataset[patient_id]["hpcodes"].nil?
            #Patients without any phenotype described in the dataset
        else
            n_phens_per_pat[pat_dataset[patient_id]["hpcodes"].length.to_s] += 1 
            pat_dataset[patient_id]["hpcodes"].each do |hpcode|
            hpcodes_occurrences[hpcode] += 1
            end
        end
    end

    std_dev = get_std_deviation(pat_dataset.map{|id, data| pat_dataset[id]["genetics"].length})
    total_hpcodes, hpcodes_frequencies = get_frequencies(hpcodes_occurrences, get_total_quantity = true, total = n_patients)
    chromosomes_frequencies = get_frequencies(chromosomes_affected, get_total_quantity = false, total = n_patients)
    n_chroms_frequencies = get_frequencies(n_chroms)
    n_phens_per_pat_frequencies = get_frequencies(n_phens_per_pat)

    return n_patients, total_hpcodes, (total_hpcodes/n_patients).round(2), hpcodes_frequencies, chromosomes_frequencies, n_chroms_frequencies, n_phens_per_pat_frequencies, std_dev
end

#Support function of prepare_cnv_data for getting the longest CNV
def get_max_and_min_cnv(patient_dataset)
    maximum_cnv = 0
    minimum_cnv = 10000000
    patient_dataset.each do |id, fields|
        cnv_length = ( (fields["genetics"][0][2]).to_i - (fields["genetics"][0][1]).to_i )
        if  cnv_length > maximum_cnv
            maximum_cnv = cnv_length
        end
        if cnv_length < minimum_cnv
            minimum_cnv = cnv_length
        end
    end
    return maximum_cnv, minimum_cnv
end 

#Function for preparing data to be show in a histogram of distribution of CNV length
def prepare_cnv_length_data(patient_dataset, nbins, binmin = 0, binmax = 0)

    if (binmin == 0) || (binmax == 0)
        cnv_max, cnv_min = get_max_and_min_cnv(patient_dataset)
        bins_maker = ( (cnv_max.to_f - cnv_min.to_f)/ nbins).round()
        bins_array = []
        actual_bin = cnv_min
        (nbins + 1).times do 
            bins_array.push(actual_bin)
            actual_bin += (bins_maker)
        end
    else
        cnv_min = binmin
        cnv_max = binmax
        bins_maker = ((cnv_max.to_f - cnv_min.to_f) / nbins).round()
        bins_array = []
        actual_bin = cnv_min
        (nbins + 1).times do 
            bins_array.push(actual_bin)
            actual_bin += (bins_maker)
        end
    end

    whole_cnv_length_data = Hash.new(0)
    chrx_cnv_length_data = Hash.new(0)
    bin_ranges = []
    for bin in (0...(bins_array.length - 1))
        bin_ranges.push("#{(bins_array[bin])/1000}-#{(bins_array[bin+1])/1000} Kbps")
    end

    patient_dataset.each do |id, fields|

        cnv_size = ( (fields["genetics"][0][2]).to_i - (fields["genetics"][0][1]).to_i )
        

        for index in (0..(nbins - 1))    
            if (bins_array[index] <= cnv_size) && (bins_array[(index + 1)] > cnv_size)
                whole_cnv_length_data["#{(bins_array[index])/1000}-#{(bins_array[(index + 1)])/1000} Kbps"] += 1
                if fields["genetics"][0][0] == "chrX"
                    chrx_cnv_length_data["#{(bins_array[index])/1000}-#{(bins_array[(index + 1)])/1000} Kbps"] += 1
                end
                break
            end
        end
    end

    whole_cnv_array = []
    chrx_cnv_array = []
    bin_ranges.each do |bin|
        whole_cnv_array.push(whole_cnv_length_data[bin])
        chrx_cnv_array.push(chrx_cnv_length_data[bin])
    end

    return [whole_cnv_array, chrx_cnv_array, bin_ranges] 
end

#For loading the hp.obo ontology file (Human Phenotype Ontology file in flat format).
#It will be needed to analyse phenotype descriptions of the patient's dataset
def get_hpo_data(hpofile="hp.obo")
    hpo_data = Hash.new { |hash, key| hash[key] = {'is_obsolete' => false,  'parents' => [], "hpterm" => "HP term without description" } }
    hpo_id = ""
    File.open(hpofile).each do |line|
        if line.start_with? "id:" 
            hpo_id = line[4..13].chomp.dup  
        elsif line.start_with? "is_a:" #Condition for saving the parent(s) term(s) id(s) of that hp term
            hpo_data[hpo_id]["parents"].push(line[6..15].chomp)
        elsif line.start_with? "name:" 
            hpo_data[hpo_id]["hpterm"] = line[6...line.length].chomp
        end
    end
    return hpo_data
end

#Saving new HPO codes for those who has obsolete codes
def get_hpo_obsolete_data(hpofile="hp.obo")
    hpo_obsolete = {}
    hpo_id = ""
    File.open(hpofile).each do |line|
        if line.start_with? "id:" 
            hpo_id = line[4..13].chomp.dup 
        elsif line.start_with? "replaced_by:" 
            hpo_obsolete[hpo_id] = (line[13..22].chomp)
        elsif line.start_with? "alt_id" 
            hpo_obsolete[line[(8...line.length)].chomp] = hpo_id
        end 
    end
    return hpo_obsolete
end

#It replaces obsolete instances of HPO codes used in the phenotypic description of patients
def replace_obsolete(hpo_data, hpo_obsolete)
    new_hpo_data = hpo_data.dup
    hpo_obsolete.each do |old_term, new_term|
        new_hpo_data[old_term] = (new_hpo_data[new_term]).merge({"is_obsolete" => true})
    end
    return new_hpo_data
end

#For setting HPO-related data to patient's phenotypic descriptions in the dataset
def get_patients_hpo_data(patients_dataset, hpo_data)
    patients_hpo_data = {} 
    patients_dataset.each do |pat_id, pat_info|
        if !(patients_dataset[pat_id]["hpcodes"].nil?)
            patients_dataset[pat_id]["hpcodes"].each do |hpcode| 
                if !(patients_hpo_data.keys.include? hpcode) #If that term is not already in the hash
                    patients_hpo_data[hpcode] = hpo_data[hpcode] 
                end
            end
        end
    end
    return patients_hpo_data
end

#Getting the number of obsolete HPO terms in our dataset
def get_n_obsolete(hpo_patient_data)
    n_obsolete = 0
    hpo_patient_data.each do |key, value|
        if value["is_obsolete"]
            n_obsolete += 1
        end
    end
    return n_obsolete
end
    
#Improved version of find_hpo_distance_to_root (now deleted), for getting the number of steps required 
#to reach the root node/term
def find_hpo_distance_to_root_improved(hpo_patient_data, hpo_data)
    hpo_distance = {}
    hpo_patient_data.each do |key, value|
        hpo_distance[key] = breadth_first(value, hpo_data)
    end
    hpo_distance.reject!{|key, value| key.nil?}
    return hpo_distance
end

#Implementation of breadth first search, a graph-related algorithm for getting the shortest path between 
#two nodes, in this case the HPO term and the root term of the ontology. 
#It is a helper function for find_hpo_distance_to_root_improved
def breadth_first(hpo_patient_data_single, hpo_data)
    search_queue = []
    steps = {}
    search_queue = search_queue.concat(hpo_patient_data_single['parents'])
    #Setting up initial steps 
    if !search_queue.empty?   
        search_queue.each do |term|
            steps[term] = 1
        end
    end
    #Here goes the graph algorithm
    while !search_queue.empty?
        parent = search_queue.shift()
        if hpo_data[parent]['parents'].empty?
            return steps[parent]
        else
            hpo_data[parent]['parents'].each do |grandparent|
                steps[grandparent] = steps[parent] + 1  #Setting up the steps of parents of the current element in the queue
            end
            search_queue= search_queue.concat(hpo_data[parent]['parents'])
        end
    end
end


#Support function for get_hpo_comparison_data
def count_number_of_steps(hp_hash)
    steps_hash = {}
    steps_hash.default = 0
    hp_hash.each do |hpkey, hpvalue|
        steps_hash[hpvalue.to_s] += 1
    end
    return steps_hash
end

#Function that prepares HPO and patients levels (distance from the root of a HPO term) data
#to be displayed in canvasxpress graph (Logarithm of terms at a given level and density distribution
#along the levels). 
def get_hpo_comparison_data(dataset_hpo, ontology_hpo)
    dataset_steps = count_number_of_steps(dataset_hpo)
    ontology_steps = count_number_of_steps(ontology_hpo).reject{|key, value| key.empty?}

    dataset_total = 0
    ontology_total = 0
    ontology_steps.each do |key, value|
        dataset_total += dataset_steps[key]
        ontology_total += value
    end
    #Creating the results array to send to canvaxpress graphics function
    results = Array.new(5)
    results[0] = []
    results[1] = []
    results[2] = []
    results[3] = []
    results[4] = []
    
    headers = []
    #Getting the variables data
    ontology_steps.sort_by{|k, v| k.to_i}.to_h.each do |key, value|
        headers.push(key)
        if dataset_steps[key] != 0 
            results[0].push(Math.log(dataset_steps[key])) #Log number of occurrencies in the dataset at that step
        else
            results[0].push(0)
        end
        results[1].push(Math.log(ontology_steps[key]))  #Log Number of occurrencies in the whole ontology at that step

        results[2].push(((dataset_steps[key] / dataset_total.to_f)*100).to_i())     #Percentage of occurrencies at that step divided by total dataset occurrencies
        results[3].push(((ontology_steps[key] / ontology_total.to_f)*100).to_i())   #Percentage of occurrencies at that step divided by total ontology occurrencies
    end

    #Getting the mean and standard deviation values of steps of dataset and ontology
    dataset_product = 0
    ontology_product = 0
    headers.each do |step_number|
        dataset_product += dataset_steps[step_number].to_f * step_number.to_f
        ontology_product += ontology_steps[step_number].to_f * step_number.to_f
    end
    dataset_mean = (dataset_product / dataset_total.to_f).round(2)
    ontology_mean = (ontology_product / ontology_total.to_f).round(2)
    means = [dataset_mean, ontology_mean]

    #Getting the standard deviation of steps of the dataset and ontology
    dataset_squared_sub = 0
    ontology_squared_sub = 0
    headers.each do |step_number|
        dataset_squared_sub += ((step_number.to_f - dataset_mean)**2) * dataset_steps[step_number].to_f
        ontology_squared_sub += ((step_number.to_f - ontology_mean)**2) * ontology_steps[step_number].to_f
    end
    dataset_std = ((dataset_squared_sub.to_f / dataset_total.to_f)**(0.5)).round(2)
    ontology_std = ((ontology_squared_sub.to_f / ontology_total.to_f)**(0.5)).round(2)
    stds = [dataset_std, ontology_std]
    return [results, headers, means, stds]    
end

#For getting GVR(Genomic Variant Regions) from patients CNV(Copy Number Variation) data, filtered by chromosome 
def prepare_data_for_gvr(pat_dataset, threshold, gvrspecific = false)

    gvrs = {}
    filtered_by_chr = {}
    ((1..22).to_a.concat(["X"])).each do |n|
        pat_dataset.each do |key , fields|
            chr_counter = 0 
            fields["genetics"].each do |chr| 
                if chr[0] == "chr#{n}"
                    chr_counter += 1
                end
            end
            if chr_counter > 0
                filtered_by_chr[key] = fields
            end
        end
        gvrs["chr#{n}"] = get_gvr(filtered_by_chr, n, threshold, gvrspecific)

        filtered_by_chr = {}
    end

    return gvrs
end

#Suport function for prepare_data_for_gvr
def get_gvr (same_chrom_data, chr_number, threshold, gvrspecific)
    
    cvn_counter = 0
    name_counter = 0
    hits = []
    current_ids = []
    hpos = []
    starts = []
    stops = []
    same_chrom_data.each do |k,fields| 
        fields["genetics"].each do |chrs|
            starts.push(chrs[1].to_i) if (chrs[0] == "chr#{chr_number}")
            stops.push(chrs[2].to_i) if (chrs[0] == "chr#{chr_number}")
        end
    end 
    
    gvrs = Hash.new { |hash, key| hash[key] = {'start' => 0, 'stop' => 0, 'patients_id' => [], 'cvn_number' => 0, 'hpos' => "Not defined" } }

    previous_gvr = ""
    
    until stops.empty?
 
        if starts.empty? && !stops.empty?
            stop_minimum = stops.min()
            
            same_chrom_data.each do |key,fields| 
                fields["genetics"].each do |chr| 
                    if (chr[2] == stop_minimum.to_s() && chr[0] == "chr#{chr_number}")
                        hits.push(key)
                    end
                end
            end

            cvn_counter -= hits.length
            hits.each do |keyname|
                current_ids = current_ids.select{|item| item != keyname}
            end

            if true
                name_counter += 1
  
                current_ids.each do |keyname|
                    if !same_chrom_data[keyname]["hpcodes"].nil?
                        same_chrom_data[keyname]["hpcodes"].each do |hp| 
                            hpos.push(hp) if !hpos.include?(hp)
                        end
                    end
                end
                     
                gvrs["gvr#{name_counter}"]["start"] = stop_minimum
                gvrs["gvr#{name_counter}"]["cvn_number"] = cvn_counter
                gvrs["gvr#{name_counter}"]["hpos"] = hpos 
                gvrs["gvr#{name_counter}"]["patients_id"] = current_ids
                gvrs[previous_gvr]["stop"] = stop_minimum
                previous_gvr = "gvr#{name_counter}"
            end
            hpos = []
            hits = [] 
            stops.delete(stop_minimum)
            next
        end

        start_minimum = starts.min()
        stop_minimum = stops.min()

        if (start_minimum < stop_minimum)
            
            same_chrom_data.each do |key,fields| 
                fields["genetics"].each do |chr| 
                    if (chr[1] == start_minimum.to_s() && chr[0] == "chr#{chr_number}")
                        hits.push(key) 
                    end
                end
            end

            current_ids = current_ids.concat(hits)
            cvn_counter += hits.length

            if true 
                name_counter += 1
                current_ids.each do |keyname|
                    if !same_chrom_data[keyname]["hpcodes"].nil?
                        same_chrom_data[keyname]["hpcodes"].each do |hp| 
                            hpos.push(hp) if !hpos.include?(hp)
                        end
                    end
                end

                gvrs["gvr#{name_counter}"]["start"] = start_minimum
                gvrs["gvr#{name_counter}"]["cvn_number"] = cvn_counter
                gvrs["gvr#{name_counter}"]["patients_id"] = current_ids
                gvrs["gvr#{name_counter}"]["hpos"] = hpos

                if previous_gvr != ""
                    gvrs[previous_gvr]["stop"] = start_minimum
                end
                previous_gvr = "gvr#{name_counter}"    
            end
            hpos = []
            hits = [] 
            starts.delete(start_minimum)


        elsif (start_minimum >= stop_minimum)

            same_chrom_data.each do |key,fields| 
                fields["genetics"].each do |chr| 
                    if (chr[2] == stop_minimum.to_s() && chr[0] == "chr#{chr_number}")
                        hits.push(key) 
                    end
                end
            end

            cvn_counter -= hits.length
            hits.each do |keyname|
                current_ids = current_ids.select{|item| item != keyname}
            end

            if true
                name_counter += 1
  
                current_ids.each do |keyname|
                    if !same_chrom_data[keyname]["hpcodes"].nil?
                        same_chrom_data[keyname]["hpcodes"].each do |hp| 
                            hpos.push(hp) if !hpos.include?(hp)
                        end
                    end
                end
                     
                gvrs["gvr#{name_counter}"]["start"] = stop_minimum
                gvrs["gvr#{name_counter}"]["cvn_number"] = cvn_counter
                gvrs["gvr#{name_counter}"]["patients_id"] = current_ids
                gvrs["gvr#{name_counter}"]["hpos"] = hpos 

                gvrs[previous_gvr]["stop"] = stop_minimum
                previous_gvr = "gvr#{name_counter}"
            end
            hpos = []
            hits = [] 
            stops.delete(stop_minimum)
        end
    end
    if gvrspecific
        gvrs.delete_if{ |key, fields| (fields["cvn_number"] > 1) && (fields["cvn_number"] == 0)}
        return gvrs
    end
    gvrs.delete_if{ |key, fields| fields["cvn_number"] < threshold }
    return gvrs
end

#It gets the number of mutations per chromosome
def get_n_mutations_per_chromosome(pat_dataset)
    chr_muts = {}
    ((1..22).to_a.concat(["X"])).each do |n|
        muts_counter = 0
        pat_dataset.each do |key , fields|

            fields["genetics"].each do |chr| 
                if chr[0] == "chr#{n}"
                    muts_counter += 1
                end
            end
        end
        chr_muts["chr#{n}"] = muts_counter
    end
    return chr_muts
end


#MAIN
#Defining the script parser options (parameters that will be received from bash script) 
options = {}
OptionParser.new do |opts|
    options[:standard_table] = nil
    opts.on( '-s', '--standard_table PATH', 'Load the standard table to the script' ) do |standard_table|
        options[:standard_table] = standard_table
    end

    options[:table_headers] = 1
    opts.on( '-S', '--table_headers INTEGER', 'Number of header rows in standard table' ) do |table_headers|
        options[:table_headers] = table_headers.to_i
    end

    options[:hp_file] = nil
    opts.on( '-H', '--hp_file PATH', 'Load phenotable file to the script' ) do |hp_file|
        options[:hp_file] = hp_file
    end

    options[:output_file] = "results.json"
    opts.on( '-o', '--output_file PATH', 'Filename for output table' ) do |output_file|
        options[:output_file] = output_file
    end

    options[:html_file] = "report.html"
    opts.on( '-R', '--report_html PATH', 'Name of the html report file' ) do |html_file|
        options[:html_file] = html_file
    end

    options[:template] = nil
    opts.on( '-t', '--template_file PATH', 'Template to use in renderize the html report' ) do |item|
        options[:template] = item
    end

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

    options[:gvr_complement] = nil
    opts.on( '-m', '--gvr_complement PATH', 'pathway indicating the output filename with gvr complement data for gene set' ) do |gvr_complement|
        options[:gvr_complement] = gvr_complement
    end

    options[:gene_set_results_go_bp] = nil
    opts.on( '-z', '--geneset_results_go_bp PATH', 'pathway indicating the result of over-representational test of GO BP' ) do |gene_set_results_go_bp|
        options[:gene_set_results_go_bp] = gene_set_results_go_bp
    end

    options[:gene_set_results_go_mf] = nil
    opts.on( '-Z', '--geneset_results_go_mf PATH', 'pathway indicating the result of over-representational test of GO BP' ) do |gene_set_results_go_mf|
        options[:gene_set_results_go_mf] = gene_set_results_go_mf
    end

    options[:gene_set_results_kegg] = nil
    opts.on( '-x', '--geneset_results_kegg PATH', 'pathway indicating the result of over-representational test of KEGG' ) do |gene_set_results_kegg|
        options[:gene_set_results_kegg] = gene_set_results_kegg
    end

    options[:gene_set_results_reactome] = nil
    opts.on( '-X', '--geneset_results_reactome PATH', 'pathway indicating the result of over-representational test of Reactome' ) do |gene_set_results_reactome|
        options[:gene_set_results_reactome] = gene_set_results_reactome
    end

    options[:sup_data] = nil
    opts.on( '-w', '--suplementary_data STRING', 'Filename for other stats' ) do |sup_data|
        options[:sup_data] = sup_data
    end
end.parse!

#Loading standard table and hpo ontology
standard_table, tableheaders = load_file(options[:standard_table], options[:table_headers])
hpo_data = get_hpo_data(options[:hp_file])
hpo_obsolete_data = get_hpo_obsolete_data(options[:hp_file])
hpo_data = replace_obsolete(hpo_data, hpo_obsolete_data)
#Loading and getting patients HPO data
patients_hpo_data = get_patients_hpo_data(standard_table, hpo_data)
patients_supl_data = load_patients_sup_data(options[:sup_data])
#Getting the number of obsolete HPO codes in the patient's dataset
n_obsolete = get_n_obsolete(patients_hpo_data)
#Getting some statistics
n_pats, n_phens, phens_per_pat, phens_freq, chroms_freq, n_chroms, n_phens_per_pat_frequencies, phens_per_pat_std_dev = get_statistics(standard_table)
#Getting patient's HPO levels (distance from the root term) and preparing the plot
hpo_distance_to_root = find_hpo_distance_to_root_improved(patients_hpo_data, hpo_data)
whole_ontology_distance_to_root = find_hpo_distance_to_root_improved(hpo_data, hpo_data)
hpo_comparison_data = get_hpo_comparison_data(hpo_distance_to_root, whole_ontology_distance_to_root)
#Getting GVR(Genomic Variant Region) from CNV(Copy Number Variation) patient's data
gvr_data_1 = prepare_data_for_gvr(standard_table, 1, gvrspecific = true)
gvr_data_2 = prepare_data_for_gvr(standard_table, 2)
gvr_data_10 = prepare_data_for_gvr(standard_table, 10)
gvr_data_x = prepare_data_for_gvr(standard_table, 1)
#Getting number of mutations per chromosome
n_muts_per_chr = get_n_mutations_per_chromosome(standard_table)
#Getting distribution of CNV length in whole dataset and ChrX affected patients
whole_cnv_length, chrx_cnv_length, bin_ranges = prepare_cnv_length_data(standard_table, 15, 1, 6000000)
#Getting GVR with statistically significant phenotype traits associated 
gvr_annots_and_hyper_values = load_gvr_and_hyper_values(options[:gvr_annots], options[:hyper_values], options[:h_threshold])
#Loading enriched GVR results from functional analysis of GO BP and MF, KEGG & Reactome
geneset_complement_gvrdata = load_enriched_complement_gvr_for_geneset(options[:gvr_complement])

geneset_results_go_bp = load_geneset_results(options[:gene_set_results_go_bp], "go_bp")
geneset_results_go_bp_complete = fusion_hashmaps(geneset_results_go_bp, geneset_complement_gvrdata)

geneset_results_go_mf = load_geneset_results(options[:gene_set_results_go_mf], "go_mf")
geneset_results_go_mf_complete = fusion_hashmaps(geneset_results_go_mf, geneset_complement_gvrdata)

geneset_results_kegg = load_geneset_results(options[:gene_set_results_kegg], "kegg")
geneset_results_kegg_complete = fusion_hashmaps(geneset_results_kegg, geneset_complement_gvrdata)

geneset_results_reactome = load_geneset_results(options[:gene_set_results_reactome], "reactome")
geneset_results_reactome_complete = fusion_hashmaps(geneset_results_reactome, geneset_complement_gvrdata)

common_gen_set_results = get_common_get_set_results(geneset_results_go_bp_complete,
     geneset_results_go_mf_complete,
     geneset_results_kegg_complete,
     geneset_results_reactome_complete)

#Generating html file and merging results with the ERB template file
template = File.open(options[:template]).read
html = ERB.new(template).result(binding)
#Writing the html report to the specified filename
File.open(options[:html_file], "w") do |file|
    file.puts(html)
end

#CODE FOR TESTING PURPOSES
puts("#{common_gen_set_results.length} got results in functional analysis")
puts("-"*70)