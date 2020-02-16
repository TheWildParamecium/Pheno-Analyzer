#! /usr/bin/env ruby

#LOADING REQUIRED LIBRARIES
require 'optparse'

#FUNCTION DEFINITIONS

#Loading the standard table.
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
                data[pat_id]["hpcodes"] = fields[3].nil? ? nil : fields[3].split(",") 
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

#Getting GVR (Genomic Variant Regions) from patients CNV(Copy Number Variation) data,
#filtered by chromosome 
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

#Support function for prepare_data_for_gvr
def get_gvr (same_chrom_data, chr_number, threshold, gvrspecific)
    starts = []
    stops = []
    name_counter = 0 
    current_ids = []
    previous_gvr = ""
    gvrs = Hash.new { |hash, key| hash[key] = {'start' => 0, 'stop' => 0, 'patients_id' => [], 'cvn_number' => 0} }
    
    same_chrom_data.each do |k,fields| 
        fields["genetics"].each do |chrs|
            starts.push(chrs[1].to_i) if (chrs[0] == "chr#{chr_number}")
            stops.push(chrs[2].to_i) if (chrs[0] == "chr#{chr_number}")
        end
    end 
    
    until stops.empty?
        hpos = []
        hits = []
        name_counter += 1
        start_minimum = starts.min()
        stop_minimum = stops.min()

        #Condicion para evitar error en la comparacion cuando starts ya esta vacio
        if start_minimum.nil?
            start_minimum = stop_minimum + 1
        end

        if (start_minimum < stop_minimum)
            
            same_chrom_data.each do |key,fields| 
                fields["genetics"].each do |chr| 
                    if (chr[1] == start_minimum.to_s() && chr[0] == "chr#{chr_number}")
                        hits.push(key)
                    end
                end
            end

            current_ids = current_ids.concat(hits)

            gvrs["gvr#{name_counter}"]["start"] = start_minimum
            gvrs["gvr#{name_counter}"]["cvn_number"] = current_ids.length
            gvrs["gvr#{name_counter}"]["patients_id"] = current_ids.dup

            if previous_gvr != ""
                gvrs[previous_gvr]["stop"] = start_minimum
            end
            previous_gvr = "gvr#{name_counter}"    
        
            starts.delete(start_minimum)


        elsif (start_minimum >= stop_minimum)

            same_chrom_data.each do |key,fields| 
                fields["genetics"].each do |chr| 
                    if (chr[2] == stop_minimum.to_s() && chr[0] == "chr#{chr_number}")
                        hits.push(key) 
                    end
                end
            end

            current_ids = current_ids - hits

            gvrs["gvr#{name_counter}"]["start"] = stop_minimum
            gvrs["gvr#{name_counter}"]["cvn_number"] = current_ids.length
            gvrs["gvr#{name_counter}"]["patients_id"] = current_ids.dup

            gvrs[previous_gvr]["stop"] = stop_minimum
            previous_gvr = "gvr#{name_counter}"

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

    options[:genphen] = "Tripartite_network.txt"
    opts.on( '-T', '--tripartite_network PATH', 'Filename for genotype, phenotype and patient associations' ) do |genphen|
        options[:genphen] = genphen
    end

    options[:gvrinfo] = "GVR_annots.txt"
    opts.on( '-G', '--gvr_data PATH', 'Filename for GVR data file' ) do |gvrinfo|
        options[:gvrinfo] = gvrinfo
    end

end.parse!

#Loading standard table of patient dataset
standard_table, tableheaders = load_file(options[:standard_table], options[:table_headers])

#Getting GVR (Genomic Variant Region) from patient's CNV mutations
gvr_data = prepare_data_for_gvr(standard_table, 2)

#Writing output files
File.open(options[:genphen], "w") do |file|
    standard_table.each do |pat, fields|
        if !(fields["hpcodes"].nil?)
            fields["hpcodes"].each do |hpcode|
                file.puts(hpcode + "\t" + "ID#{pat}")
            end
        end
    end
    gvr_data.each do |chr, patients|
        gvr_data[chr].each do |gvr, fields|
            fields["patients_id"].each do |pat|
                file.puts("#{chr}.#{gvr}.#{fields["cvn_number"]}" + "\t" + "ID#{pat}")
            end
        end
    end
end

File.open(options[:gvrinfo], "w") do |file|
    gvr_data.each do |chr, patients|
        gvr_data[chr].each do |gvr, fields|
            file.puts("#{chr}.#{gvr}.#{fields["cvn_number"]}" + "\t" + fields["start"].to_s + "\t" + fields["stop"].to_s + "\t" + chr + "\t" + fields["cvn_number"].to_s)
        end
    end
end