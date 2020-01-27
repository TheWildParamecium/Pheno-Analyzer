#! /usr/bin/env ruby

#LOADING LIBRARIES
require 'optparse'

#FUNCTION DEFINITIONS
#For loading a file(in a suitable format)
def load_genotable_file(file_path, initial_line = 0)
    data = {}
    headers = []
    count = 0
    File.open(file_path).each do |line|
        if count >= initial_line
            fields = line.chomp.split("\t")
            pat_id = fields.shift
            chr = fields.shift
            start = fields.shift
            stop = fields.shift
                if data[pat_id].nil?
                    data[pat_id] = [[chr, start, stop]]
                else
                    data[pat_id] << [chr, start, stop]
                end
        else
            headers << line.chomp.split("\t")
        end
        count += 1
    end
    return data, headers
end

def load_phenotable_file(file_path, initial_line = 0)
    data = {}
    headers = []
    count = 0
    File.open(file_path).each do |line|
        fields = line.chomp.split("\t")
        pat_id = fields.shift
        if count >= initial_line
            data[pat_id] = fields
        else
            headers << fields
        end
        count += 1
    end
    return data, headers
end

#Transform an each column phenotype structure to an array of arrays structure 
def transform_traits_to_array(inputfile, headers)
    traits_inputs, hpo_inputs = headers
    traits_outputs = []
    hpos_outputs = []
    inputfile.each_with_index do |value, i|
        if value.downcase == "x"
            traits_outputs << traits_inputs[i] 
            hpos_outputs << hpo_inputs[i]
        end
    end
    return traits_outputs, hpos_outputs
end


#MAIN
#Defining the script parser options 
options = {}
OptionParser.new do |opts|
    options[:genotable_file] = nil
    opts.on( '-g', '--genotable_file STRING', 'Load genotable file to the script' ) do |genotable_name|
        options[:genotable_file] = genotable_name
    end

    options[:gheaders] = 1
    opts.on( '-G', '--genotable_header INTEGER', 'Number of header rows in genotable' ) do |g_headers|
        options[:gheaders] = g_headers.to_i
    end

    options[:phenotable_file] = nil
    opts.on( '-p', '--phenotable_file STRING', 'Load phenotable file to the script' ) do |phenotable_name|
        options[:phenotable_file] = phenotable_name
    end

    options[:pheaders] = 1
    opts.on( '-P', '--phenotable_header INTEGER', 'Number of header rows in phenotable' ) do |p_headers|
        options[:pheaders] = p_headers.to_i
    end

    options[:output_file] = nil
    opts.on( '-o', '--output_file STRING', 'Filename for output table' ) do |output_file|
        options[:output_file] = output_file
    end

    options[:sup_data] = nil
    opts.on( '-s', '--suplementary_data STRING', 'Filename for other stats' ) do |sup_data|
        options[:sup_data] = sup_data
    end
end.parse!

#Loading genotype and phenotype tables
genotable, genoheaders = load_genotable_file(options[:genotable_file], options[:gheaders])
phenotable, phenoheaders = load_phenotable_file(options[:phenotable_file], options[:pheaders])

#Formatting phenotype values
formatted_phenotable = {}
phenotable.each do |pat_id, pat_attrs|
    formatted_phenotable[pat_id] = transform_traits_to_array(pat_attrs, phenoheaders)
end

#Combining both tables
combination = []
genotable.each do |pat_id, geno_attrs|
    hp_data = formatted_phenotable[pat_id] 
    if !hp_data.nil?
        hp_codes = hp_data[1]  #Choosing only hp codes and leaving hp terms behind
        combination << [pat_id].concat(geno_attrs).concat([hp_codes])
    end
end

puts("-"*70)
puts("  Cantidad total de pacientes de la tabla de fenotipos: #{phenotable.length}")
puts("  Cantidad total de pacientes de la tabla de genotipos: #{genotable.length}")
puts("  Cantidad de pacientes que tienen datos de genotipado pero no de fenotipado: #{genotable.length - combination.length}")
puts("  Cantidad de pacientes que tienen datos de fenotipado pero no de genotipado: #{phenotable.length - combination.length}")
puts("  Cantidad de pacientes que tienen datos tanto de fenotipado como de genotipado: #{combination.length}")
puts("-"*70)



#Writing the standard table on a file (tsv format)
File.open(options[:output_file], "w") do |file|
    file.puts(["patID", "chr", "start", "stop", "hpcodes"].join("\t")) #PatID, Chr, start, stop, hpCodes 
    combination.each do |row|
        id = row[0].to_s
        hp_codes = row[-1].kind_of?(Array) ? row[-1].join(",") : row[-1]
        #################RELLENO NUEVO
        row[1...-1].each do |chrom_record|
            chrom_record = chrom_record.join("\t")
            file.puts(id + "\t" + chrom_record + "\t" + hp_codes)
        #####################RELLENO NUEVP
        end
    end
end

File.open(options[:sup_data], "w") do |file|
    file.puts("  Cantidad total de pacientes de la tabla de fenotipos: #{phenotable.length}")
    file.puts("  Cantidad total de pacientes de la tabla de genotipos: #{genotable.length}")
    file.puts("  Cantidad de pacientes que tienen datos de genotipado pero no de fenotipado: #{genotable.length - combination.length}")
    file.puts("  Cantidad de pacientes que tienen datos de fenotipado pero no de genotipado: #{phenotable.length - combination.length}")
    file.puts("  Cantidad de pacientes que tienen datos tanto de fenotipado como de genotipado: #{combination.length}")
end