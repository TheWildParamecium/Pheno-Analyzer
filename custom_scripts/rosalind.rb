#!/usr/bin/env ruby

#LOADING DEPENDENCIES
require 'open-uri'
require 'benchmark'


#FUNCTION DEFINITIONS

#Utilities

#Binary search function. Given a sorted array and a item, it returns the position where the item is in the array
#(if it exist), otherwise it returns nil(null). It runs at O(log n) speed.
def use_binary_search(list, item)
    low = 0
    high = list.length - 1
    while low <= high
        mid = (low + high)
        guess = list[mid]
        if guess == item
            return mid
        end
        if guess > item
            high = mid - 1
        else
            low = mid + 1
        end
    end
    return nil
end

#Algoritmic function to find the smallest item in a given array. It runs at O(n) speed as it assumes 
# the array is unsorted. It works as a helper function for use_selection_sort
def find_smallest(arr)
    smallest = arr[0] #Stores the smallest value
    smallest_index = 0 #Stores the index of the smallest value
    for i in (1...arr.length)
        if arr[i] < smallest
            smallest = arr[i]
            smallest_index = i
        end
    end
    return smallest_index
end

#Algoritmic function to sort an unsorted array. It runs at O(n x n) speed.
def use_selection_sort(arr)
    newArr = []
    for i in (0...arr.length)
        smallest = find_smallest(arr)
        puts(arr[smallest])
        newArr.push(arr.delete_at(smallest))
    end
    return newArr
end

#Algoritmic function to sort an unsorted array. It runs at O(n x n) speed in the worst case, O(n x log n) on 
#the average case.
def use_quick_sort(arr)
    if arr.length < 2
        return arr #Base case: arrays with 0 or 1 element are already “sorted.”
    else
    pivot = arr[0] #Recursive case
    less =  arr[(1...arr.length)].select{|x| x <= pivot} #Sub-array of all the elements less than the pivot
    greater = arr[(1...arr.length)].select{|x| x > pivot} #Sub-array of all the elements greater than the pivot
    return use_quick_sort(less) + [pivot] + use_quick_sort(greater)
    end
end

#Algoritmic function to be used with unweigthed graphs. It searches from "item" until it finds another item that 
#fulfills some lambda logic_function (be a certain item, be greater than, etc) and returns some of the following options:
#1)true if searched item is found, (then use returned = "found");
#2)the number of steps required, (then use returned = "steps")
def use_breadth_first(item, graph, logic_function = ->(x){graph[x].empty?} , returned = "steps" )
    search_queue = []
    steps = {}
    search_queue = search_queue.concat(graph[item])
    searched = []
    #Setting up initial steps 
    if !search_queue.empty?
        search_queue.each do |term|
            steps[term] = 1
        end
    end
    #Here goes the graph algorithm
    while !search_queue.empty?
        person = search_queue.shift()
        if !( searched.include?(person) )

            if logic_function.call(person)
                if returned == "steps"
                    return steps[person]
                end
                if returned == "found"
                    return true
                end
            else
                if !(graph[person].nil?) 
                    graph[person].each do |related|
                        steps[related] = steps[person] + 1  #Setting up the steps of parents of the current element in the queue
                    end
                    search_queue = search_queue.concat(graph[person])
                end
            end

        end
    end
    return false
end

#Implementation of the Dijkstra algorithm. It is an function to be used with weigthed graphs. It searches the fastest 
#way between two nodes or between a given node and the final node, taking into account the weigths of each edge, in
#order to find the "cheapest" way. Recall this algorithm only work with DAGs (Direct Acyclic Graphs), and only if 
#the weigths are equal/greater than 0


grafo = {}
grafo["start"] = {}
grafo["start"]["a"] = 6
grafo["start"]["b"] = 2
grafo["a"] = {}
grafo["a"]["fin"] = 1
grafo["b"] = {}
grafo["b"]["a"] = 3
grafo["b"]["fin"] = 5
grafo["fin"] = {}

    infinity = Float::INFINITY
    costs = {}
    costs["a"] = 6
    costs["b"] = 2
    costs["fin"] = infinity

    parents = {}
    parents["a"] = “start”
    parents["b"] = “start”
    parents["fin"] = None

    processed = []

    node = find_lowest_cost_node(costs)
    while node != nil
        cost = costs[node]
        neighbors = graph[node]
        neighbors.keys.each do |n|
            new_cost = cost + neighbors[n]
            if costs[n] > new_cost
                costs[n] = new_cost
                parents[n] = node
            end
        processed.append(node)
        node = find_lowest_cost_node(costs)
        end
    end
end

def find_lowest_cost_node(costs)
    lowest_cost = float(“inf”)
    lowest_cost_node = None
    for node in costs  #Go through each node.
        cost = costs[node]
        if cost < lowest_cost and node not in processed:
            lowest_cost = cost  #set it as the new lowest-cost node.
            lowest_cost_node = node
        end
    end
    return lowest_cost_node
end
#Probabilistic and mathematical functions

#Fibonacci sequence. Months represent the number of generations after we want to retrieve a results,
#and offspring the offspring produced by each pair at each generation
def fibonacci_rabbits(months, offspring, sequence = false)
    pairs = [1, 1]
    x = 2
    while x < months
        pairs.push( (pairs[x-1]) + (pairs[(x-2)]*offspring)  )
        x += 1
    end
    if sequence
        return pairs
    end
    return pairs[months - 1]
end

#Mortal rabbits Fibonacci sequence. Months represent the number of generations after we want to retrieve a results,
#and lived months the number of months that each pair of rabbits produced at a given generation live
def mortal_fibonacci_rabbits(months, lived_months)
    pairs = [1, 1]
    x = 2
    while x < months
        total = 0
        if lived_months > x
            pairs.push( (pairs[x-1]) + (pairs[(x-2)])  )
        else
            for n in (2..lived_months)
                total += pairs[x - n]
            end
            pairs.push(total)
        end
        x += 1
    end
    return pairs[months - 1]
end

#Mendel firt Law. Given k,m,n as the number of homozygous dominant, heterozygous and homozygous recessive indiviuals,
#it gives the probability that two random organisms selected will produce a dominant allele (either homozygous or 
# heterozygous)
def get_dominant_allele_probability(k, m, n)
    poblation = k + m + n
    probs = {}
    probs["mn"] = ( (m.to_f/poblation)*(n.to_f/(poblation - 1)) + (n.to_f/poblation)*(m.to_f/(poblation -1)) ) * 0.5
    probs["kn"] = ( (k.to_f/poblation)*(n.to_f/(poblation - 1)) + (n.to_f/poblation)*(k.to_f/(poblation - 1)) ) * 1
    probs["km"] = ( (k.to_f/poblation)*(m.to_f/(poblation - 1)) + (m.to_f/poblation)*(k.to_f/(poblation - 1)) ) * 1
    probs["kk"] = ( (k.to_f/poblation)*((k.to_f - 1)/(poblation - 1))) * 1
    probs["mm"] = ( (m.to_f/poblation)*((m.to_f - 1)/(poblation - 1))) * 0.75
    return (probs.values.sum()).round(5)
end

#Given the number of couples with a certain genetic background at a given allele and the number of childs for
#each couple, it returns the expected number organisms with the dominant phenotype 
def calculate_expected_dominant_offspring(xAA_AA, xAA_Aa, xAA_aa, xAa_Aa, xAa_aa, xaa_aa, offspring_each = 2)
    xAA_AA_prob = offspring_each.to_f * xAA_AA * 1 
    xAA_Aa_prob = offspring_each.to_f * xAA_Aa * 1 
    xAA_aa_prob = offspring_each.to_f * xAA_aa * 1 
    xAa_Aa_prob = offspring_each.to_f * xAa_Aa * 0.75 
    xAa_aa_prob = offspring_each.to_f * xAa_aa * 0.5 
    xaa_aa_prob = offspring_each.to_f * xaa_aa * 0 
    return xAA_AA_prob + xAA_Aa_prob + xAA_aa_prob + xAa_Aa_prob + xAA_aa_prob + xAa_aa_prob + xaa_aa_prob
end


#Given a n value it returns all the possible permutations of length n with that n numbers
def get_permutations(n)
    return (1..n).to_a.permutation.to_a
end

#Given numbers N and K, counts the total number of partial permutations, also known as variations
#of k objects that can be formed from a collection of n objects
def get_partial_permutations(n,k)
    result = 1
    while (k > 0)
        result *= n
        n -= 1
        k -= 1
    end
    return (result % 1000000)
end

#Given a number N it returns every permutations of length N of the N numbers between -N and N, without repeating
#the same number inside every permutation
def get_signed_permutations(n)
    numbers = (-n..n).to_a
    numbers.delete(0)
    perms = numbers.permutation(n).to_a
    perms = perms.select{|item| item.map{|subitem|subitem.abs()}.uniq.length == item.length}
    return perms.uniq
end

#Returns the factorial of a given number
def get_factorial(n)
    total = 1
    while n > 1
        total *= n
        n -= 1
    end
    return total
end

#Given a permutation of length N at a certain orden, it finds the longest increasing and decreasing numerical
#subsequences
def get_longest_increasing_subsequence(numeric_array)
end


def random_string(dna_string, gc_prob)
    probs = {}
    result = 1
    probs["G"] = gc_prob.to_f / 2
    probs["C"] = gc_prob.to_f / 2
    probs["A"] = (1 - gc_prob.to_f) / 2
    probs["T"] = (1 - gc_prob.to_f) / 2
    dna_string.each_char do |nt|
        result *= probs[nt]
    end
    return Math.log(result)
end

puts(random_string("ACGATACAA", 0.129))


def load_graph(numeric_file)
    n_graph = {}
    File.open(numeric_file).each do |line|
        if !(line.nil?)
            if line.chomp.split(" ").length == 2
                item = line.chomp.split(" ").map{|item| item.to_i}
                n_graph[item[0]] = item[1]
            end
        end
    end
    return n_graph
end

#Given a direct graph as a adjacency hashmap returns the minimum number of edges needed to form a tree 
def calculate_needed_edges_to_form_tree (adjacency_hashmap)
    minimum_array = []
    minimum2_array = []
    test1 = adjacency_hashmap.values
    test2 = adjacency_hashmap.keys

    test1.each do |node|
        if !(adjacency_hashmap.keys.include?(node))
            minimum_array.push(node)
        end
    end
    test2.each do |node2|
        if !(adjacency_hashmap.values.include?(node2))
            minimum2_array.push(node2)
        end
    end
    return [minimum_array.uniq.length, minimum2_array.uniq.length].min()
end
testgraph = load_graph("rosa.txt")
puts(calculate_needed_edges_to_form_tree(testgraph))
#Read a text file and store all the numbers in a numeric array
def read_numbers(numeric_file)
    numbers = []
    File.open(numeric_file).each do |line|
        if !(line.nil?)
            numbers.concat(line.chomp.split(" ").map{|item| item.to_i})
        end
    end
    return numbers
end

#Read a fasta file and returns a hash map with the sequence ID and the sequence itself
def read_fasta(fasta_file)
    fasta_name = ""
    fasta_seqs = {}
    seq = ""
    File.open(fasta_file).each do |line|
        if !(line.nil?)
            if line.start_with? ">"
                seq = ""
                fasta_name = line.chomp.split(">")[1]
            elsif
                seq = seq + line.chomp
                fasta_seqs[fasta_name] = seq
            end
        end
    end
    return fasta_seqs
end

def read_online_fasta(url)
    fasta_name = ""
    fasta_seqs = {}
    seq = ""
    open(url) do |file|
        file.read.split("\n").each do |line|
            if !(line.nil?)
                if line.start_with? ">"
                    seq = ""
                    fasta_name = line.split(">")[1]
                elsif
                    seq = seq + line
                    fasta_seqs[fasta_name] = seq
                end
            end
        end
    end
    return fasta_seqs
end

#Returns the positions of a given motif inside a DNA sequence
def get_motif (dna_seq, substring)
    pos = []
    for index in (0..(dna_seq.length - substring.length))
        if dna_seq[index...(index + substring.length)] == substring
            pos.push(index + 1)
        end
    end   
    return pos
end

#Returns the GC percent of a DNA sequence
def get_gc_percent(dna_seq)
    gc_counter = 0
    dna_seq.upcase.each_char do |nt|
        if nt == "C" || nt == "G"
            gc_counter += 1
        end
    end
    gc_perc = (((gc_counter.to_f / dna_seq.length)) * 100).round(6)
    return gc_perc
end 

#Given two DNA sequences of equal length, it returns the hamming distance, that is the number of differences
#between the two sequences, as a method to meassure the number of point mutations
def count_point_mutations(dna_seq1, dna_seq2)
    n_muts = 0
    for nt_pos in (0...(dna_seq1.length))
        if dna_seq1[nt_pos] != dna_seq2[nt_pos]
            n_muts += 1
        end
    end
    return n_muts
end

#Returns the number of each nucleotid of a given DNA sequence 
def count_nt(dna_seq)
    nts_array = dna_seq.upcase.split("")
    nts = {"A" => 0 , "G" => 0, "C" => 0, "T" => 0}
    
    nts_array.each do |nt|
        if nt == "A" 
            nts["A"] += 1
        end
        if nt == "C" 
            nts["C"] += 1
        end
        if nt == "T" 
            nts["T"] += 1
        end
        if nt == "G" 
            nts["G"] += 1
        end
    end
    return nts
end

#Returns the RNA transcript of a DNA template (but not his complement)
def get_rna_transcript(dna_seq)
    nts_array = dna_seq.upcase.split("")
    transcribed = []
    nts_array.each do |nt|
        if nt == "T" 
            transcribed.concat(["U"])
        else
            transcribed.concat([nt])
        end
    end
    transcribed = transcribed.join("")
    return transcribed
end

#Returns the reverse complement chain of a DNA sequence
def get_reverse_complement(dna_seq)
    nts_array = dna_seq.upcase.split("")
    reverse = nts_array.reverse
    complement = []
    reverse.each do |nt|
        if nt == "T"
            complement.concat(["A"])
        elsif nt == "A"
            complement.concat(["T"])
        elsif nt == "C"
            complement.concat(["G"])
        elsif nt == "G"
            complement.concat(["C"])
        end
    end
    complement = complement.join("")
    return complement            
end

#Given an array of DNA sequences it return the consensus sequence and the nucleotides statistics at each position, and the profile matrix.
#Sequences need to be the same length in order to work
def get_consensus_sequence(seqs_array)
    a_array = Array.new(seqs_array[0].length).fill(0)
    t_array = Array.new(seqs_array[0].length).fill(0)
    c_array = Array.new(seqs_array[0].length).fill(0)
    g_array = Array.new(seqs_array[0].length).fill(0)
    consensus = []
    for pos in (0...seqs_array[0].length)
        for seq in (0...seqs_array.length)
            if seqs_array[seq][pos] == "A" || seqs_array[seq][pos] == "a"
                a_array[pos] += 1
            elsif seqs_array[seq][pos] == "T" || seqs_array[seq][pos] == "t"
                t_array[pos] += 1
            elsif seqs_array[seq][pos] == "C" || seqs_array[seq][pos] == "c"
                c_array[pos] += 1
            elsif seqs_array[seq][pos] == "G" || seqs_array[seq][pos] == "g"
                g_array[pos] += 1
            end
        end
    end
    for pos in (0...a_array.length)
        consensus_nt = [a_array[pos], t_array[pos], g_array[pos], c_array[pos]].max()
        if consensus_nt == a_array[pos]
            consensus[pos] = "A"
        elsif consensus_nt == t_array[pos]
            consensus[pos] = "T"
        elsif consensus_nt == c_array[pos]
            consensus[pos] = "C"
        elsif consensus_nt == g_array[pos]
            consensus[pos] = "G"
        end
    end
    consensus = consensus.join("")
    a_array = a_array.join(" ")
    t_array = t_array.join(" ")
    g_array = g_array.join(" ")
    c_array = c_array.join(" ")
    profile_matrix = {"A" => a_array, "T" => t_array, "C" => c_array, "G" => g_array}
    return [consensus, profile_matrix]
end

#Given a list of sequences it returns the biggest shared motif in all the sequences
def get_shared_motif (seqs_array)
    template = seqs_array.min()
    seqs_array.delete(template)
    temp_len = template.length
    consensus = ""
    for stop in (0..temp_len).to_a.reverse
        start = temp_len - stop
        for start_iter in (0..start)
            if start_iter == (start_iter + stop)
                break
            end
            #To see how it works
            #print("range => #{start_iter}-#{start_iter+stop}, seq =>  ")
            #puts("#{template[start_iter...(start_iter + stop)]} ")
            if (seqs_array.all? { |seq| seq.include?(template[start_iter...(stop + start_iter)]) }) && consensus.length < ((template[start_iter...(stop + start_iter)]).length)
                #If you want to see other hits, place this return at the end of
                #the function and activate de line below
                #puts("hit with =>  #{template[start_iter...(start_iter + stop)]}")
                consensus = template[start_iter...(stop + start_iter)]
                return consensus
            end
        end
    end
end

#Given a DNA sequence it return the position and length of every palindromic sequence inside
def get_palindromic_sequences(dna_seq)
    seq_len = dna_seq.length
    min_len = 3
    max_len = 12
    pos_and_len = []
    for palin_len in (min_len..max_len)
        rango = seq_len - palin_len
        for nt in (0...rango)
            temp = dna_seq[nt..(nt + palin_len)]
            if temp == get_reverse_complement(temp)
                puts("#{nt + 1} #{palin_len + 1}")
                puts("#{temp} - #{get_reverse_complement(temp)}")
                puts()
                pos_and_len.push([nt + 1, palin_len + 1])
            end
        end
    end
    return pos_and_len
end

#Translate a mRNA sequente into a protein sequence. If stop_required is true, then if the sequence doesnt have
# a stop codon, it returns nil instead
def get_prot_seq(mrna, stop_required = false)
    prot = ""
    codons = {"GCU" => "A","GCC" => "A","GCA" => "A","GCG" => "A","GUU" => "V","GUC" => "V",
             "GUA" => "V","GUG" => "V","GAU" => "D","GAC" => "D","GAA" => "E","GAG" => "E",
             "GGU" => "G","GGC" => "G","GGA" => "G","GGG" => "G","AUU" => "I","AUC" => "I",
             "AUA" => "I","AUG" => "M","ACU" => "T","ACC" => "T","ACA" => "T","ACG" => "T",
             "AAU" => "N","AAC" => "N","AAA" => "K","AAG" => "K","AGU" => "S","AGC" => "S",
             "AGA" => "R","AGG" => "R","CUU" => "L","CUC" => "L","CUA" => "L","CUG" => "L",
             "CCU" => "P","CCC" => "P","CCA" => "P","CCG" => "P","CAU" => "H","CAC" => "H",
             "CAA" => "Q","CAG" => "Q","CGU" => "R","CGC" => "R","CGA" => "R","CGG" => "R",
             "UUU" => "F","UUC" => "F","UUA" => "L","UUG" => "L","UCU" => "S","UCC" => "S",
             "UCA" => "S","UCG" => "S","UAU" => "Y","UAC" => "Y","UAA" => "STOP",
             "UAG" => "STOP","UGU" => "C","UGC" => "C","UGA" => "STOP","UGG" => "W"}
    for i in ((0...mrna.length).step(3))
        if codons[mrna[i...(i+3)]] == "STOP"
            return prot
        elsif
            prot += codons[mrna[i...(i+3)]]
        end
    end
    if stop_required
        return nil
    end
    return prot
end

#Given a DNA sequence and a list of introns, it returns the DNA sequence with de ORF exons
def eliminate_introns(dna_seq, introns)
    template = dna_seq
    introns.each do |intron|
        template = template.split(intron)
        template = template.join("")
    end
    return template
end

#Given a DNA sequence and a DNA motif, it return the coordinates where the 
#the spliced motif can be found inside the DNA sequence, or returns Spliced motif
#not found if no full asociation could be done
def get_spliced_motif (dna_seq, spliced_motif)
    dna_seq = dna_seq.split("")
    spliced_motif = spliced_motif.split("")
    act_pos = 0
    motif_coordinates = []
    spliced_motif.each do |motif_nt|
        for seq_nt in (0...dna_seq.length)
            if (motif_nt == dna_seq[seq_nt]) && (seq_nt > act_pos)
                motif_coordinates.push(seq_nt + 1)
                act_pos = seq_nt
                break
            end
        end
    end
    
    if motif_coordinates.length == spliced_motif.length
        return motif_coordinates
    else
        return "Spliced motif not found"
    end
end




#Given an alphabet(as an array of letters) and a integer "N" obtain k-mers (permutations of every substring of length n) 
#lexicographically ordered
def get_kmers(letters, n, perm="")
    letters = letters.sort
    perms_array = []
    if n == 1
        letters.each do |letter|
            perms_array.push((perm + letter))
        end
        return perms_array
    elsif n > 1
        letters.each do |letter|
            perms_array = perms_array.concat(get_kmers(letters, n - 1, (perm + letter)))
        end
    end
    return perms_array.sort
end

#Given a genetic string and a K value (representing the substring length) it returns the number of times that
#every substring appears in the genetic string. Set kmers_key true to get also the key-value pairs
def get_kmers_composition(dna_string, k=4, kmers_keys = false)    
    kmers_keys = get_kmers(["A","T","C","G"], 4)
    kmers_values = [0] * kmers_keys.length
    kmers_comp = Hash[kmers_keys.zip(kmers_values)]
    for i in (0...(dna_string.length - k + 1))
        kmers_comp[(dna_string[i...(i+k)])] += 1
    end
    if kmers_keys
        return kmers_comp
    end
    return kmers_comp.values
end

#Given a protein sequence it return every location where the N-Glycosylation motif is found
def find_nglycosylation_motif(prot_seq)
    pos = []
    for ac_pos in (0...(prot_seq.length - 4))
        if prot_seq[ac_pos] == "N"
            if prot_seq[ac_pos + 1] != "P"
                if prot_seq[ac_pos + 2] == "S" || prot_seq[ac_pos + 2] == "T"
                    if prot_seq[ac_pos + 3] != "P"
                        pos = pos.concat([(ac_pos + 1)])
                    end
                end
            end
        end
    end
    if !(pos.nil?)
        return pos
    end
end

#Given a protein IDs array it searches the aminoacidic sequence at uniprot and give a hash with all the protein ids
#who has N-Glycosylation motifs in his sequence and the locations of the motifs 
def get_prots_with_nglyc_motif(prots_id_array)
    prots_hash = {}
    motifs = {}
    prots_id_array.each do |id|
        prots_hash[id] = (read_online_fasta("http://www.uniprot.org/uniprot/#{id}.fasta")).values[0]
    end
    prots_hash.each do |id, seq|
        pos = find_nglycosylation_motif(seq)
        if !(pos.empty?)
            motifs[id] = pos
        end
    end
    return motifs
end

#Given a protein sequence it returns the monoisotopic mass in dalton
def get_prot_mass(prot_seq)
    weights = {"A" =>71.03711,"C"  => 103.00919,"D"  => 115.02694,"E"  => 129.04259,"F"  => 147.06841,
        "G"  => 57.02146,"H"  => 137.05891,"I"  => 113.08406,"K"  => 128.09496,"L"  => 113.08406,
        "M"  => 131.04049,"N"  => 114.04293,"P"  => 97.05276,"Q"  => 128.05858,"R"  => 156.10111,
        "S"  => 87.03203,"T"  => 101.04768,"V"  => 99.06841,"W"  => 186.07931,"Y" =>  163.06333,
        "H2O" => 18.01056 }
    total_w = 0
    prot_seq.each_char do |ac|
        total_w += weights[ac]
    end
    return (total_w.round(3))
end

#Given a list of DNA IDs and sequences, it returns the directed graph as an adjacency list
#obtained by verifying that the tail of sequence X overlap with the head of sequence Y with 
#length given by window's variable
def get_direct_graph_associations(dna_ids_array, dna_seqs_array, window )
    assoc = []
    for tail in (0...dna_seqs_array.length)
        for head in (0...dna_seqs_array.length)
            if tail != head
                if dna_seqs_array[tail][((dna_seqs_array[tail].length - window)...(dna_seqs_array[tail].length))] == dna_seqs_array[head][(0...window)]
                    assoc.push([dna_ids_array[tail],dna_ids_array[head]])
                end
            end
        end
    end
    return assoc
end

#Given a certain DNA sequence, return all the theorical ORF (from ATG to STOP codon)
#found in that sequence and his reverse complement, in the 3 different offsets
def search_orf(dna_seq)
    rna_seq = get_rna_transcript(dna_seq)
    rna_revcom = get_rna_transcript(get_reverse_complement(dna_seq))
    orf1 = []
    orf2 = []
    orf3 = []
    orf4 = []
    orf5 = []
    orf6 = []
    
    for i in (0...(dna_seq.length - 2))
        if ((i+3)%3) == 0     
            if (rna_seq[(i...i+3)]).length == 3
                orf1.push(rna_seq[(i...i+3)])
            end
            if (rna_seq[(i+1...i+4)]).length == 3
                orf2.push(rna_seq[(i+1...i+4)])
            end
            if (rna_seq[(i+2...i+5)]).length == 3
                orf3.push(rna_seq[(i+2...i+5)])
            end
        end
    end
    for i in (0...(dna_seq.length - 2))
        if ((i+3)%3) == 0
            if (rna_revcom[(i...i+3)]).length == 3
                orf4.push(rna_revcom[(i...i+3)])
            end
            if (rna_revcom[(i+1...i+4)]).length == 3
                orf5.push(rna_revcom[(i+1...i+4)])
            end
            if (rna_revcom[(i+2...i+5)]).length == 3
                orf6.push(rna_revcom[(i+2...i+5)])
            end
        end
    end

    orf1 = orf1.join("")
    orf2 = orf2.join("")
    orf3 = orf3.join("")
    orf4 = orf4.join("")
    orf5 = orf5.join("")
    orf6 = orf6.join("")
    orfs = [orf1, orf2, orf3, orf4, orf5, orf6]
    prots = []
    orfs.each do |orf|
        for j in (0...(orf.length))
            if ((j+3)%3) == 0
                if orf[(j...j+3)] == "AUG" 
                    prot = get_prot_seq(orf[(j...(orf.length))], stop_required = true)
                    if !(prot.nil?) && !(prot.empty?)
                        prots.push(prot)
                    end
                end
            end
        end
    end
    return prots.uniq
end

#Given a list of sequences(reads) of the same length and a threshold of nucleotides length for the overlap
#(a percentage from 1 to 100), i.e. the amount of nucleotides length to validate the overlaps, returns 
#a superstring containing all the overlapping sequences, that is, the genome assembly of the reads, with a 
#maximum error(percentage from 1 to 100)to validate each pair of overlaps
def assamble_reads(dna_seqs_array, min_overlap = 50, max_error = 40)
    similarity_index = (100-max_error).to_f / 100
    total_read_length = dna_seqs_array[0].length
    matches = {}
    superstring = ""
    for head in (0...dna_seqs_array.length)
        skip = false
        for tail in (0...dna_seqs_array.length)
            if skip
                break
            end
            if head != tail && !(matches.values.include?(tail))
                if check_overlap(dna_seqs_array[head], dna_seqs_array[tail], min_overlap, similarity_index)
                    matches[head] = tail
                    skip = true
                end
            end
        end
    end

    start = nil
    stop = nil
    for pos in (0...dna_seqs_array.length)
        if !(start.nil?()) && !(stop.nil?())
            break
        end
        if matches.keys.include?(pos) && !(matches.values.include?(pos))
            start = pos
        elsif !(matches.keys.include?(pos)) && (matches.values.include?(pos))
            stop = pos
        end
    end

    while start != stop
        indexx = check_overlap(dna_seqs_array[start], dna_seqs_array[matches[start]], min_overlap, similarity_index, return_index = true)
        superstring += dna_seqs_array[start][(0...indexx)]
        start = matches[start]
    end
    superstring += dna_seqs_array[stop]
    return superstring
end

#Support function for assamble_reads
def check_overlap(head, tail, min_overlap, sim_index, return_index = false)
    min_overlap_length = ((([head, tail].min()).length) * ((min_overlap.to_f)/100)).to_i
    for window in (0...min_overlap_length)
        if get_similarity(head[(window...head.length)], tail) >= sim_index
            if return_index
                return window
            end
            return true
        end
    end
    return false
end

#Given two DNA sequences it returns a percentage of similarity between the two sequences. It can be used alone
#or as a support function for assamble_reads
def get_similarity(seq1, seq2)
    total = 0
    for nt in (0...(([seq1,seq2].min()).length))
        if seq1[nt] == seq2[nt]
            total += 1
        end
    end
    return (total.to_f / (([seq1,seq2].min()).length))
end

#Given two DNA sequences, it calculates the transition/transversion ratio, that is, the number of mutations of the type
#from purine to purine (ex: G -> A) or pyrimidine to pyrimidine (T -> C) versus the number of transversion, that is, from
#purine to pyrimidine or viceversa (ex: G -> C) 
def calculate_transition_transversion_ratio(seq1, seq2)
    transi = 0
    transve = 0
    for nt in (0...(([seq1,seq2].max()).length))
        if seq1[nt] != seq2[nt]
            if (seq1[nt] == "G" && seq2[nt] == "A") || (seq1[nt] == "A" && seq2[nt] == "G")
                transi += 1
            elsif (seq1[nt] == "G" && seq2[nt] == "T") || (seq1[nt] == "T" && seq2[nt] == "G")
                transve += 1
            elsif (seq1[nt] == "G" && seq2[nt] == "C") || (seq1[nt] == "C" && seq2[nt] == "G")
                transve += 1
            elsif (seq1[nt] == "T" && seq2[nt] == "C") || (seq1[nt] == "C" && seq2[nt] == "T")
                transi += 1
            elsif (seq1[nt] == "T" && seq2[nt] == "A") || (seq1[nt] == "A" && seq2[nt] == "T")
                transve += 1
            elsif (seq1[nt] == "C" && seq2[nt] == "A") || (seq1[nt] == "A" && seq2[nt] == "C")
                transve += 1
            end
        end
    end
    return (((transi.to_f) / transve).round(11))
end

#Given a protein sequence it return the number of theorical mRNAs that may
#produce the same protein
def calculate_number_of_possible_RNAs_from_prot(protseq)
    result = 1  #Variable for returning product of each aac number of codons
    stop = 3 #Variable for remembering that we have to take account of the stop codons too
    codons ={"GCU"=>"A","GCC"=>"A","GCA"=>"A","GCG"=>"A","GUU"=>"V","GUC"=>"V",
             "GUA"=>"V","GUG"=>"V","GAU"=>"D","GAC"=>"D","GAA"=>"E","GAG"=>"E",
             "GGU"=>"G","GGC"=>"G","GGA"=>"G","GGG"=>"G","AUU"=>"I","AUC"=>"I",
             "AUA"=>"I","AUG"=>"M","ACU"=>"T","ACC"=>"T","ACA"=>"T","ACG"=>"T",
             "AAU"=>"N","AAC"=>"N","AAA"=>"K","AAG"=>"K","AGU"=>"S","AGC"=>"S",
             "AGA"=>"R","AGG"=>"R","CUU"=>"L","CUC"=>"L","CUA"=>"L","CUG"=>"L",
             "CCU"=>"P","CCC"=>"P","CCA"=>"P","CCG"=>"P","CAU"=>"H","CAC"=>"H",
             "CAA"=>"Q","CAG"=>"Q","CGU"=>"R","CGC"=>"R","CGA"=>"R","CGG"=>"R",
             "UUU"=>"F","UUC"=>"F","UUA"=>"L","UUG"=>"L","UCU"=>"S","UCC"=>"S",
             "UCA"=>"S","UCG"=>"S","UAU"=>"Y","UAC"=>"Y","UAA"=>"STOP",
             "UAG"=>"STOP","UGU"=>"C","UGC"=>"C","UGA"=>"STOP","UGG"=>"W"}
    
    protseq.each_char do |aac|
        count = 0 #it saves the number of codons that produce a certain aminoacid
        codons.each do |codon, translated_codon|
            if translated_codon == aac
                count += 1
            end
        end
        result *= count
    end
    return (result * stop) % 1000000
end
