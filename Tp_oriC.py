import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def get_complementary_sequence(dna_template_sequence):
    dna_complements = {'A':'T','T':'A','C':'G','G':'C'}
    dna_complementary_sequence = ''
    for nucleotide in dna_template_sequence:
       comp_nucleotide = dna_complements[nucleotide]
       dna_complementary_sequence = dna_complementary_sequence + str(comp_nucleotide)

    #Flip the sequence to read 5' -> 3'
    reversed_dna_complementary_sequence = dna_complementary_sequence[::-1]

    return reversed_dna_complementary_sequence

def calculate_hamming_distance(string1, string2):
    hamming_dist = 0
    try:
        for i in range(0,len(string1)):
           if string1[i] != string2[i]:
              hamming_dist +=1
        return hamming_dist
    except:
       print('String1 and string2 do not have the same lenght')

def neighbors(pattern, d):
    nucleotides = {'A','C','G','T'}
    if d == 0:
       return [pattern]
    if len(pattern)==1:
       return nucleotides
    neighborhood = []
    suffixNeighbors = neighbors(pattern[1:],d)
    for text in suffixNeighbors:
        if calculate_hamming_distance(pattern[1:],text)<d:
            for nucleotide in nucleotides:
               neighborhood.append(np.concatenate((nucleotide,text),axis=None))
        else:
           neighborhood.append(np.concatenate((pattern[0],text),axis=None))
    txt_neighborhood = []

    for i in range(0,len(neighborhood)):
        txt_neighborhood.append(''.join(neighborhood[i]))

    return txt_neighborhood

def skew_counter_func(genome, plot=None):
    skew_i = 0
    c_counter = 0
    g_counter = 0
    skew_val = [skew_i]
    for i in range(0,len(genome)):
         if genome[i] == 'C':
             c_counter +=1
         elif genome[i] == 'G':
             g_counter +=1
         skew_i = g_counter-c_counter
         skew_val.append(skew_i)
    if plot==True:
        plt.plot(np.arange(len(genome)+1),skew_val,'-o')
        plt.xlabel('Genome position')
        plt.ylabel('skew (G-C)')
        plt.show()
    return skew_val

def max_freq_patterns_missmatches_reversecomplements(sequence, k, d):
    freq_map = {}
    len_pattern = k
    max_patterns = []
    for i in range(0,len(sequence)-len_pattern+1):
        actual_pattern = sequence[i:i+len_pattern]
        reverse_pattern = get_complementary_sequence(actual_pattern)
        neighborhood_pattern = neighbors(actual_pattern,d)
        neighborhood_reverse_pattern = neighbors(reverse_pattern,d)

        for j in range(0,len(neighborhood_pattern)):
            neighbor = neighborhood_pattern[j]
            if neighbor in freq_map:
               freq_map[neighbor] += 1
            else:
               freq_map[neighbor] = 1

        for k in range(0,len(neighborhood_reverse_pattern)):
            neighbor_ = neighborhood_reverse_pattern[k]
            if neighbor_ in freq_map:
               freq_map[neighbor_] += 1
            else:
               freq_map[neighbor_] = 1

    m = max(freq_map.values())
    for key_pattern in freq_map:
        if freq_map[key_pattern] == m:
           max_patterns.append(key_pattern)

    return max_patterns


def remove_spaces(file_path):
    # Read the contents of the file
    with open(file_path, 'r') as file:
        content = file.read()

    # Remove spaces from the content
    content = content.replace(' ', '')

    # Overwrite the file with the updated content
    with open(file_path, 'w') as file:
        file.write(content)


data = pd.read_csv('./Data/Salmonella_enterica.txt',delimiter='\t',header=None)
sequence2 = data.values.flatten()
genome = ''
for i in range(0,len(sequence2)):
   genome += sequence2[i]
print(len(genome))

genome_skew = skew_counter_func(genome,True)

min_skew_value = min(genome_skew)
min_sekw_index = genome_skew.index(min_skew_value)
len_window = 425
oriA_window = genome[min_sekw_index-1-len_window:min_sekw_index-1+len_window]

print(max_freq_patterns_missmatches_reversecomplements(oriA_window,9,1))