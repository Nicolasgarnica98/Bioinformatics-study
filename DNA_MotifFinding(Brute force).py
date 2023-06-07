import numpy as np

# Hamming distance problem
def calculate_hamming_distance(string1, string2):
    hamming_dist = 0
    try:
        for i in range(0,len(string1)):
           if string1[i] != string2[i]:
              hamming_dist +=1
        return hamming_dist
    except:
       print('String1 and string2 do not have the same lenght')


# Neighborhood of a string
'''
    pattern - sequence
    d - hamming distance max.
    Output: The collection of strings Neighbors(Pattern, d)
    
'''

#By recursion
def Neighbors(pattern, d):
    nucleotides = {'A','C','G','T'}
    if d == 0:
       return [pattern]
    if len(pattern)==1:
       return nucleotides
    neighborhood = []
    suffixNeighbors = Neighbors(pattern[1:],d)
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


#Motif enumeration function

'''''
    dna_array: list of dna strings
    k: lenght of the k-mer
    d: max number of mismatches

'''''

def MotifEnumeration(dna_array,k,d):
    patterns = []
    #No need to iterate on each dna string. If there is a motif, it should be present in all
    #dna strings!
    dna_string = dna_array[0]
    #Check each possible k-mer on that dna string
    for i in range(0,len(dna_string)-k+1):
        act_pattern = dna_string[i:i+k]
        #get the neighborhood of that ith k-mer (Mutations)
        act_pattern_neighborhood = Neighbors(act_pattern,d)
        #Iterate in all the dna strings. Find if any k-mer of the ith k-mer neighborhood is present
        #in all the dna strings with at least d mismatches.
        for pattern in act_pattern_neighborhood:
            count = 0
            for dna_string_inner in dna_array:
                for j in range(0,len(dna_string_inner)-k+1):
                    act_pattern_inner = dna_string_inner[j:j+k]
                    if calculate_hamming_distance(pattern,act_pattern_inner)<=d:
                        count += 1
                        break
            if count == len(dna_array):
                patterns.append(pattern)

    motifs = np.unique(patterns)
    return motifs

# test_dna_list = 'ATTTGGC TGCCTTA CGGTATC GAAAATT'
# test_dna_list = test_dna_list.split(' ')
# result_test = MotifEnumeration(test_dna_list,3,1)
# result = ' '.join(str(element) for element in result_test)
# print(result_test)
# print(result)

excercise_dna_list = 'TGCCCTATTGAAGGGTTGAGTGGTG CCACACATTGTCCTACATACCTTAT GATTGCTGGGGACCCGCGGCAGGCC CTCCTCATTGCACGATCTAATGTTT CAGCCGGAGGGTGACCATTGACTGT TTAGCGATTGCGGCGCCGGATCGAC'
excercise_dna_list = excercise_dna_list.split(' ')
result_excercise = MotifEnumeration(excercise_dna_list,5,1)
result = ' '.join(str(element)for element in result_excercise)
print(result)