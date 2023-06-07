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


#Motif enumeration function

'''''
    dna_array: list of dna strings
    k: lenght of the k-mer
    d: max number of mismatches

'''''

def MotifEnumeration(dna_array,k,d):
