import numpy as np
from tqdm import tqdm

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

#Median string problem
'''''
    The aim of this function is to find the median string (pattern) that minimizes the distance d(pattern,dna) by getting the motifs.
    1st - we have to get the pattern by sliding window trhough all possible patterns(k-mers) in the dna strings.
    2nd - check the minimum Hamming distance between that pattern and the rest of the possible k-mers in all the strings. 
          Choose the one k-mer that minimizes the Hamming distance by each dna string.
    3rd - Sum the minimum hamming distances from each dna string to get d(pattern, dna)
    4th - Among all possible patterns, check which of them have the smallest d(pattern, dna)
    5th - That pattern will be the Median string, and the k-mers will be the collection of motifs.


    dna_array: Collection of DNA sequences of same lenght
    k: lenght of the k-mers to form patterns and motifs

'''''
def MedianString(dna_array,k):
    median_string = ''
    possiblePatternsArray = []
    for dna_string in dna_array:
      for i in range(0,len(dna_string)-k+1):
        kmer_neighborhood = neighbors(dna_string[i:i+k],int(len(dna_string[i:i+k])/2))
        for kmer in kmer_neighborhood:
            possiblePatternsArray.append(kmer)
    
    d_pattern_dna = []
    for pattern in possiblePatternsArray:
        possibleMotifCollection = []
        possibleMotifCollection_minValues = []
        for dna_string in dna_array:
            kmer_collection_temp = []
            kmer_hd_temp = []
            for i in range(0,len(dna_string)-k+1):
                kmer_collection_temp.append(dna_string[i:i+k])
                kmer_hd_temp.append(calculate_hamming_distance(pattern,dna_string[i:i+k]))
            min_hd_kmer = min(kmer_hd_temp)
            min_hd_kmer_str = kmer_collection_temp[kmer_hd_temp.index(min_hd_kmer)]
            possibleMotifCollection.append(min_hd_kmer_str)
            possibleMotifCollection_minValues.append(min_hd_kmer)
        
        d_pattern_dna.append(np.sum(possibleMotifCollection_minValues))
    
    min_d = min(d_pattern_dna)
    median_string = possiblePatternsArray[d_pattern_dna.index(min_d)]

    return median_string 


excercise_dna_list = 'GTATCTGCTACCGAAAACTGTAGTCAGTATGGCATTTATGGT TCTCTTATATCAGCAAACGATTTAGCCCGTTCTGACTGAGGA GGAAACGTATTGGATGCACCGGGGCAGCTGTGGGAGAATTCC CAAATACGCAATGTGGGTGTACAGAGAGATGAAAACTTTTTT CTGGTCTGATTCTATCTAACGACCGGTGGAGAAAACAGGCAG AGTAAATAAACGGCAAACGATTTTGTATGACTCAAACAACCT CTGAAAGCTTAGTCCCCACCCCCCGAAAACCACCTACCTAGT CTATTAGTAAACGGGCCTCAACTAATAGTGCGCCGCTAAACC GTAGGTCAAGGCCGCCCGAGATTGTTAAATCATGAGGAAAAC CGGGAACATCACGGAAACACTGTGATGTTAACCTTTCATTTT'
excercise_dna_list = excercise_dna_list.split(' ')

print(MedianString(excercise_dna_list,6))