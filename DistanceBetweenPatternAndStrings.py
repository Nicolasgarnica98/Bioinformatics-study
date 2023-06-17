
def calculate_hamming_distance(string1, string2):
    hamming_dist = 0
    try:
        for i in range(0,len(string1)):
           if string1[i] != string2[i]:
              hamming_dist +=1
        return hamming_dist
    except:
       print('String1 and string2 do not have the same lenght')

def distanceBetweenPatternAndStrings(pattern,dna_array):
    k = len(pattern)
    for dna_string in dna_array:
        motif_array = []
        kmer_hd = 400
        for i in range(0,len(dna_string)-k+1):
            temp_kmer_hd = calculate_hamming_distance(dna_string[i:i+k])
            if temp_kmer_hd <= kmer_hd:
                kmer_hd = temp_kmer_hd
        motif_array.append(kmer_hd)

