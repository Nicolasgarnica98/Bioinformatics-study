import numpy as np
from tqdm import tqdm
#Most probable kmer on a string given a profile matrix of probabilities.

def Profile_Most_Probable_kmer(dna_string,k,profile_matrix):
    '''''
    This function will calculate the most probable k-mer present on the dna sequence given a matrix of probabilities of nucleotide occurrence(profile)
    Parameters: 
    - dna_string -> String of dna to check for the most probable k-mer.
    - k -> size of the k-mers to form
    - profile_matrix -> matrix of probabilities of nucleotide occurrence.
    '''''
    MostProbableKmer = ''
    MostProbableKmer_prob = -1
    temp_prob = 0
    for i in range(0,len(dna_string)-k+1):
        actual_kmer = dna_string[i:i+k]
        prob_array = []
        for j in range(0,len(actual_kmer)):
            if actual_kmer[j] == 'A':
                prob_array.append(profile_matrix['A'][j])
            elif actual_kmer[j] == 'C':
                prob_array.append(profile_matrix['C'][j])
            elif actual_kmer[j] == 'G':
                prob_array.append(profile_matrix['G'][j])
            elif actual_kmer[j] == 'T':
                prob_array.append(profile_matrix['T'][j])
        temp_prob = np.prod(prob_array)
        if temp_prob > MostProbableKmer_prob:
            MostProbableKmer_prob = temp_prob
            MostProbableKmer = actual_kmer
    return MostProbableKmer, MostProbableKmer_prob

sequence = 'GATGCGAATTCAGTTTGACCAGATGTTATTGCATACCTCACAGCCCCACGCCCCCTTCTGCCCCGCCTCGTACAACGGGCCAGTGCTACATGGTCTAGCCCGCTGAGCAGGTAGGACAGGTAGTACCTCCGGGTAGATGGTTCGTGTCTGGTATTTAGAGTGAGCCTGTGTAAGATAGGCGAGTGATGTCTCAAGCCCCGAGGCATGGGGGGTTTCGGTCTAGCTCTGTCGCAAGAGATATAGCTAGCGTCGGCGTGGACATTATCACTGCGTGTGCTTGGTCACACTACGGACGAGGACTCGATAAGTGGCGCATCTGGCGCGCTTAAGTATATATCCTTGAAATCAAAGGGCCATCATCGGCCCTCTTAATCTATTATGGCGCATATACATGTACAAATGTAGTGTCTAAGCGATTCGTCAATTCGCTAGACTCCGTGACTATTGACCGCGTGAGCGGGCGCGCGCGTAAAAAAATTTGCCCGCAAATTAGTCCGTCCTGACACGAGACAACCTGCACTTCACTATGACATACAGTGTCAGGGCCGATAGGACATAGTGGGCATCGGCTTGCCTGCTGGAAAAATCGAAGGTTGAGCCACTCCTTTTGTCGTGATCTCGTGATGCGGATGCCGTAAAATGCTCTCAGTCACCGTGTATACGTGTCCCATCAACGCTATTGAATCCCGATAACTTTTTCCTGAGTACTAGGCCAATCATAAATCCTAATACGGTTGCCTAGTATCTCGCTCATATAATCAGCAGTCTTGCACTAGAGATTAATAACCAGTTGCGCGCCCTGTCTCGTCTATATGGTACAACGGTGCTGAGACAGTCGTACCGTCTGCGCCGTATCAGATCACGACAGGGCGCTAGCACAATGGGGTAGCAAGGGAAATCTGCATCATAGATCAGATAACGACGCCTCCGCTACAATTGAGTACTAGAATGGGCTACTACCCCCTCATCCAGTTGACGGGCTCTTGCCTATGGCAGTTTC'
profile_matrix = {'A':[0.303,0.227,0.273,0.333,0.212,0.273,0.242,0.182,0.227,0.197,0.273,0.258,0.197,0.318,0.242],
                  'C':[0.197,0.348,0.167,0.167,0.258,0.197,0.227,0.273,0.273,0.227,0.318,0.227,0.258,0.197,0.288],
                  'G':[0.227,0.212,0.273,0.227,0.288,0.303,0.227,0.258,0.227,0.333,0.227,0.258,0.318,0.167,0.212],
                  'T':[0.273,0.212,0.288,0.273,0.242,0.227,0.303,0.288,0.273,0.242,0.182,0.258,0.227,0.318,0.258]}

def calculate_hamming_distance(string1, string2):
    hamming_dist = 0
    try:
        for i in range(0,len(string1)):
           if string1[i] != string2[i]:
              hamming_dist +=1
        return hamming_dist
    except:
       print('String1 and string2 do not have the same lenght')

#Greedy Motif Search algorithm

def greede_motif_search(dna_array,k):
    '''''
    This algorithm will return the best motifs by trying all the possible profile matrix formed by the kmers in the dna strings.
    We need to maximize the p value until finding the best kmers that will form the deffinitive motif matrix.
    Parameters:
    - dna_array: Strings of dna
    - k: Size of the k-mers to form -> size of the motifs
    '''''
    #Get the profile matrix
    def get_profile_matrix(motif_matrix):
        for i in range(len(motif_matrix)):
            motif_matrix[i] = list(motif_matrix[i])
    
        motif_matrix = np.array(motif_matrix)

        profile_matrix_dict = {'A':np.zeros(len(motif_matrix[0])),
                              'C':np.zeros(len(motif_matrix[0])),
                              'G':np.zeros(len(motif_matrix[0])),
                              'T':np.zeros(len(motif_matrix[0]))}

        for i in range(0,motif_matrix.shape[1]):
            actual_column = motif_matrix[:,i]
            A_count = np.count_nonzero(actual_column == 'A')
            profile_matrix_dict['A'][i] = np.round(A_count/motif_matrix.shape[0],2)
            C_count = np.count_nonzero(actual_column == 'C')
            profile_matrix_dict['C'][i] = np.round(C_count/motif_matrix.shape[0],2)
            G_count = np.count_nonzero(actual_column == 'G')
            profile_matrix_dict['G'][i] = np.round(G_count/motif_matrix.shape[0],2)
            T_count = np.count_nonzero(actual_column == 'T')
            profile_matrix_dict['T'][i] = np.round(T_count/motif_matrix.shape[0],2)

            A_count, C_count, G_count, T_count = (0,0,0,0)
        profile_matrix = np.array([profile_matrix_dict['A'],profile_matrix_dict['C'],profile_matrix_dict['G'],profile_matrix_dict['T']])
        return profile_matrix_dict

    def calculate_motif_matrix_score(motif_matrix, return_consensus_motif = None):
        motif_matrix_copy = motif_matrix
        for i in range(len(motif_matrix_copy)):
            motif_matrix_copy[i] = list(motif_matrix_copy[i])
    
        motif_matrix_copy = np.array(motif_matrix_copy)
        consensus_motif = []
        for i in range(0,motif_matrix_copy.shape[1]):
            actual_column = motif_matrix_copy[:,i]
            unique_elements, counts = np.unique(actual_column, return_counts=True)
            most_frequent_index = np.argmax(counts)
            most_frequent_element = unique_elements[most_frequent_index]
            consensus_motif.append(most_frequent_element)
        score = 0
        for dna_string in motif_matrix:
            score = score + calculate_hamming_distance(dna_string,consensus_motif)
        if return_consensus_motif == True:
            return consensus_motif
        else:
            return score
    
    #
    best_motifs = [dna_string[0:k] for dna_string in dna_array]
    for i in tqdm(range(0,len(dna_array[0])-k+1),'Analizing: '):
        kmer1 = dna_array[0][i:i+k]
        motif_matrix = [kmer1]
        for dna_string in dna_array[1:]:
            profile_matrix =  get_profile_matrix(motif_matrix)
            kmer_i = Profile_Most_Probable_kmer(dna_string,k,profile_matrix)[0]
            motif_matrix.append(kmer_i)
        
        if calculate_motif_matrix_score(best_motifs) > calculate_motif_matrix_score(motif_matrix):
            best_motifs = motif_matrix

    return best_motifs, calculate_motif_matrix_score(best_motifs,return_consensus_motif=True)

sequences = 'TCCTGATCGTACCTACAGCGAGCCTAAACTGAAATAGCAGACGCTGTCCCTCCAGCTATGAACTTAAGGCGTTTTGAAAACTAGCAACAGACGTATATCACTTCGACCTTGCGCTTGTTGTACGACTACCAGAACTCAACTATGTAGTATAGATGG CTAATTTTGTCGAACGAGTATCTGTAACCGGCGAGTGTGCATACCACGCCCCGTTGGTTTCAGTGGCTGGCCCGAAACTACACGGAGAGCTGGACGACGGGGCCGCCTTAGCCAGACCGTTGTACAACAACGGTACCGGCGCATGAACAAGCTGTA CGGGCGCCTGTTCAAATGCTACTACACAAGGTACCCCTCGTTCAGTCGGGACCCCGATCTCTTGATTAGTAGGACGAATATCTGTAAGCCGGCACGCACGAATTGCGTGCACTCAGCTGAAACAGCCGGAAGAGCCGTGACAGGTCCTGTCTGAGA TGCTACTCCTACACAACTCTGTAGCAAAGGGCGGATCGTGCGGTGCAAGAGCCCTTGTCGACCTTTGTTTTGTACGAGTACCCGCGGTGGGAAAATACTACACGTACGCCCTTTTCATCGATTGACATCGTATGGCGATCGCTCGATATGGTCCCT GACGAATACCTGTTCTTTGACGATGAGGCACTTAAGAGTTCTCATATCCTTAACGCGGTTTCGTAGGCCTACTACACGAGGTCATTCGAGAGTGCTTTAAGAACGTGTGGTAGTCCGATGAAGCTAGGAGACCGCCACACCCGAAATGGGACGCCC TCTAAGCTTGTCTCTATGGTGCCTCACGAATAACGGGGGCCTCGAGTTTCGATCTTGCCCATTAGAGAATTGCATTGCACGATGTGAGCCGCCCCTGTTACGCATCCGTGGACTGGCAGAACAACTCCCTAGCCATTTACTGAGTCTGAGTTCCAT TAAGCGAGCTTACTATCAGGTGTTCTACTCTTCATGTTGGTGTTGGGCTACGATTATCGGACCAGCCATAGCTAGATCCTCGCTGACGCGTGTCGTAGACAGCTCAAACCATAGTCATCTGTGAACAGGATCACACGTCATAAGTCGGTTGTACCC CGGGCGTCGGTTGCGTTTAATAGCTACGGCCACTATTCTAGCGTCACTCACCCGAATCGAAACGAGTATCTGGCTAATCCACTTTGATCACTCGGGGAGGTAATGTCGAAATTCGGTCTCTACCGCAATTACTCGCTCCAGGTATTGTTGGGCGAC TCAATTTAAGGTGATTATACCGTTTACGAGTATCGGTGCCATCTCAGGAACGTAGCCGTGTTGTCCCTGTTATCCAAGGGGTTCGTGCTGCATCCTCCGAACGAGGTATGACGACCATTTTGTCATACTGTTGCCACGGTGCCTTAGCGGTAGACA GGAGGCCCTTGTCCTGAAGTAACTCTTTCGAAGATCCGCTTACGCACTAGGCATAAAATTAGGTTTCCCTCCGGGTAAAGGCATCAACCACTAACGTGGTCGCATGGATATCACTGATGATAGCAAGACAGCAGCGGAACGGTCAACGACTAACAG GATCGGTTAGCTCCCAGGTGTGCCCTCTTCGGGCTTCTGTCGATCCCCACGGAAAACTCCGACGAATACCTGTCTTCTGTTTATACTGCGCCGCGAAGTGGCTCGTCTTATTCAAAACAGTCTGTGCGCGACGGAACCAGAAGTGCTCAAGGCTGC CCAGGCATGCTCCGGGGGTAACCAGTACCATCTTGTACTGAACAGAACATTGGGAAAAGCAGTTCTACTTTGTACGAGTAACGGCCCTGCTACATTGTGTGTGTATGTGTCGAGTTGCTCGCGCCGGCTAGCCTCCGCCGAGTTAATTGGACTCGA CGCCTACAACGCCCCAGCCCCCCAACTCCCACGATAGCTCGCCTACCCAACCCTGGCGCTTCTAACGCATGTTGACATGAATCTCTGCTGAAACATACATTCATGGTCACGGGACCGGACCATTTCTGCTGGTACGAGTACCGGGGTTCCTGTGAC CTCTACGTACAGGGACCAGGACAACGGTAAAGGAAGACAGTCAAGGCGAGAATCCCTCCCTCCGACATCGTCATGCAGGTGAGCACAGTGAACGCTATGAAGTCCTATGACGAATAACAGCCTAGTGCTTCCGCAGGGACAGATGTACAAATCTTC AAGAACGTCTAACCTAGAACACTTTTGAAACGGTTCTGGCGCCACTACCACGAGTATCAGCAGCGGGAGCTTTCCGCTCAACTCTGCTAGATCTTAGTGCACGCTCAGGCGCTGTCCGATCCATGTCGTCAGCGATGTAAGATTGTGGGTAAACTA CCTCTACTGATGAATCGCAGGGGTGGTCGCGTATCAAACGATTACCGGTACTGGCGGATTAGGCTCGTTTAGGAGTGTGGACTTGAATGTGCAATTCCGCTGGAAGCAGCAAACTCCCCTATGCGCCTGTAAGGACCTGCCTCTGCACGGCATCGA TGGAGTGCGCCCAGGTCAGCGCAGCCCCAATTGAATCTAACGCTGAGCTACGAATACCGGCGCACAAGCGTCCAGGCAGCGTAAAAACATGTCCAGTGTGACGAGGAAATTTTAATTGATTGGTAGATGCTCAGATTGTGGACTGAAATCGAACAT AACAGACTCCTAGGCAGCGTATCGAGCAATAATGGCGATATTTTTCCATACCTTCCTATCAACAATGAGACGAACGAGGTCTGATGCGCGGAGTAGTAGTGCGGCTCGCCCCGTGCGAAAAAGTGCTACATAACGGCTGCGGATCACGATTAGCCG CTTGCTGCTAATAAAATCAAGTTCACTAGAATTTGTAGATGATTAATCCACGAATATCTGCTCTTGACTAGGTGGTATCTTGCTCCAGGTAGTGGATCGGTGTAACTCATATTTAGGAAGGCGCGCGAAGTGGTTCACCCACAGTTTGTTAGCGTA ATGCCCTGTTCTAGAAATGCCCTAGCTGGGGTGGGACCCAATAATGTTCAGATGAAAAATAGCTTCGTCTCGGGAGCGTGCCTAGACGACTAACAGCACACTAAAAAATCATGACCGTTCGCTGCGGCGGAGAGATCTACACCACGCGATGGTACG GAGCGTTCGTCATTGCGTGAGCTATTCGCACAAAACCATACTCCTAAATAATGTTGTCAGTGTTGTGTACTTGTTTGACCGAGCATTGTAGACTTGACTTGTGTCATTCCCTCGGGAGCGTCTGTTGCGTACGACGAATAACAGGGAGGGCCCGCG GCTTGCATTGTAGTGAAATTGAGCTCGATAGAGTTTCGAACAATTCTCGCAAATTACGGGTACGATTATCCGCTTTGACCCCCGTATTCACAAGGGCTTCGGCGTGCGCCCTTCGTGCACCTGGTCGACTGTTCTAAATTCAGTAGACCAATAGTT GAGTCAAGAACGTGTTGGGCTGTCTCAATGAGAACATTGCTTGAGTTATAGAATTACTAAGGGTATCGGGCGATCAGCGCGGCAAGTTTTGTGCAGTTTATTACTCTTACGCTAGGTCCATACGATTAGCTGAATCACGAGGATCTGAATCTAGGA CACACGCGCTTTAATGTGCGTTAGCCTACACGGGCATAGTTATCAAAAACGCATCGCCTAGATAGTAATTCTGCTCAACGACTGGCCTGGTTTAACGCTTTGGTCGGTGTTTAAGTCTTGTACGACTACCCGATAGTGCGCGGAATATGCGCTGTA CCGTACAATTAGCTCTCGACGCATTGATCCCTGTGGTTGGACGCTCTCGAATTATTGGGGGTTTTGTGTAACTGATATCCACACTTACGGTCAGGCCACGAATAACCGGCCCGCACGGAAATAAGATAATCCCCCAATCCGCAATGGCCATTAGAA'
sequences = sequences.split(' ')
result = greede_motif_search(sequences,12)

for i in range(0,len(result[0])):
    result[0][i] = ''.join(result[0][i])
print('Best motifs:' ,' '.join(result[0]))
print('Consensus motif (Implanted motif): ', ''.join(result[0][1]))