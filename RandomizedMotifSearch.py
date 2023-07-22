import numpy as np

#Hamming distance function
def calculate_hamming_distance(string1, string2):
    hamming_dist = 0
    try:
        for i in range(0,len(string1)):
           if string1[i] != string2[i]:
              hamming_dist +=1
        return hamming_dist
    except:
       print('String1 and string2 do not have the same lenght')

#Calculate the profile matrix score
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

#Get the profile matrix
def get_profile_matrix(motif_matrix, LaPlace_rule_succession = False):
    for i in range(len(motif_matrix)):
        motif_matrix[i] = list(motif_matrix[i])
    motif_matrix = np.array(motif_matrix)
    profile_matrix_dict = {'A':np.zeros(len(motif_matrix[0])),
                           'C':np.zeros(len(motif_matrix[0])),
                           'G':np.zeros(len(motif_matrix[0])),
                           'T':np.zeros(len(motif_matrix[0]))}
    
    if LaPlace_rule_succession == True:
        pseudocount_value = 1
    else:
        pseudocount_value = 0
    for i in range(0,motif_matrix.shape[1]):
        actual_column = motif_matrix[:,i]
        A_count = np.count_nonzero(actual_column == 'A') + pseudocount_value
        C_count = np.count_nonzero(actual_column == 'C') + pseudocount_value
        G_count = np.count_nonzero(actual_column == 'G') + pseudocount_value
        T_count = np.count_nonzero(actual_column == 'T') + pseudocount_value
        total_count = A_count + C_count + G_count + T_count
        profile_matrix_dict['A'][i] = np.round(A_count/total_count,2)            
        profile_matrix_dict['C'][i] = np.round(C_count/total_count,2)            
        profile_matrix_dict['G'][i] = np.round(G_count/total_count,2)            
        profile_matrix_dict['T'][i] = np.round(T_count/total_count,2)
        A_count, C_count, G_count, T_count = (0,0,0,0)
    profile_matrix = np.array([profile_matrix_dict['A'],profile_matrix_dict['C'],profile_matrix_dict['G'],profile_matrix_dict['T']])
    return profile_matrix_dict

def Profile_Most_Probable_kmer(dna_string,k,profile_matrix):
    '''''
    This function will calculate the most probable k-mer present on the dna sequence given a matrix of probabilities of nucleotide occurrence(profile)
    Parameters: 
    - dna_string -> String of dna to check for the most probable k-mer.
    - k -> size of the k-mers to form
    - profile_matrix -> matrix of probabilities of nucleotide occurrence.
    '''''
    MostProbableKmer = dna_string[0:k]
    MostProbableKmer_prob = 0
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
i = 0
 
def RandomizedMotifSearch(dna_array,k,max_iter):
    i = 0
    best_motifs_array = []
    scores_array = []
    while i<max_iter:
        best_motifs = []
        for dna_string in dna_array:
            random_index = np.random.randint(len(dna_string)-k+1)
            random_motif = dna_string[random_index:random_index+k]
            best_motifs.append(random_motif)

        temp_motifs = []
        counter = 0
        while 1<2:
            profile_matrix = get_profile_matrix(best_motifs, LaPlace_rule_succession=True)
            for dna_string in dna_array:
                temp_motifs.append(Profile_Most_Probable_kmer(dna_string,k,profile_matrix)[0])
            if calculate_motif_matrix_score(temp_motifs) < calculate_motif_matrix_score(best_motifs):
                best_motifs = temp_motifs
                temp_motifs = []
            else:
                break

        best_motifs_array.append(best_motifs)
        scores_array.append(calculate_motif_matrix_score(best_motifs))
        index = scores_array.index(min(scores_array))
        final_best_motifs = best_motifs_array[index]

        i+=1
    
    print('Min score: ', min(scores_array))
    return final_best_motifs


sequences = 'CCGGTCCGTTGCTTCTCGTACAAAGGGACGAAGAGCTTCAACACATAATCAATAGATGTCATCGAAATTATCTGAAGTAGGAGTACCTAGCCTGATCGGGGAGCACTACAGAATATTTAGAAATGCAGGATTCCATGACAATGGCCTAATAGAGGCCTGTGCTTGGAGAAGTCTATAAGGCCAAACAGTAGTTTAACCGGTCCGTTGCTTC TCGTACAAAGGGACGAAGAGCTTGCAAACTTGAAGGCTCAACACATAATCAATAGATGTCATCGAAATTATCTGAAGTAGGAGTACCTAGCCTGATCGGGGAGCACTACAGAATATTTAGAAATGCAGGATTCCATGACAATGGCCTAATAGAGGCCTGTGCTTGGAGAAGTCTATAAGGCCAAACAGTAGTTTAACCGGTCCGTTGCTTC CGGTGTAGTTTGTGCTCTCGCACCGTTGAAGTTCAGTATAGAAAAATCAGGCATCCTCGGAGTGTGCTAGATTTGATTTATCGGTTAAAAAATGGTGCTAGACCACGAGTGGAACAACACGAAGCGAAAGACGCTCGTGATTGGTTACGCGTTCGAGTGCATTGTTCTGTTAGATACAGAGTCCGAACCGACCAATTATATAACTTTGTCG CTGAACGCTAGCGTTGAACGGAGAGGGATGAATGGAATCGCTTGGCTAGATGAGGGCGTCCCTTAATTTGATGTTGGACGGGAGGTATGGTCATTAGAGCTTACGCCACCCGACATGAACCTCATCCACCCTTAATTCTACATTGAGCTAAGTTGAAGGCTATATGTAAAGCGCCGAGGTCACGGGGACGGCCTGTGAGTGAAAACGTGTC GTACATCCAGTAGATAAAGGCCCCCAGACAATAGAAACTTTCCTGCACAGTGTTAAGCTAGCTTCGTAAGATAGGTTGCAAAGTAATTTTAGTCAAGTGGATATAGCGGCAGTCACTCGCAAAAAGGGACAGTGCACCAAAGAAGGCTGATGTAGGGGCACACCTGATGCTAACCGCTATACATAGAGAATTGCGCTATTGAGTTCGTGCC TGTGCAGTACCGAAGACGTCATTTACCTGTGTGGTCACCAAAGGGAATATCCCCAGCGAACTGCCTGGCGGGACTGTAGCCGCGGCGGCGGGTTGTCAAGCACTCACGAGCTACGAGGAGTGCCGCACCGTTGTGTGCTATCTACCAGCACAATTTTAGGTCGTGCCATCTACCACAAACGCGGTTTAGACGGCCTGATTTCCTTTCACTT AAGCTCTAGCTCGACGTACACGCGGGAAGCTAAATTATACTGCAACGCGGGATAGATTCCATGGACCTCGAGAATTCTTCACTCGGTAAAGGGCTTATCTCGCAGCCACAAGATAGCCCATTTGGCTTTAATATACGGATGCGGTAGTCAGGACGTTGCAACAGTGTGGCCGTTGAAGGCTGTGCACGAATTCGGTGCGAATTAACGTTGG AATTACTATAAGGTGCTTTGGGGAAAGCTCAGACTCGACATGGTCCTTGGAATTTGTGGCACTAGGACTAAGCAAATTTTGTGCCTACTGCACTACCGGGGACCGTACTACCTTGACGTGAAGAAAAAACTCATCTATTACTCCTGTTGGTCAGACCGTAATATTGGGGCACCGTTGAATATTAACCTCCATCTACGAAGTGGGATTTACG TGCCTCGCGACCGGATCTAAGCGGGTAAGAGGCGTTGAAGGCTGCGAGGGGTAAACTTAGAGGTAGGCACTTAGATGCAGAGACTCTATTCACGTAATTCCCCGCTACCCATCGCACAGGCGAGCTCCTCGCCACCAGACATGCGCTGGGTTAAACGTCACAGACGGTAGCCTATTCTGCTCCCGTTAGACGCAACATCCACATTCCACAC GACGTATAAGATTCCGGGATCCTTGTATCCGATCGATAGCGACCTACGTACGTCAGCGCCCTGAAACGGAGCCCGGGTCTAGGCTGACAACCTGTTTGGGCGGCTATCATATAAACTCTAGTCAGCGCGGGTATTAGACAAATCAGCACGCGTGAAGGCTCCATCAGCCAGGTCGGACAACCGCGTATGAGATACGACTCGGAAACCACGA GCGACTTAGCTAATGATTATTGGATGAGTGGTCATACGGGTTGCTACTTTCTTGTGTCCGACGAGTACCGAGTGCCGTCGGCACCTCAGAAGGCTAAGAGACTCTAACATCATCCCAAGACTCATTGGGCCGGGACGTGATGTTCCTTTTCTGGAGTCAGAGAGCCCTTTCGAGACTGCCTCCATGAGGCAAGTAAGTCAAACTTACGTAC TGCCGCGCTGTTACCGGGTCTGCACCGTTTTCGGCTAGACTAGGCGTCGCGACAATCGTGGCGCACACGAGGCCCCGCTCTCTATGCTAGAGCTCTATGGGCAAATGGGCTGCTAGTGCTCGGATCAACACTTTCATGGGTCCAGGATCCGTTGTATACCACCCAGGAACAAGCACAACTTGCGTTATATCAGTTTCTATCCAATAATTGG ATCAAGATAATACGGAGCCAAGAAGGATAGCGGACGGGCACCGTTGAGAACTGAAAGTAGAATCGGAGTTAATAAATGATTCAGATTCAGACAGGACCGACGGCGGCCCCAATCCTATCCTTAAGTCTGATCGAAGATTAAGTGGATGGATATGGCCCTATAGACTAATTAAGCAACGCAAAGAAGGGGAAGTATCGGGATAAGCTATCCG ATGCGTTTCTTTTCGACGGTGTAGCAAAAAAATCAGCACGCAGTTCGAGGTAATTGTTTGGCCCTCCTGCGGTAACTAGGCTAGCATCCTGGCGCTCAAGTCTCAAGACTGGGCACCACTCTGAATCGTCACTCGAGTCAGGTCAAATAATGCACGCATGTTTGAAGGCTGCACTTCCCTGCCCGTCACCGGTCCTTGCCTGTATCTTGAG AACAGCGTGGACGCCTCTCCGGGCAACCGTAGCCCGTCAAGAGTAATGGCTGCTACCGCTTTTATGATTGATCTATGCCTTTACCCCCGCCCTCGCACCGGAAAAGGCTCTAGCGTACCTTTTATCCTTGACGAACTACTGACTCCGGGAGCGTGGTGGGGTTCCTCCTGCGACCATCGTTTGCGGCGACGCAAGTAAAGGACGCCCGGGC CGGGCCGCTTAAGAACGTGTAGAGACAACCGGGAGCTGGGTTCGGACGACAGGATGATTCGTCACTACCGTTGAAGGCGCCGCAAAACCCCCTCTCCAAGTAGCTTCCTGAGTTCATCCATACAGGAGGAGGCGCCTGCCTATGATGTATCTTTCCTTTTATAACAGGTGACTGAAATAGTTCATTCTAAAAAGTTTTGACTTGCAGTGCA AGGGACTATATTCATCCTAGTATGACTATGGGCTCACACGGAGCGACTTAAGTTGTTCGCGGCACCGGCAAAGGCTCCATTACGACCTCCGAAATTCTCGCGTCATCTGACTAGGGTGATAAGTTGTTACGGGGAGGTTACTAGCTCTCGGCCATAGGGGATGCGTGCCTTTATAACAAGGGTTATATAGGCCTCACCGATCGCAGTTGGG CCCGGGGTCGTTCTGATCACACCGTTGAAGGAAACAATCCAGCACCTTGCAACATTACTGCTAGGTGTCCGACAGTCGTCTTCGGGGGCTCGCAACGTGATGTTAAGAGAACTTACACATTCAGATCGTTTATGTCGTCGCCAGGTTAATATAGCGGCCAGCCCGCTTAAGGGCATGAAGCGGCTCTACGAGGTCAAAGAACTCAATCCCA GATTTGTTTGGAGTACTGTGCGAGGAACATGAGGCTTACTGGGAAGGAGATCCCCAGCTTCCGCCGAAAACAACGCAGGCTGTACAAGTTAACGCACGAATGAAGGCTACGCGTGAGCGTTAATCAGCATGAGCCGGTGCAATCTCGCAATGTTAGATTGAGCCATGCTCCAATAATTCACTATACTTTTATCGTGAGCCTATAGCTAATC CTTGCTATATAGTTTGGGCATACACACCTACCAAAGGGATCTGCACCGTCTGAGGCTGACTTATGCCCCGGTCACGGGCTGTGGTGCAAATCCCACGCCTCCTACGCATTCGAGAATACCGGCACGATACGCTCACCAAGTACGAAACCGAGGTTTGACAGGACTTTCGTACGGGAACAAATCTGTTTAGAGATGAAACGATGCCTTATCT'
sequences = sequences.split(' ')
result = RandomizedMotifSearch(sequences,15,1000)
for i in range(0,len(result)):
    result[i] = ''.join(result[i])
print('Best motifs:' ,' '.join(result))