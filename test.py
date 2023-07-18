import numpy as np

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
    return profile_matrix

test_input = ['ACTCATCGAACTGG']



print(' ')
print(test_input)
print(' ')
print(get_profile_matrix(test_input))
print(' ')

#%%
import numpy as np
def calculate_hamming_distance(string1, string2):
    hamming_dist = 0
    try:
        for i in range(0,len(string1)):
           if string1[i] != string2[i]:
              hamming_dist +=1
        return hamming_dist
    except:
       print('String1 and string2 do not have the same lenght')

def calculate_motif_matrix_score(motif_matrix):
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
            print(calculate_hamming_distance(dna_string,consensus_motif))
            score = score + calculate_hamming_distance(dna_string,consensus_motif)
        print(consensus_motif)
        return score

test_input = ['ACTAGT',
              'ATTGTA',
              'GCTTGA']
print(calculate_motif_matrix_score(test_input))
# %%
