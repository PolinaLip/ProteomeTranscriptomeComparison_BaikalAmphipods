#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from Bio.Align import substitution_matrices # to get blosum62 matrix for distance calculation
#from itertools import product # to perform pairwise comparison of keys from dictionary
import scipy
from scipy import spatial, cluster

''' To select regions inside all sequences (dictionary of seq names and aligned regions) '''
def select_region(dict_w_seqs, index, klength):
    region_dict = {}
    for seq in dict_w_seqs:
        region_dict[seq] = dict_w_seqs[seq][0][index:index + klength]
    return region_dict

''' To find distance between two given sequances '''
def pairwise_comparison(seqname1, seqname2, region_dict_w_seqs, score_dict):
    #print(seqname1)
    #print(region_dict_w_seqs[seqname1])
    seq1 = region_dict_w_seqs[seqname1]
    #print(seqname2)
    #print(region_dict_w_seqs[seqname2])
    seq2 = region_dict_w_seqs[seqname2]
    #print('-----------------')
    score = 0
    for i in range(len(seq1) - 1):
        letter_seq1 = seq1[i]
        letter_seq2 = seq2[i]
        score += score_dict[(letter_seq1, letter_seq2)]
    return score 

def main():
    parser = argparse.ArgumentParser(description='To look for the regions in the sequences that determine the inducibility of proteins (ex. Hsp/Hsc70)')
    parser.add_argument('--fasta', type=argparse.FileType(), help='fasta file with aligned sequences of interest')
    parser.add_argument('--features', type=argparse.FileType(), help='file with wanted topology of sequences (based on features)')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), help='output file with regions and associated features')
    parser.add_argument('-k', type=int, help='kmer length (= 6, by default)', default=6)
    args = parser.parse_args()
    
    outp = args.output
    kmer_length = args.k

    ''' Prepare blosum62 dictionary or dictionary with simple scoring system (penalty = 1 if two letters are different) '''
    blosum62_mtrx = substitution_matrices.load('BLOSUM62')
    aas = blosum62_mtrx.alphabet.replace('*', '-')
    blosum62_dict = {}
    for index_row, aa in enumerate(aas):
        for index_col, value in enumerate(blosum62_mtrx[index_row]): 
           blosum62_dict[(aa, aas[index_col])] = value
    
    simple_score_dict = {}
    for letter1 in aas:
        for letter2 in aas:
            if letter1 == letter2:
                simple_score_dict[letter1, letter2] = 0
            else:
                simple_score_dict[letter1, letter2] = 1

    ''' Prepare dictionary with sequences and their names '''
    name_seq_dict = {}
    for record in SeqIO.parse(args.fasta, 'fasta'):
        name_seq_dict[record.id] = [str(record.seq)]
    
    aln_length = len(str(record.seq))

    ''' Prepare dictionary with wanted clusters '''
    wanted_clusters = {}
    for row in args.features:
        row = row.strip()
        values = set()
        if row[0] == '>':
            cluster_name = row[1:]
        else:
            seqnames = row.split(' ')
            values.update(seqnames)
            wanted_clusters[cluster_name] = values

    '''  '''
    #seq_combinations = product(name_seq_dict) # Creates all the possible combinations of the dictionary keys (seq names)
    for step in range(aln_length - kmer_length - 1):
        seq_names = []
        dist_matrix = []
        region_aln = select_region(name_seq_dict, step, kmer_length)
        
        for protein_seq_name1 in name_seq_dict:
            current_row = []
            seq_names.append(protein_seq_name1)
            for protein_seq_name2 in name_seq_dict:
                score_value = pairwise_comparison(protein_seq_name1, protein_seq_name2, region_aln, simple_score_dict)
                current_row.append(score_value)
            dist_matrix.append(current_row)

        #print(dist_matrix)
        condensed_dist = scipy.spatial.distance.squareform(dist_matrix)
        linkage = scipy.cluster.hierarchy.linkage(condensed_dist, method='complete')
        clusters_dict = {i:set([seqname]) for i, seqname in enumerate(seq_names)} # start the dict with clusters
        #print(linkage)
        for new_cl in linkage:
            #print(clusters_dict[20])
            #print(int(new_cl[0]))
            clusters_dict[max(clusters_dict)+1] = clusters_dict[int(new_cl[0])] | clusters_dict[int(new_cl[1])] # overlap two exited clusters
        
        for cluster_name_real, cluster_real in clusters_dict.items():
            for cluster_name_wanted, cluster_wanted in wanted_clusters.items():
                #print('%s\n----------------\n' % (cluster_wanted))
                if cluster_real == cluster_wanted:
                    if len(cluster_wanted) == 1:
                        if clusters_dict[max(clusters_dict)] - cluster_wanted not in clusters_dict.values():
                            break
                    outp.write('@%s\n' % (cluster_name_wanted))
                    for seq in region_aln:
                        outp.write('%s\t%s\n' % (seq, region_aln[seq]))

if __name__ == '__main__':
    main()        

