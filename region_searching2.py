#!/usr/bin/env python3

import argparse
import numpy as np
from Bio import SeqIO
import biotite.sequence.align # to get blosum62 matrix for distance calculation
#from itertools import product # to perform pairwise comparison of keys from dictionary
import scipy
from scipy import spatial, cluster

''' To get distance value for the pair of aas '''
def get_distance(similarities, i, j):
    s_max = (similarities[i,i] + similarities[j,j]) / 2
    return s_max - similarities[i,j]

''' To select regions inside all sequences (dictionary of seq names and aligned regions) '''
def select_region(dict_w_seqs, index, klength):
    region_dict = {}
    for seq in dict_w_seqs:
        region_dict[seq] = dict_w_seqs[seq][0][index:index + klength]
    return region_dict

''' To find distance between two given sequences '''
def pairwise_comparison(seq1, seq2, score_dict):
    score = 0
    for i in range(len(seq1) - 1):
        letter_seq1 = seq1[i]
        letter_seq2 = seq2[i]
        score += score_dict[(letter_seq1, letter_seq2)]
    return score

def main():
    parser = argparse.ArgumentParser(
        description='To look for the regions in the sequences that determine the inducibility of proteins (ex. Hsp/Hsc70)')
    parser.add_argument('--fasta', type=argparse.FileType(),
        help='fasta file with aligned sequences of interest')
    parser.add_argument('--features', type=argparse.FileType(),
        help='file with wanted topology of sequences (based on features)')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'),
        help='output file with regions and associated features')
    parser.add_argument('-k', type=int, help='kmer length (= 6, by default)', default=6)
    parser.add_argument('--partial', help='Output partial matches', action='store_true')
    args = parser.parse_args()

    outp = args.output
    kmer_length = args.k

    ''' Prepare blosum62 dictionary or dictionary with simple scoring system (penalty = 1 if two letters are different) '''
    #blosum62_mtrx = substitution_matrices.load('BLOSUM62')
    #aas = blosum62_mtrx.alphabet.replace('*', '-')
    blosum62_matrix = biotite.sequence.align.SubstitutionMatrix.std_protein_matrix() # obtain blosum62
    aas = list(blosum62_matrix.get_alphabet1()) # only aas in the right order
    aas.append('-')
    similarities = blosum62_matrix.score_matrix() # only values of similarities
    distances = np.zeros(similarities.shape) 
    for i in range(distances.shape[0]):
        for j in range(distances.shape[1]):
            distances[i,j] = get_distance(similarities, i, j)
    
    # add penalty for the indel (the maximum of the penalty from blosum62 matrix)
    dist_with_gap_penalty = np.empty((0, len(distances[0]) + 1), int)
    for row in distances:
        row = list(np.append(row, distances.max()))
        dist_with_gap_penalty = np.append(dist_with_gap_penalty, np.array([row]), axis=0)

    indel_row = np.append(np.array([[distances.max()] * len(distances[0])]), 0)
    dist_with_gap_penalty = np.append(dist_with_gap_penalty, np.array([indel_row]), axis=0)

    blosum62_dict = {}
    for index_row, aa in enumerate(aas):
        for index_col, value in enumerate(dist_with_gap_penalty[index_row]):
           blosum62_dict[(aa, aas[index_col])] = value

    simple_score_dict = {}
    for letter1 in aas:
        for letter2 in aas:
            simple_score_dict[letter1, letter2] = int(letter1 != letter2)

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
        seq_names = list(name_seq_dict)
        dist_matrix = []
        region_aln = select_region(name_seq_dict, step, kmer_length)

        unique_seqs = list(set(region_aln.values()))
        unique_seqs_ixs = { seq: i for i, seq in enumerate(unique_seqs) }
        clusters_dict = { i: set() for i in range(len(unique_seqs)) }
        for name, seq in region_aln.items():
            clusters_dict[unique_seqs_ixs[seq]].add(name)

        for seq1 in unique_seqs:
            current_row = []
            for seq2 in unique_seqs:
                score_value = pairwise_comparison(seq1, seq2, blosum62_dict)
                current_row.append(score_value)
            dist_matrix.append(current_row)

        #print(dist_matrix)
        condensed_dist = scipy.spatial.distance.squareform(dist_matrix)
        #print(condensed_dist)
        if len(condensed_dist) == 0:
            continue
        linkage = scipy.cluster.hierarchy.linkage(condensed_dist, method='complete')
        #print(linkage)
        for new_cl in linkage:
            #print(clusters_dict[20])
            #print(int(new_cl[0]))
            clusters_dict[max(clusters_dict)+1] = clusters_dict[int(new_cl[0])] | clusters_dict[int(new_cl[1])]
            # overlap two existing clusters

        for cluster_name_real, cluster_real in clusters_dict.items():
            has_full_match = False
            for cluster_name_wanted, cluster_wanted in wanted_clusters.items():
                #print('%s\n----------------\n' % (cluster_wanted))
                if cluster_real == cluster_wanted:
                    if len(cluster_wanted) == 1:
                        if clusters_dict[max(clusters_dict)] - cluster_wanted not in clusters_dict.values():
                            break
                    outp.write('@%s, pos = %d, full\n' % (cluster_name_wanted, step + 1))
                    write_clusters(cluster_real, cluster_wanted, region_aln, outp)
                    #write_linkage(linkage, clusters_dict, outp)
                    has_full_match = True
                    break

            if has_full_match or not args.partial:
                continue
            for cluster_name_wanted, cluster_wanted in wanted_clusters.items():
                if len(cluster_real) > 1 and len(cluster_wanted) > 1 and len(cluster_real ^ cluster_wanted) == 1:
                    outp.write('@%s, pos = %d, partial\n' % (cluster_name_wanted, step + 1))
                    write_clusters(cluster_real, cluster_wanted, region_aln, outp)
                    #write_linkage(linkage, clusters_dict, outp)
                    break


def write_clusters(cluster_real, cluster_wanted, region_aln, outp):
    outp.write('# observed: %s\n' % ' '.join(cluster_real))
    outp.write('# wanted:   %s\n' % ' '.join(cluster_wanted))
    for seq in region_aln:
        outp.write('%-130s  %s\n' % (seq, region_aln[seq]))


def write_linkage(linkage, clusters_dict, outp):
    n = linkage.shape[0] + 1
    for i in range(n):
        outp.write('# Cluster %2d:           %s\n' % (i, ' '.join(clusters_dict[i])))
    for i in range(n - 1):
        outp.write('# Cluster %2d (%2d + %2d): %s\n' % (n + i, linkage[i, 0], linkage[i, 1],
            ' '.join(clusters_dict[n + i])))


if __name__ == '__main__':
    main()

