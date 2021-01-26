#!/usr/bin/env python3

import argparse
from collections import defaultdict
from Bio import SeqIO

def align_peptides(protein_group_name, protein_seq_dict, peptide_seq_dict, peptide_ids):
    protein_group_name = protein_group_name.split(';')
    peptide_ids = peptide_ids.split(';') 
    protein_and_its_peptides = defaultdict(list)
    for protein in protein_group_name:
        for pep_id in peptide_ids:
            if peptide_seq_dict[pep_id] in protein_seq_dict[protein]:
                protein_and_its_peptides[protein].append(peptide_seq_dict[pep_id]) # save protein name (i.e., contig name) and the list of peptides aligned to it to the dict
    return(protein_and_its_peptides)

def main():
    parser = argparse.ArgumentParser(description='To extract target proteins, their sequences and observed peptides')
    parser.add_argument('--prot_groups', type=argparse.FileType(), help='proteinGroups.txt file')
    parser.add_argument('--evidence', type=argparse.FileType(), help='evidence.txt file')
    parser.add_argument('--annot', type=argparse.FileType(), help='proteinGroups annotation file')
    parser.add_argument('--proteins', type=argparse.FileType(), help='fasta file with protein database used in MQ run')
    parser.add_argument('--target', help='key words to find the target protein (comma separator)')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), help='the name of the output file')

    args = parser.parse_args()
    outp = args.outp
    prot_intr = args.target.split(',') # names of the protein of interest
    
    ''' Choose protein groups contained key words in annotation '''
    prot_group_annot = {}
    for pg in args.annot:
        annotation = pg.strip().split('\t')
        all_eggnog = annotation[5]
        all_diamond = annotation[6]
        if any(name in all_eggnog or name in all_diamond for name in prot_intr):
            prot_group_annot[annotation[0]] = [all_eggnog, all_diamond] # save protein groups and their annotation to the dict
    
    ''' Find peptides ids for the target proteins to connect them with evidence file '''
    for protein_group in args.prot_groups:
        protein_group = protein_group.strip().split('\t')
        if protein_group[0] in prot_group_annot:
            prot_group_annot[protein_group[0]].append(protein_group[102]) # protein_group[102] - a string with peptide ids through semicolon
    
    ''' Make a dict with peptide ids as keys and peptide sequences as values '''
    peptide_id_seq = {}
    evi = args.evidence
    next(evi)
    for peptide in evi:
        peptide = peptide.strip().split('\t')
        peptide_id_seq[peptide[82]] = peptide[0] # peptide[82] - peptide id

    ''' Extract target sequences '''
    protein_seq = {}
    for record in SeqIO.parse(args.proteins, "fasta"):
        for prot_group in prot_group_annot:
            proteins = prot_group.split(';')
            if any(protein == record.id for protein in proteins):
                protein_seq[record.id] = [record.seq]

    protein_group_number = 0
    for pr_group in prot_group_annot:
        protein_group_number += 1
        proteins_aligned_peptides = align_peptides(pr_group, protein_seq, peptide_id_seq, prot_group_annot[pr_group][2])
        

