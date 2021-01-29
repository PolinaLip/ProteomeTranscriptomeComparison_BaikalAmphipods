#!/usr/bin/env python3

import argparse
from collections import defaultdict
from Bio import SeqIO

''' To find peptides aligned to the protein '''
def align_peptides(protein_group_name, protein_seq_dict, peptide_seq_dict, peptide_ids):
    protein_group_name = protein_group_name.split(';')
    peptide_ids = peptide_ids.split(';') 
    protein_and_its_peptides = defaultdict(set)
    for protein in protein_group_name:
        if 'REV__' in protein:
            continue
        for pep_id in peptide_ids:
            #print(protein, pep_id, protein_seq_dict[protein][0], peptide_seq_dict[pep_id], sep='\t')
            if peptide_seq_dict[pep_id] in protein_seq_dict[protein][0]:
                protein_and_its_peptides[protein].add(peptide_seq_dict[pep_id]) # save protein name (i.e., contig name) and the list of peptides aligned to it to the dict
    
    return protein_and_its_peptides 

''' To choose the best protein from the protein group - based on the follow: 
    1. The protein with the maximum number of peptides; 
    2. Sometimes there are several proteins with the same set of peptides. In this case, if proteins are with different lengths 
    the protein with the minimal sequence length will be choosen; 
    3. If proteins are with the same lengths (but some differences in the sequence) 
    the protein with the investigated species name will be choosen (we can rely on transcriptomic assembly) '''
def choose_the_best(aligned_peptides_dict, protein_seq_dict, species_name):
    the_best_proteins = []
    max_len = 0
    for pr in aligned_peptides_dict:
        if len(aligned_peptides_dict[pr]) > max_len:
            max_len = len(aligned_peptides_dict[pr])
    for pr in aligned_peptides_dict:
        if len(aligned_peptides_dict[pr]) == max_len:
            the_best_proteins.append(pr)

    if len(the_best_proteins) == 1:
        return the_best_proteins
    else:
        best_protein = []
        if all(len(protein_seq_dict[the_best_proteins[0]][0]) == len(protein_seq_dict[pr][0]) for pr in the_best_proteins[1:]):
            for protein_name in the_best_proteins:
                if species_name in protein_name:
                    best_protein.append(protein_name)
        else:
            min_len = len(protein_seq_dict[the_best_proteins[0]][0])
            for pr_name in the_best_proteins[1:]:
                if len(protein_seq_dict[pr_name][0]) < min_len:
                    min_len = len(protein_seq_dict[pr_name][0])
            for pr_name in the_best_proteins:
                if len(protein_seq_dict[pr_name][0]) == min_len:
                    best_protein.append(pr_name)
        return best_protein

''' To find the indices of the substring (peptides) in the string (protein) '''
def locate_peptide(protein, protein_pept_dict, protein_seq_dict):
    peptide_locations = []
    for peptide in protein_pept_dict[protein]:
        protein_seq = protein_seq_dict[protein]
        peptide_start = protein_seq[0].find(peptide)
        if peptide_start == -1:
            print('Peptide %s is not in the %s' % (peptide, protein))
        else:
            peptide_locations.append((peptide_start, peptide))
    
    return peptide_locations

''' To compose strings with peptides which will be used for the protein-peptides alignment '''
def peptide2print(pep_indices):
    to_print = []
    for pep in pep_indices:
        if not to_print:
            align_str = '-' * pep[0] + pep[1]
            to_print.append(align_str)
            continue
        intersect = True
        for i, string in enumerate(to_print):
            if len(string) < pep[0]:
                to_print[i] += '-' * (pep[0] - len(string)) + pep[1]
                intersect = False
                break
        if intersect:
            align_str = '-' * pep[0] + pep[1]
            to_print.append(align_str)

    return to_print

''' To launch locate_peptide function -> sorting of tuples (index, peptide sequence) -> launch peptide2print '''
def render_peptides(protein, proteins_aligned_peptides_dict, protein_seq_dict):
    peptides_indices = locate_peptide(protein, proteins_aligned_peptides_dict, protein_seq_dict)
    peptides_indices.sort()
    to_print = peptide2print(peptides_indices)
    
    return to_print

def main():
    parser = argparse.ArgumentParser(description='To extract target proteins, their sequences and observed peptides')
    parser.add_argument('--prot_groups', type=argparse.FileType(), help='proteinGroups.txt file')
    parser.add_argument('--evidence', type=argparse.FileType(), help='evidence.txt file')
    parser.add_argument('--annot', type=argparse.FileType(), help='proteinGroups annotation file')
    parser.add_argument('--proteins', type=argparse.FileType(), help='fasta file with protein database used in MQ run')
    parser.add_argument('--target', help='key words to find the target protein (comma separator); case-insensitive')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), help='the name of the output file')
    parser.add_argument('--alignment', type=argparse.FileType('w'), help='file with proteins and peptides aligned to them')
    parser.add_argument('--alignment_best', type=argparse.FileType('w'), help='file with the best proteins (from protein group) and peptides aligned to them')
    parser.add_argument('--species', help='the species name; in the format used in the end of contig (protein) names (ex.: Ecy)')
    parser.add_argument('--exonerate_outp', type=argparse.FileType(), help='exonerate output file only with sugar format')
    parser.add_argument('--swiss_db', type=argparse.FileType(), help='SWISS-Prot database with target sequences (needed to extracting annotations)')

    args = parser.parse_args()
    outp1 = args.output
    outp1.write('protein\tpeptides\tprotein_group\tpg_number\tfinal_annotation\teggnog_annotations\tdiamond_annotations\tbest_protein\thomology_annot\tproteins_found_by_homology\n')
    outp2 = args.alignment
    outp3 = args.alignment_best
    prot_intr = args.target.upper().split(',') # list of alternative names of the protein of interest
    
    ''' Prepare swiss database and the output data from exonerate '''
    swiss_db = {}
    for s in args.swiss_db:
        if s[0] == '>':
            pr_annotation = s.split(' ', 1)[1].split('OS=')[0]
            pr_name = s[1:].split(' ')[0]
            swiss_db[pr_name] = pr_annotation # dict with fasta header (cutted by the first space) as a key and annotation as a value

    found_by_homology = {}
    for p in args.exonerate_outp:
        p = p.split(' ')
        if p[0] == 'sugar:':
            db_target = p[5]
            gene_id = db_target.split('|')[1]
            contig = p[1]
            found_by_homology[contig] = '%s:%s' % (gene_id, swiss_db[db_target]) # dict with contig (protein) name as a key and annotation as a value 

    ''' Choose protein groups contained key words in annotation or protein groups proteins in which were found by homology with the target database'''
    prot_group_annot = {}
    for pg in args.annot:
        annotation = pg.strip().split('\t')
        pg_name = annotation[0].split(';')
        found_proteins = []
        annot_info = None
        for pg_protein in pg_name:
            if pg_protein in found_by_homology:
                found_proteins.append(pg_protein)
                annot_info = found_by_homology[pg_protein]
        if found_proteins:
            prot_group_annot[annotation[0]] = [annotation[4], annotation[5], annotation[3], ';'.join(found_proteins), annot_info]
            continue
        if annotation[3] == '*':
            continue
        all_eggnog = annotation[4]
        all_diamond = annotation[5]
        final_annot = annotation[3]
        if any(name in all_eggnog.upper() or name in all_diamond.upper() for name in prot_intr):
            prot_group_annot[annotation[0]] = [all_eggnog, all_diamond, final_annot, '*', '*'] # save protein groups and their annotation to the dict
    
    ''' Find peptides ids for the target proteins to connect them with evidence file '''
    for protein_group in args.prot_groups:  
        protein_group = protein_group.strip().split('\t')
        if protein_group[0] == 'Protein IDs':
            column_with_ids = protein_group.index('Evidence IDs')
        if protein_group[0] in prot_group_annot:
            prot_group_annot[protein_group[0]].append(protein_group[column_with_ids]) # protein_group[column_with_ids] - a string with peptide ids through semicolon
    
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
                protein_seq[record.id] = [str(record.seq)]
    
    ''' Combine all info about protein together and create three outputs: 
    1. proteins and info, 2. proteins with aligned peptides, 3. best proteins with aligned peptides '''
    protein_group_number = 0
    for pr_group in prot_group_annot:
        protein_group_number += 1
        proteins_aligned_peptides = align_peptides(pr_group, protein_seq, peptide_id_seq, prot_group_annot[pr_group][5])
        #print(proteins_aligned_peptides)
        the_best_proteins_list = choose_the_best(proteins_aligned_peptides, protein_seq, args.species)

        for prot in proteins_aligned_peptides:
            best = 'no'
            if prot in the_best_proteins_list:
                best = 'yes'
                to_print_best = render_peptides(prot, proteins_aligned_peptides, protein_seq)
                protein_sequence_best = protein_seq[prot]
                outp3.write('>%s;protein_group=%i\n' % (prot, protein_group_number))
                outp3.write('%s\n' % (protein_sequence_best[0]))
                for s in to_print_best:
                    outp3.write('%s\n' % (s))
                outp3.write('\n')

            outp1.write('%s\t%s\t%s\t%i\t%s\t%s\t%s\t%s\t%s\t%s\n' % (prot, ';'.join(proteins_aligned_peptides[prot]), pr_group, protein_group_number, prot_group_annot[pr_group][2], prot_group_annot[pr_group][0], prot_group_annot[pr_group][1], best, prot_group_annot[pr_group][4], prot_group_annot[pr_group][3]))

            to_print = render_peptides(prot, proteins_aligned_peptides, protein_seq)
            
            protein_sequence = protein_seq[prot]
            
            outp2.write('>%s;protein_group=%i\n' % (prot, protein_group_number))
            outp2.write('%s\n' % (protein_sequence[0]))
            for s in to_print:
                outp2.write('%s\n' % (s))
            outp2.write('\n')

if __name__ == '__main__':
    main()
