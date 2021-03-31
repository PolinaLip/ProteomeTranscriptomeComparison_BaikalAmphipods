#!/usr/bin/env python3

import argparse

def make_info_dict(file_after_fish, species_name):
    protein_info_dict = {}
    for protein in file_after_fish:
        protein = protein.split('\t')
        protein_name = protein[0]
        pg_name = protein[2]
        annotation = protein[4]
        n_peptides = len(protein[1].split(';'))
        protein_info_dict[protein_name] = [pg_name, annotation, species_name, n_peptides]
    
    return protein_info_dict

def make_exoner_dict(file_after_exonerate, species_name):
    exonerate_result = {}
    for pr in file_after_exonerate:
        pr = pr.split(' ')
        if pr[0] == 'sugar:':
            protein_name = pr[1].rsplit(';', 1)
            if species_name in protein_name[1]:
                exonerate_result[protein_name[0]] = pr[5].split(';')[1]
    return exonerate_result

def find_species_column(column_names_list, species_name_pattern):
    index_ = [i for i, s in enumerate(column_names_list) if species_name_pattern in s]
    return index_

def collect_all_info(protein_list_orthofinder, exonerate_dict, info_dict, orthogroup_name, output_file):
    for protein in protein_list_orthofinder:
        #print(protein)
        protein_type = '*'
        if protein.split(';species')[0] in exonerate_dict:
            protein_type = exonerate_dict[protein.split(';species')[0]]
        protein = protein.split(';')[0]
        if protein in info_dict:
            output_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%i\n' % (protein, info_dict[protein][0], info_dict[protein][2], orthogroup_name, info_dict[protein][1], protein_type, info_dict[protein][3]))

def main():
    parser = argparse.ArgumentParser(description='To gather information regarding found orthologues and paralogues as well as type of Hsp70 (constitutive or inducible)')
    parser.add_argument('--outp_fish_ecy', type=argparse.FileType(), help='ecy output tsv file from protein_fishing.py')
    parser.add_argument('--outp_fish_eve', type=argparse.FileType(), help='eve output tsv file from protein_fishing.py')
    parser.add_argument('--outp_fish_gla', type=argparse.FileType(), help='gla output tsv file from protein_fishing.py')
    parser.add_argument('--outp_ortho', type=argparse.FileType(), help='tsv output file from OrthoFinder (in the folder Orthogroups)')
    parser.add_argument('--outp_exonerate', type=argparse.FileType(), help='output file in sugar format after exonerate')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), help='generated output file')

    args = parser.parse_args()
    outp = args.output
    outp.write('protein\tprotein_group\tspecies\torthogroup\tannotation\thsp70_type\tn_peptides\n')

    ecy_info_dict = make_info_dict(args.outp_fish_ecy, 'Ecy')
    eve_info_dict = make_info_dict(args.outp_fish_eve, 'Eve')
    gla_info_dict = make_info_dict(args.outp_fish_gla, 'Gla')
    
    exonerate_list = []
    for row in args.outp_exonerate:
        exonerate_list.append(row)
    
    ecy_exo_dict = make_exoner_dict(exonerate_list, 'Ecy')
    eve_exo_dict = make_exoner_dict(exonerate_list, 'Eve')
    gla_exo_dict = make_exoner_dict(exonerate_list, 'Gla')
    
    orthg_name = 0
    for orthogroup in args.outp_ortho:
        orthogroup = orthogroup.split('\t')
        if orthogroup[0] == '# Species':
            ecy_column = find_species_column(orthogroup, 'ecy')
            eve_column = find_species_column(orthogroup, 'eve')
            gla_column = find_species_column(orthogroup, 'gla')
            continue
        orthg_name += 1
        proteins_ecy = orthogroup[ecy_column[0]].strip().split(',')
        collect_all_info(proteins_ecy, ecy_exo_dict, ecy_info_dict, orthg_name, outp)
        proteins_eve = orthogroup[eve_column[0]].strip().split(',')
        collect_all_info(proteins_eve, eve_exo_dict, eve_info_dict, orthg_name, outp)
        proteins_gla = orthogroup[gla_column[0]].strip().split(',')
        collect_all_info(proteins_gla, gla_exo_dict, gla_info_dict, orthg_name, outp)

if __name__ == '__main__':
    main()
