#!/usr/bin/env python3

import argparse
from collections import Counter, defaultdict
import itertools

def protein_annotation(db, db_type):
    prot_annot = defaultdict(str)
    for row in db:
        row = row.split('\t')
        if db_type == 'eggnog':
            if row[5]:
                prot_annot[row[0]] = [row[5], row[4]] # updated
        elif db_type == 'diamond':
            annot = row[17]
            annot = annot.split(maxsplit=1)[1]
            annot = annot.rsplit('[', maxsplit=1)[0].strip()
            organism = row[16] # new
            prot_annot[row[0]] = [annot, organism] # updated
        else:
            raise ValueError('wrong db type')
    return prot_annot

def keep_only_characterized(pg, db_dict):
    is_characterized = []
    for protein in pg:
        if protein in db_dict and ('uncharacterized' in db_dict[protein][0] or 'hypothetical' in db_dict[protein][0]): # updated
            is_characterized.append(False)
        else:
            is_characterized.append(True)

    if not any(is_characterized):
        return pg
    return list(itertools.compress(pg, is_characterized))

def prot_group_annot(pg, db_dict):
    prot_annots = []
    organisms = set()
    pg_new = keep_only_characterized(pg, db_dict)
    for protein in pg_new:
        if protein in db_dict:
            annot, organism = db_dict[protein]
            prot_annots.append(annot)
            organisms.add(organism)
    #print(prot_annots)
    return Counter(prot_annots), organisms

def get_most_common(annot_count):
    if annot_count:
        most_common_annotation = annot_count.most_common(1)[0][0]
    else:
        return '*'
    return most_common_annotation

def main():
    parser = argparse.ArgumentParser(description='Annotate protein groups')
    parser.add_argument('-p', type=argparse.FileType(),  help='protein groups file')
    parser.add_argument('-e', '--eggnog', type=argparse.FileType(), help='eggnog annotation of the predicted protein groups')
    parser.add_argument('-d', '--diamond', type=argparse.FileType(), help='diamond annotation of the predicted protein groups')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), help='the name of the output file')

    args = parser.parse_args()
    outp = args.output
    outp.write('protein_group\teggnog_annot\tdiamond_annot\tannotation\teggnog_proteins_annot\tdiamond_proteins_annot\t'
        'organism_eggnog\torganism_diamond\n')

    eggnog_db = protein_annotation(args.eggnog, 'eggnog')
    diamond_db = protein_annotation(args.diamond, 'diamond')
    
    #print(eggnog_db)

    for protein_group in args.p:
        if protein_group.startswith('Protein IDs'):
            continue
        protein_group = protein_group.split('\t')[0].split(';')
        eggnog_annot_count, eggnog_organisms = prot_group_annot(protein_group, eggnog_db)
        diamond_annot_count, diamond_organisms = prot_group_annot(protein_group, diamond_db)
        eggnog_whole_annot = ['%s:%i' % (name, count) for name, count in eggnog_annot_count.items()]
        diamond_whole_annot = ['%s:%i' % (name, count) for name, count in diamond_annot_count.items()]

        eggnog_most_common = get_most_common(eggnog_annot_count)
        diamond_most_common = get_most_common(diamond_annot_count)
        
        enriched_annot = eggnog_most_common
        if enriched_annot == '*':
            enriched_annot = diamond_most_common

        outp.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (';'.join(protein_group), eggnog_most_common, diamond_most_common, enriched_annot, 
            ';'.join(eggnog_whole_annot), ';'.join(diamond_whole_annot), ';'.join(eggnog_organisms), ';'.join(diamond_organisms)))

if __name__ == '__main__':
    main()
    

