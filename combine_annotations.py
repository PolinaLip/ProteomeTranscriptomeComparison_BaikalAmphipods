#!/usr/bin/env python3

import argparse

def main():
    parser = argparse.ArgumentParser(description='combine DA and FA annotations, FA annotation is in the priority')
    parser.add_argument('-d', type=argparse.FileType(),  help='diamond annotation')
    parser.add_argument('-e', type=argparse.FileType(), help='eggnog annotation')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), help='the names of output file')

    args = parser.parse_args()
    outp = args.output
    outp.write('contig\tannot\tdiamond_annot\n')
    diamond = args.d.readlines()

    diamond_dict = {}
    for row in diamond:
        row = row.strip().split("\t")
        annot = row[19]
        if '<>' in annot:
            annot = annot.split('<')[1]
        diamond_dict[row[1]] = ' '.join(annot.rsplit('[', 1)[0].split(' ')[1:])

    eggnog_dict = {}
    for row in args.e:
        if row[0] != '#':
            row = row.split('\t')
            eggnog_dict[row[0]] = row[5]
    
    for row in diamond:
        contig_name = row.split('\t')[1]
        diamond_annot = diamond_dict[contig_name]
        eggnog_annot = None
        if contig_name in eggnog_dict:
            eggnog_annot = eggnog_dict[contig_name]
        
        if eggnog_annot:
            outp.write('%s\t%s\t%s\n' % (contig_name, eggnog_annot, diamond_annot))
        else:
            outp.write('%s\t%s\t%s\n' % (contig_name, diamond_annot, diamond_annot))

if __name__ == '__main__':
    main()

