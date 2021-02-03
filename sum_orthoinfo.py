#!/usr/bin/env python3

import argparse

def main():
    parser = argparse.ArgumentParser(description='To gather information regarding found orthologues and paralogues as well as type of Hsp70 (constitutive or inducible)')
    parser.add_argument('--outp_fish', type=argparse.FileType(), help='output tsv file from protein_fishing.py')
    parser.add_argument('--outp_orthofinder', type=argparse.FileType(), help='tsv output file from OrthoFinder (in the folder Orthogroups)')
    parser.add_argument('--outp_exonerate', type=argparse.FileType(), help='output file in sugar format after exonerate')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), help='generated output file')

    args = parser.parse_args()
    outp = args.output

    
