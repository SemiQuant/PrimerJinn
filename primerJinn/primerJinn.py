#!/bin/python
import os
import argparse

def main():
    parser = argparse.ArgumentParser(description='PrimerJinn: A tool for primer design and in silico PCR')
    parser.add_argument('--output_fasta', action='store_true', help='Output amplicon sequences in FASTA format', default=False)
    parser.add_argument('--exclude_primers', action='store_true', help='Exclude primer sequences from FASTA output', default=False)
    args = parser.parse_args()

    print("Welcome to primerJinn")
    print("")
    
    cmd_args = ""
    if args.output_fasta:
        cmd_args += " --output_fasta"
    if args.exclude_primers:
        cmd_args += " --exclude_primers"
    
    os.system(f'getMultiPrimerSet{cmd_args}')
    print("")
    os.system(f'PCRinSilico{cmd_args}')

if __name__ == '__main__':
    main()