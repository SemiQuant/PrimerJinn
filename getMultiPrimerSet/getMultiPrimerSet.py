#!/bin/python
from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import primer3
import pandas as pd
import itertools
import argparse
import numpy as np
import concurrent.futures
import sys
import random

parser = argparse.ArgumentParser(description='PCR primer design')

from getMultiPrimerSet.getMultiPrimerSet import main

# Define input arguments
# parser.add_argument('--cpus', type=int, help='number of threads (even though it says cpus)', default=1000)
parser.add_argument('--region_file', type=str, help='The path to the primer design region file. Two columns, start position and end position (1-based). tsv or xlxs', required=True)
parser.add_argument('--input_file', type=str, help='The path to the input FASTA file.', required=True)
parser.add_argument('--target_tm', type=float, help='The desired melting temperature (Tm) for the primers.', default=65)
parser.add_argument('--primer_len', type=int, help='The desired length of the primers.', default=20)
parser.add_argument('--product_size_min', type=int, help='The desired min size for the PCR product. Default is 400.', default=400)
parser.add_argument('--product_size_max', type=int, help='The desired max size for the PCR product. Default is 800.', default=800)
parser.add_argument('--ret', type=int, help='The maximum number of primer pairs to return.', default=100)
parser.add_argument('--Q5', action='store_true', help='Whether to use Q5 settings for primer3.', default=True)
parser.add_argument('--background', type=str, help='The path to the mispriming library FASTA file.', default='')
parser.add_argument('--output', type=str, help='Output name.', default='MultiPlexPrimerSet')
parser.add_argument('--eval', type=int, help='The maximum number of primer sets to evaluate.', default=10000)
args = parser.parse_args()

product_size_range = (args.product_size_min, args.product_size_max)

def design_primers(input_file, region_start, region_end, target_tm=65, primer_len=20, product_size_range=(400, 800), name='', ret=100, Q5=True, background=''):
    """
    This function takes a FASTA file, start and end positions of a target region, and desired melting temperature
    and primer length, and returns the best primer set for that region using Biopython and primer3.

    Args:
    input_file (str): The path to the input FASTA file.
    region_start (int): The start position of the target region (1-based).
    region_end (int): The end position of the target region (1-based).
    target_tm (float): The desired melting temperature (Tm) for the primers.
    primer_len (int): The desired length of the primers.
    product_size_range (tuple): The desired size range for the PCR product.
    name (str): The name of the primer set.
    ret (int): The maximum number of primer pairs to return.
    Q5 (bool): Whether to use Q5 settings for primer3.
    background (str): The name of a mispriming library to use for primer3.

    Returns:
    list: A list of dicts containing information about the primer pairs.
    """

    # Parse the input FASTA file and extract the target region sequence.
    print("Picking initial primers for " + name)
    if region_start - 10000 > 0:
        target_seq = input_file.seq[region_start-10000:region_end+10000]
        seq_target = (10000-1, region_end-region_start)
    else:
        target_seq = input_file.seq[0:end]
        seq_target = (region_start-1, region_end-region_start)
    
    # target_seq = input_file.seq
    # seq_target = (region_start-1, region_end-region_start+2)

    # Set up the primer3 input parameters.
    input_params = {
    'SEQUENCE_ID': record.id,
    'SEQUENCE_TEMPLATE': str(target_seq),
    'SEQUENCE_TARGET': seq_target,
    'PRIMER_OPT_SIZE': primer_len,
    'PRIMER_MIN_SIZE' : primer_len - 10,
    'PRIMER_MAX_SIZE' : primer_len + 20,
    'PRIMER_PICK_INTERNAL_OLIGO': 0,
    }

    if background != '':
        input_params.update({'PRIMER_MISPRIMING_LIBRARY': background,})

    global_args = {
    'PRIMER_MIN_GC': 40,
    'PRIMER_MAX_GC': 80,
    'PRIMER_NUM_RETURN': ret,
    'PRIMER_EXPLAIN_FLAG': 1,
    'PRIMER_MIN_TM': target_tm-4,
    'PRIMER_OPT_TM': target_tm,
    'PRIMER_MAX_TM': target_tm+4,
    'PRIMER_PRODUCT_SIZE_RANGE': product_size_range,
    }
        
    if Q5:
        global_args.update({
        'PRIMER_SALT_DIVALENT': 2,
        'PRIMER_SALT_MONOVALENT': 80,
        'PRIMER_DNA_CONC': 5800,
        'PRIMER_TM_FORMULA': 1,
        'PRIMER_SALT_CORRECTIONS': 2,
        })

    # Run primer3 to design the primers.
    primers = primer3.bindings.designPrimers(input_params, global_args)
    try:
        # Check if any primer pairs were returned.
        if 'PRIMER_PAIR_NUM_RETURNED' not in primers or primers['PRIMER_PAIR_NUM_RETURNED'] == 0:
            # print("No primers found for: " + name)
            return []
        else:
            # Extract information about the primer pairs.
            primer_pairs = []
            for i in range(primers['PRIMER_PAIR_NUM_RETURNED']):
                forward_seq = primers[f'PRIMER_LEFT_{i}_SEQUENCE']
                reverse_seq = primers[f'PRIMER_RIGHT_{i}_SEQUENCE']
                forward_tm = primers[f'PRIMER_RIGHT_{i}_TM']
                reverse_tm = primers[f'PRIMER_LEFT_{i}_TM']
                product_size = primers[f'PRIMER_PAIR_{i}_PRODUCT_SIZE']

                # Check that the product size is within the specified range.
                # if product_size < int(product_size_range[0]) or product_size > int(product_size_range[1]):
                primer_pairs.append({
                    'Forward Primer': forward_seq,
                    'Reverse Primer': reverse_seq,
                    'Forward tm': forward_tm,
                    'Reverse tm': reverse_tm,
                    'Product Size': product_size,
                    'name' : name
                    })
    except Exception as e:
        print("An error occurred:", e)

    return primer_pairs


def is_valid_combination(combination, names):
    if set(entry['name'] for entry in combination) == names:
        return combination
    return []

def evaluate_combination(comb):
    # Create a list of primer pairs
    primer_pairs = [(p1, p2) for p1, p2 in itertools.combinations(comb, 2)]
    # Extract the primer sequences from the primer pairs
    primer_pairs_sequences = [(p1['Forward Primer'], p2['Reverse Primer']) for p1, p2 in primer_pairs]
    # Calculate heterodimer formation energy
    heterodimer_scores = []
    for p in primer_pairs_sequences:
        heterodimer_scores.append(primer3.calcHeterodimer(p[0], p[1], temp_c=args.target_tm).dg) #gibbs free energy
        
    score = {'size': product_size_range,
             'tmp': tm_range,
             'mean': np.mean(heterodimer_scores),
             'range': abs(np.min(heterodimer_scores)) - abs(np.max(heterodimer_scores))
            }
    
    return score, comb

# Read the file into a pandas DataFrame
file_path = args.region_file  # or 'file.xlsx'
if file_path.endswith('.tsv'):
    df = pd.read_csv(file_path, sep='\t')
else:
    df = pd.read_excel(file_path)

seen_names = set()
modified_names = {}

primers_all = []

record = SeqIO.read(args.input_file, "fasta")

# Iterate over each row
for index, row in df.iterrows():
    # Access the "start" and "end" values using the column names
    # name=str(row['name'])
    start = int(row['start'])
    end = int(row['end'])
    name = 'Target ' + str(start) + "-" + str(end)

    if end - start > int(product_size_range[1]):
        print("Region " + str(start) + ", ", str(end) + " larger than product_size_range!")
        break

    # Check if name is already in seen_names
    if name in seen_names:
        # If so, modify the name to make it unique
        if name in modified_names:
            modified_name = modified_names[name]
        else:
            modified_name = name + '_' + str(random.randint(1, 100))
            modified_names[name] = modified_name
        # Update name in row
        df.at[index, 'name'] = modified_name
        seen_names.add(modified_name)
        name = modified_name
    else:
        seen_names.add(name)

    # Call the design_primers function to design primers for the target region.
    primer_tmp = []
    primer_tmp = design_primers(
        input_file=record,
        region_start=start,
        region_end=end,
        target_tm=args.target_tm,
        primer_len=args.primer_len,
        product_size_range=product_size_range,
        name=name,
        ret=args.ret,
        Q5=args.Q5,
        background=args.background
        )

    if primer_tmp is not None and len(primer_tmp) > 0:
        primers_all.append(primer_tmp)
    else:
        print("No primer found for: " + name)


#Typically, a ΔG value below -9 kcal/mol is considered to indicate a high likelihood of dimer formation, while values between -6 and -9 kcal/mol indicate a moderate likelihood. Values above -6 kcal/mol are considered unlikely to result in dimer formation.
# Define initial placeholder score
best_score = {'size': float('inf'), 'tmp': float('inf'), 'range': float('inf'), 'mean': float('inf')}
best_comb = []
# get a set of unique names
flat_list = []
for sublist in primers_all:
    for item in sublist:
        flat_list.append(item)

primers_all = flat_list

names = set([primer['name'] for primer in primers_all])  # get unique names
if len(names) <= 1:
    print("not enough targets! " + ' '.join(map(str, names)))
    exit(1)
else:
    print("Finding best primer set for the following targets: " + ', '.join(map(str, names)))
    best_score = {'mean': float('inf'), 'range': float('inf'), 'tmp': float('inf'), 'size': float('inf')}
    best_comb = None
    valid_combinations = []
    pos = 0
    name_primers_dict = {}
    for name in names:
        name_primers_dict[name] = [primer for primer in primers_all if primer['name'] == name]  # group primers by name

    # This nees to be fixed
    combinations = []
    for i in range(args.eval):
        combination = []
        for name in names:
            name_primers = name_primers_dict[name]  # get primers for the current name
            combination.append(random.choice(name_primers))
        combinations.append(combination)

    print("Picking best set from: ", len(combinations), " combinations")
    # tm and size
    product_sizes = [d['Product Size'] for d in primers_all]
    tm_values = [d['Forward tm'] for d in primers_all] + [d['Reverse tm'] for d in primers_all]
    # product_size_range = max(product_sizes) - min(product_sizes)
    tm_range = max(tm_values) - min(tm_values)
    best_score = {'mean': float('inf'), 'range': float('inf'), 'tmp': float('inf')}
    best_comb = None

    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [executor.submit(evaluate_combination, comb) for comb in combinations]

        for future in concurrent.futures.as_completed(futures):
            score, comb = future.result()
            if best_score['mean'] > score['mean']:
                if best_score['range'] * 0.9 > score['range']:
                    best_score = score
                    best_comb = comb
                elif best_score['tmp'] > 4 or score['tmp'] > 4 and best_score['tmp'] > score['tmp']:
                    best_score = score
                    best_comb = comb

    print("")
    print("writing file: " + args.output + '.xlsx')
    # convert best_comb tuple to DataFrame object
    best_comb = pd.DataFrame(list(best_comb))
    best_comb.to_excel(args.output + '.xlsx', sheet_name='PrimerSet', index=False)

    