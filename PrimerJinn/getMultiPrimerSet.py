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
import requests
from sklearn.cluster import AgglomerativeClustering
import tempfile
import subprocess
import shutil
import os
import openpyxl
import matplotlib.pyplot as plt
from collections import defaultdict
from tabulate import tabulate

parser = argparse.ArgumentParser(description='PCR primer design')
# Define input arguments
# parser.add_argument('--cpus', type=int, help='number of threads (even though it says cpus)', default=1000)
parser.add_argument('--region_file', type=str, help='The path to the primer design region file. Two columns, start position and end position (1-based). tsv or xlxs', required=True)
parser.add_argument('--input_file', type=str, help="Reference FASTA file", required = True)
parser.add_argument('--target_tm', type=float, help='The desired melting temperature (Tm) for the primers.', default=65)
parser.add_argument('--primer_len', type=int, help='The desired length of the primers.', default=20)
parser.add_argument('--product_size_min', type=int, help='The desired min size for the PCR product. Default is 400.', default=400)
parser.add_argument('--product_size_max', type=int, help='The desired max size for the PCR product. Default is 800.', default=800)
parser.add_argument('--ret', type=int, help='The maximum number of primer pairs to return.', default=100)
parser.add_argument('--Q5', action='store_true', help='Whether to use Q5 approximation settings for Tm calculations.', default=True)
parser.add_argument('--background', type=str, help='The path to the mispriming library FASTA file.', default='')
parser.add_argument('--no_self_background', action='store_true', help='If specified, primers are NOT checked for mispriming against the input genome as for large genomes this is very slow')
parser.add_argument('--output', type=str, help='Output name.', default='MultiPlexPrimerSet')
parser.add_argument('--eval', type=int, help='The maximum number of primer sets to evaluate.', default=10000)
parser.add_argument('--wiggle', type=int, help='Half the region around the optimal Tm', default=3)
parser.add_argument('--ill_adapt', action='store_true', help='Add Illumina partial adapters', default=False)
parser.add_argument('--clamp', type=int, help='Require GC clamp', default=0)
parser.add_argument('--poly', type=int, help='Maximum allowable length of a mononucleotide repeat (poly-X) in the primer sequence', default=3)
args = parser.parse_args()

product_size_range = (args.product_size_min, args.product_size_max)

def design_primers(input_file, region_start, region_end, target_tm=args.target_tm, primer_len=20, product_size_range=(400, 800), name='', ret=100, Q5=True, background='', wiggle=args.wiggle, clamp=args.clamp, poly=args.poly, no_self_background=False):
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

    print("Picking initial primers for " + name)

    # Initialize mispriming library with minimal content if empty
    if os.path.exists(background) and os.path.getsize(background) == 0:
        with open(background, "w") as fasta_file:
            fasta_file.write(">dummy_sequence\n")
            fasta_file.write("GCTAGCTAGCTAGCTAGCTAG\n")

    # Calculate flanking region size based on max product size
    flank_size = min(500, product_size_range[1])
    
    # Extract sequence with flanking regions, ensuring we don't go out of bounds
    seq_start = max(0, region_start - flank_size)
    seq_end = min(len(input_file.seq), region_end + flank_size)
    target_seq = input_file.seq[seq_start:seq_end]
    
    # Adjust target coordinates relative to extracted sequence
    target_start = region_start - seq_start
    target_length = region_end - region_start
    seq_target = (target_start, target_length)

    # Add non-target regions to mispriming library if requested
    if not no_self_background and len(input_file.seq) > len(target_seq):
        with open(background, "a") as fasta_file:
            if seq_start > 0:
                fasta_file.write(">upstream_sequence\n")
                fasta_file.write(str(input_file.seq[0:seq_start]) + "\n")
            if seq_end < len(input_file.seq):
                fasta_file.write(">downstream_sequence\n")
                fasta_file.write(str(input_file.seq[seq_end:]) + "\n")

    # Set up primer3 parameters
    input_params = {
        'SEQUENCE_ID': input_file.id,
        'SEQUENCE_TEMPLATE': str(target_seq),
        'SEQUENCE_TARGET': seq_target,
        'PRIMER_OPT_SIZE': primer_len,
        'PRIMER_MIN_SIZE': 10,
        'PRIMER_MAX_SIZE': 36,
        'PRIMER_PICK_INTERNAL_OLIGO': 0,
    }

    if not no_self_background:
        input_params['PRIMER_MISPRIMING_LIBRARY'] = background

    global_args = {
        'PRIMER_MIN_GC': 40,
        'PRIMER_MAX_GC': 80,
        'PRIMER_NUM_RETURN': ret,
        'PRIMER_EXPLAIN_FLAG': 1,
        'PRIMER_MIN_TM': target_tm-wiggle,
        'PRIMER_OPT_TM': target_tm,
        'PRIMER_MAX_TM': target_tm+wiggle,
        'PRIMER_PRODUCT_SIZE_RANGE': product_size_range,
        'PRIMER_INTERNAL_MAX_POLY_X': poly,
    }

    if Q5:
        global_args.update({
            'PRIMER_SALT_DIVALENT': 2,
            'PRIMER_SALT_MONOVALENT': 70,
            'PRIMER_DNA_CONC': 3300,
            'PRIMER_TM_FORMULA': 1,
            'PRIMER_SALT_CORRECTIONS': 2,
        })

    if clamp > 0:
        global_args.update({'PRIMER_GC_CLAMP': clamp})

    # Run primer3 with increasingly relaxed constraints if needed
    primers = primer3.bindings.design_primers(input_params, global_args)
    
    if 'PRIMER_PAIR_NUM_RETURNED' not in primers or primers['PRIMER_PAIR_NUM_RETURNED'] == 0:
        print(f"No primers found for: {name}; trying with relaxed constraints.")
        del global_args['PRIMER_INTERNAL_MAX_POLY_X']
        global_args.update({'PRIMER_PICK_ANYWAY': 1})
        primers = primer3.bindings.design_primers(input_params, global_args)
        
        if 'PRIMER_PAIR_NUM_RETURNED' not in primers or primers['PRIMER_PAIR_NUM_RETURNED'] == 0:
            if not no_self_background and 'PRIMER_MISPRIMING_LIBRARY' in input_params:
                del input_params['PRIMER_MISPRIMING_LIBRARY']
            global_args.update({
                'PRIMER_MIN_TM': target_tm - wiggle * 1.5,
                'PRIMER_MAX_TM': target_tm + wiggle * 1.5
            })
            primers = primer3.bindings.design_primers(input_params, global_args)

    try:
        if 'PRIMER_PAIR_NUM_RETURNED' not in primers or primers['PRIMER_PAIR_NUM_RETURNED'] == 0:
            return []

        # Extract and process primer pairs
        primer_pairs = []
        for i in range(primers['PRIMER_PAIR_NUM_RETURNED']):
            forward_seq = primers[f'PRIMER_LEFT_{i}_SEQUENCE']
            reverse_seq = primers[f'PRIMER_RIGHT_{i}_SEQUENCE']
            forward_tm = primers[f'PRIMER_LEFT_{i}_TM']
            reverse_tm = primers[f'PRIMER_RIGHT_{i}_TM']
            product_size = primers[f'PRIMER_PAIR_{i}_PRODUCT_SIZE']
            
            # Adjust binding sites to global coordinates
            start = primers[f'PRIMER_LEFT_{i}'][0] + seq_start
            end = primers[f'PRIMER_RIGHT_{i}'][0] + seq_start

            if args.ill_adapt:
                forward_seq = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT" + forward_seq
                reverse_seq = "GACTGGAGTTCAGACGTGTGCTCTTCCGATCT" + reverse_seq

            primer_pairs.append({
                'Forward Primer': forward_seq,
                'Reverse Primer': reverse_seq,
                'Forward tm': forward_tm,
                'Reverse tm': reverse_tm,
                'Product Size': product_size,
                'Binding Start': str(start),
                'Binding End': str(end),
                'name': name
            })

    except Exception as e:
        print(f"An error occurred while processing primers for {name}:", e)
        return []

    return primer_pairs

def q5_melting_temp(seq1, seq2="", salt=0.5):
    url = "https://tmapi.neb.com/tm/q5/%s/%s/%s?&fmt=short" % (salt, seq1, seq2)
    response = requests.get(url)
    json_string = response.json()
    tm1 = json_string["data"]["tm1"]
    # tm2 = json_string["data"]["tm2"]
    return tm1 #[tm1, tm2]

def is_valid_combination(combination, names):
    if set(entry['name'] for entry in combination) == names:
        return combination
    return []

def evaluate_combination(comb, Q5=args.Q5):
    # Create a list of primer pairs
    primer_pairs = [(p1, p2) for p1, p2 in itertools.combinations(comb, 2)]
    # Extract the primer sequences from the primer pairs
    primer_pairs_sequences = [(p1['Forward Primer'], p2['Reverse Primer']) for p1, p2 in primer_pairs]
    # Calculate heterodimer formation energy
    heterodimer_scores = []
    for p in primer_pairs_sequences:
        if Q5:
            heterodimer_scores.append(primer3.calc_heterodimer(p[0], p[1], temp_c=args.target_tm, dv_conc = 2, mv_conc = 70, dna_conc = 3300).dg) #gibbs free energy
        else:
            heterodimer_scores.append(primer3.calc_heterodimer(p[0], p[1], temp_c=args.target_tm).dg) #gibbs free energy
    product_size_weight=0.2
    heterodimer_score_weight =0.8
    
    # Calculate tm_range
    tm_values = []
    for primer in comb:
        tm_values.append(primer['Forward tm'])
        tm_values.append(primer['Reverse tm'])
    tm_range = max(tm_values) - min(tm_values)
    
    score = {'size': product_size_range,
             'tmp': tm_range,
             'mean': np.mean(heterodimer_scores),
             'range': abs(np.min(heterodimer_scores)) - abs(np.max(heterodimer_scores))
             # 'weighted_score': product_size_weight * (1 - product_size_range) + heterodimer_score_weight * (1 - abs(np.mean(heterodimer_scores)) / 100)
            }    
    return score, comb


# Define the features to consider for clustering
def get_features(primer, target_tm=args.target_tm):
    fwd_tm = primer['Forward tm']
    rev_tm = primer['Reverse tm']
    avg_tm = (fwd_tm + rev_tm) / 2
    product_size = primer['Product Size']
    tm_difference = abs(fwd_tm - rev_tm)
    dimer_prob = primer3.calc_heterodimer(primer['Forward Primer'], primer['Reverse Primer'], temp_c=target_tm, dv_conc = 2, mv_conc = 70, dna_conc = 3300).dg
    return (avg_tm, product_size, tm_difference, dimer_prob)
    

def cluster_primers(primers, target_tm=args.target_tm, distance_threshold=10):
    # Cluster the primers based on their features
    X = [get_features(primer) for primer in primers]
    clusterer = AgglomerativeClustering(n_clusters=None, metric='euclidean', linkage='ward', distance_threshold=distance_threshold)
    clusters = clusterer.fit_predict(X)

    # Select the largest cluster
    name_counts = {}
    for i, cluster in enumerate(clusters):
        name = primers[i]['name']
        if cluster in name_counts:
            name_counts[cluster].add(name)
        else:
            name_counts[cluster] = {name}

    most_unique_cluster = max(name_counts, key=lambda k: len(name_counts[k]))
    largest_cluster = [primers[i] for i, cluster in enumerate(clusters) if cluster == most_unique_cluster]

    # Add missing names from the next closest cluster until all names are included
    all_names = set(primer['name'] for primer in primers)
    selected_names = set(primer['name'] for primer in largest_cluster)
    missing_names = all_names - selected_names

    while missing_names:
        # Find the closest cluster that contains missing names
        closest_cluster_idx = None
        closest_distance = float('inf')
        for i, cluster in enumerate(clusters):
            if cluster == most_unique_cluster:
                continue
            cluster_names = set(primers[j]['name'] for j in range(len(primers)) if clusters[j] == cluster)
            common_names = cluster_names & missing_names
            if common_names:
                distance = np.linalg.norm(np.array(X[i]) - np.array(X[0]))
                if distance < closest_distance:
                    closest_cluster_idx = cluster
                    closest_distance = distance

        # Add the missing names to the largest cluster
        closest_cluster = [primers[i] for i, cluster in enumerate(clusters) if cluster == closest_cluster_idx]
        for primer in closest_cluster:
            if primer['name'] in missing_names:
                largest_cluster.append(primer)
                missing_names.remove(primer['name'])

    # Remove duplicates by name
    selected_primers = []
    selected_names = set()
    for primer in largest_cluster:
        name = primer['name']
        if name not in selected_names:
            selected_primers.append(primer)
            selected_names.add(name)
    return selected_primers


def plot_cluster(clusterer, primers):
    # Get features for each primer
    X = [get_features(primer) for primer in primers]
    # Convert to 2D for visualization using first two features
    X_2d = [[x[0], x[1]] for x in X]
    labels = clusterer.labels_
    
    # Create the plot
    fig, ax = plt.subplots()
    ax.scatter([x[0] for x in X_2d], [x[1] for x in X_2d], c=labels, cmap='viridis')
    ax.set_xlabel('Average Tm')
    ax.set_ylabel('Product Size')
    
    # Add labels for each point
    for i, primer in enumerate(primers):
        ax.annotate(primer['name'], (X_2d[i][0], X_2d[i][1]))
    
    plt.title('Primer Clusters')
    
    # Create a table of cluster assignments
    rows = []
    cluster_dict = defaultdict(list)
    for i, primer in enumerate(primers):
        cluster_dict[labels[i]].append(primer['name'])
    
    for cluster_id, primer_names in cluster_dict.items():
        rows.append([', '.join(primer_names), cluster_id])
    
    return (fig, tabulate(rows, headers=['Primer Name', 'Clusters']))


def main():
    #########################################
    # Read the file into a pandas DataFrame #
    file_path = args.region_file  # or 'file.xlsx'
    if file_path.endswith('.tsv'):
        df = pd.read_csv(file_path, sep='\t')
    else:
        df = pd.read_excel(file_path, engine='openpyxl') 

    seen_names = set()
    modified_names = {}

    primers_all = []

    record = SeqIO.read(args.input_file, "fasta")


    # check if regions overlap
    # Sort the DataFrame by 'start' column
    df = df.sort_values(by='start')
    overlap = df.loc[(df['start'] < df['end'].shift()) & (df['end'] > df['start'].shift())]
    if not overlap.empty:
        print("Overlaps found, please edit input")
        print(overlap)
        exit(1)


    # Iterate over each row
    # from tqdm import tqdm
    # for index, row in tqdm(df.iterrows(), total=len(df)):
    for index, row in df.iterrows():
        # Access the "start" and "end" values using the column names
        # name=str(row['name'])
        start = int(row['start'])
        end = int(row['end'])
        name = 'Target ' + str(start) + "-" + str(end)

        if end - start > int(product_size_range[1]):
            print("Region " + str(start) + ", ", str(end) + " larger than product_size_range!")
            pass

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

        # copy background fasta
        tmp_dir = tempfile.TemporaryDirectory()
        tmp_back_n = tmp_dir.name
        if args.background != '':
            shutil.copy(args.background, tmp_back_n)
            tmp_back = os.path.join(tmp_back_n, os.path.basename(args.background))
        else:
            tmp_back = os.path.join(tmp_back_n, 'no_back.fasta')
            # with open(tmp_back, "a") as fasta_file:
            #     fasta_file.write(">random_testing\n")
            #     fasta_file.write("gcaggcaggcaggcag" + "\n")
            #     fasta_file.write(">random_testing2\n")
            #     fasta_file.write("gcaggcaggcaggccag" + "\n")
            #     # Create a temporary file for writing
            #     with tempfile.NamedTemporaryFile(delete=False) as tmp_back:
            #         # You can write to the temporary file as needed
            #         tmp_back.write(b'')

            # # Use the 'tmp_back.name' as the path to the temporary file

            # # copy background fasta
            # tmp_dir = tempfile.TemporaryDirectory()
            # tmp_back_n = tmp_dir.name
            # if args.background != '':
            #     shutil.copy(args.background, tmp_back_n)
            #     tmp_back = os.path.join(tmp_back_n, os.path.basename(args.background))
            # else:
            #     tmp_back = os.path.join(tmp_back_n, "tmp.fasta")
            #     open(tmp_back_n, 'w').close()
            #     tmp_back.write(b'')

            # print(tmp_back)
            # if args.background != '':
            #     with tempfile.TemporaryDirectory() as tmp_back:
            #         # copy the file to the temporary directory
            #         shutil.copy(args.background, tmp_back)
            #         # construct the destination path in the temporary directory
            #         tmp_back = os.path.join(tmp_back, os.path.basename(args.background))
            # else:
            #     # create a temporary file with a .fasta extension
            #     with tempfile.NamedTemporaryFile(suffix='.fasta', delete=False) as tmp_back:
            #         # write an empty string to the file
            #         tmp_back.write(b'')

            # add two random lines for initialization
            with open(tmp_back, "a") as fasta_file:
                fasta_file.write(">appended_rdm_sequence1\n")
                fasta_file.write("gcagcagtcgctcgatccgat" + "\n")
                fasta_file.write(">appended_rdm_sequence2\n")
                fasta_file.write("gcagcagtcggagatagacctcgatccgat")

            # # Read the content of tmp_back and filter out blank lines, in case there were in the input file
            # with open(tmp_back, "r") as fasta_file:
            #     lines = fasta_file.readlines()
            #     filtered_lines = [line.strip() for line in lines if line.strip()]

            # # Write the filtered content back to tmp_back
            # with open(tmp_back, "w") as fasta_file:
            #     fasta_file.write("\n".join(filtered_lines))

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
                background=tmp_back,
                no_self_background=args.no_self_background
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
    cluster_method = True

    names = set([primer['name'] for primer in primers_all])  # get unique names
    if len(names) <= 1:
        print("not enough targets! " + ' '.join(map(str, names)))
        exit(1)
    else:
        print("Finding best primer set for the following targets: " + ', '.join(map(str, names)))
        if cluster_method:
            best_comb = cluster_primers(primers_all)
        else:
            best_score = {'mean': float('inf'), 'range': float('inf'), 'tmp': float('inf'), 'size': float('inf')}
            best_comb = None
            
            with concurrent.futures.ThreadPoolExecutor() as executor:
                # Pass primers_all to evaluate_combination
                futures = [executor.submit(evaluate_combination, comb, primers_all) 
                          for comb in combinations]
                
                for future in concurrent.futures.as_completed(futures):
                    score, comb = future.result()
                    if best_score['mean'] > score['mean']:
                        if best_score['range'] * 0.9 > score['range']:
                            best_score = score
                            best_comb = comb
                        elif best_score['tmp'] > 4 or score['tmp'] > 4 and best_score['tmp'] > score['tmp']:
                            best_score = score
                            best_comb = comb
            
        ######
        print("")
        print("writing file: " + args.output + '.xlsx')
        # convert best_comb tuple to DataFrame object
        best_comb = pd.DataFrame(list(best_comb))
        # add in NEB calc.
        best_comb['Forward Primer TM NEB'] = best_comb['Forward Primer'].apply(q5_melting_temp)
        best_comb['Reverse Primer TM NEB'] = best_comb['Reverse Primer'].apply(q5_melting_temp)

        best_comb.to_excel(args.output + '.xlsx', sheet_name='PrimerSet', index=False, engine='openpyxl') 


    tmp_dir.cleanup()

# (plt, tab) = plot_cluster(clusterer, primers)
# print(tab)

if __name__ == '__main__':
    main()