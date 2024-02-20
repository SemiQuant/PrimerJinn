#!/bin/python3
import os
import shutil # we're using f"" so we must be on a recent enough Python for shutil.which
import subprocess
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp
from Bio import SeqIO
import argparse
import itertools
import json
import requests
import primer3

# Define the arguments to be parsed
parser = argparse.ArgumentParser(description='PCR parameter inputs')
parser.add_argument('--target_tm', type=int, help='Annealing temperature (in Celsius)', default = 55.0)
parser.add_argument('--salt_concentration', type=float, help='Salt concentration (in nM)', default = 50)
parser.add_argument('--product_size_max', type=float, help='maximum length of PCR products in nucleotides', default = 2000)
parser.add_argument('--req_five', type=str, help="Require the 5' end of the primer to bind?", default = True)
parser.add_argument('--primer_seq', type=str, help="Primer sequences, one per line", required = True)
parser.add_argument('--input_file', type=str, help="Reference FASTA file", required = True)
parser.add_argument('--output', type=str, help="Output file name", default = "in_silico_PCR")
parser.add_argument('--Q5', action='store_true', help='Whether to use Q5 approximation settings for Tm calculations.', default=False)

# check if blastn is in the system path
if not shutil.which("blastn"):
    print('Error: BLAST is not installed.')
    exit(1)

# Check if required arguments are missing and print help message
args = parser.parse_args()
if not args.primer_seq or not args.input_file:
    parser.print_help()
    exit()

# Define the inputs
annealing_temp = args.target_tm - 4
salt_conc = args.salt_concentration
max_amplicon_len = args.product_size_max
req_five = args.req_five
primer_seq = args.primer_seq
ref_fasta_file = args.input_file
out_file = args.output


def q5_melting_temp(seq1, seq2, salt=0.5):
    url = "https://tmapi.neb.com/tm/q5/%s/%s/%s?&fmt=short" % (salt, seq1, seq2)
    response = requests.get(url)
    json_string = response.json()
    tm1 = json_string["data"]["tm1"]
    tm2 = json_string["data"]["tm2"]
    return [tm1, tm2]

def calc_tm(seq, seq2='', Q5=args.Q5):
    if Q5:
        # tm = MeltingTemp.Tm_NN(seq, c_seq = seq2, nn_table=MeltingTemp.DNA_NN4, dnac1=2500, dnac2=2500, Na=40, K=0, Tris=0, Mg=0.4, dNTPs=0, saltcorr=5.0)
        # tm = q5_melting_temp(str(seq), str(seq2), salt_conc/1000)
        if seq2 == '':
            tm = primer3.bindings.calcTm(seq, dv_conc=2, mv_conc=70, dna_conc=3300)
        else:
            primer3.calcHeterodimer(seq, seq2, temp_c=100, dv_conc = 2, mv_conc = 70, dna_conc = 3300)
    else:
        tm = MeltingTemp.Tm_NN(seq, c_seq = seq2, Na=salt_conc, dnac1=250, dnac2=0, saltcorr=7, nn_table=MeltingTemp.DNA_NN4, selfcomp=False, check=True, shift=0.0)
    return tm

def find_binding_positions(primer_seq, ref_fasta_file, annealing_temp, req_five, salt_conc):
    # Run BLAST to search for primer sequence against reference FASTA
    try:
        r = subprocess.run(["blastn","-query"]+primer_seq.split()+["-subject",ref_fasta_file,"-task","blastn-short","-outfmt",'6 qseqid sseqid qstart qend sstart send qseq sseq mismatch length'],capture_output=True,check=True)
    except subprocess.CalledProcessError as e:
        print(f"BLAST execution failed with exit code {e.returncode}.")
        print(e.stderr.decode('utf-8'))
        raise
    blast_df = pd.DataFrame(columns=['Query ID', 'Subject ID', 'Query Start', 'Query End', 'Subject Start', 'Subject End', 'Query Sequence Match', 'Direction', 'Binding Position', 'Mismatches', 'Binding Length'])
    blast_output = r.stdout.decode('utf-8').splitlines()
    for line in blast_output:
            fields = line.strip().split('\t')
            qseqid, sseqid, qstart, qend, sstart, send, qseq, sseq, mismatch, length = fields
            direction = '+' if int(sstart) < int(send) else '-'
            binding_pos = sstart if direction == '+' else send
            # Check if primer binds at the specified annealing temperature
            if req_five and int(qstart) >= 2:
                continue
            # tm = MeltingTemp.Tm_NN(Seq(sseq), Na=salt_conc, dnac1=250, dnac2=0, saltcorr=7, nn_table=MeltingTemp.DNA_NN4, selfcomp=False, check=True, shift=0.0)
            tm = calc_tm(Seq(sseq))
            if annealing_temp <= tm:
                blast_df.loc[len(blast_df)] = [qseqid, sseqid, qstart, qend, sstart, send, sseq, direction, binding_pos, mismatch, length]
    return blast_df

def find_compatible_pairs(blast_df, max_len):
    pairs = itertools.combinations(blast_df.index, 2)
    compatible_pairs = []
    for pair in pairs:
        row1 = blast_df.loc[pair[0]]
        row2 = blast_df.loc[pair[1]]
        amp_size = int(row1['Binding Position']) - int(row2['Binding Position'])
        if abs(amp_size) <= max_len and row1['Direction'] != row2['Direction'] and row1['Subject ID'] == row2['Subject ID']:
            if row1['Direction'] == '+' and amp_size < 0:
                compatible_pairs.append({'qseq1': row1['Query Seq'], 'qseq1_input': row1['Sequence'], 'qstart1': row1['Query Start'], 'qend1': row1['Query End'], 'direction1': row1['Direction'], 'mismatch1': row1['Mismatches'],
                                         'qseq2': row2['Query Seq'], 'qseq2_input': row2['Sequence'], 'qstart2': row2['Query Start'], 'qend2': row2['Query End'], 'direction2': row2['Direction'], 'mismatch2': row2['Mismatches'],
                                         'binding_pos_diff': abs(amp_size), 'reference': row1['Subject ID'], 'ref_region': str(row2['Binding Position']) + ", " + str(row1['Binding Position'])})
            elif row1['Direction'] == '-' and amp_size > 0:
                compatible_pairs.append({'qseq1': row2['Query Seq'], 'qseq2_input': row2['Sequence'], 'qstart1': row2['Query Start'], 'qend1': row2['Query End'], 'direction1': row2['Direction'], 'mismatch2': row2['Mismatches'],
                                         'qseq2': row1['Query Seq'], 'qseq1_input': row1['Sequence'], 'qstart2': row1['Query Start'], 'qend2': row1['Query End'], 'direction2': row1['Direction'], 'mismatch1': row1['Mismatches'],
                                         'binding_pos_diff': abs(amp_size), 'reference': row1['Subject ID'], 'ref_region': str(row1['Binding Position']) + ", " + str(row2['Binding Position'])})
            else:
                continue
    return pd.DataFrame(compatible_pairs)

def find_oligo_dimers(fasta_file, temp, salt_conc):
    # Read in the fasta sequences
    seq_records = list(SeqIO.parse(fasta_file, "fasta"))

    # Create a list to hold the dimer pairs and melting temperatures
    dimers = []
    melting_temps = []

    # Loop through all possible pairs of sequences
    for seq1, seq2 in itertools.combinations(seq_records, 2):
        # Calculate the melting temperature of the dimer
        # dimer_tm = MeltingTemp.Tm_NN(str(seq1.seq + seq2.seq), nn_table=MeltingTemp.DNA_NN4, Na=salt_conc, saltcorr=7)
        dimer_tm = calc_tm(str(seq1.seq + seq2.seq))

        # Check if the dimer's melting temperature is above the given temperature
        if dimer_tm > temp:
            # Add the dimer's names and melting temperature to the lists
            dimers.append((seq1.id, seq2.id))
            melting_temps.append(round(dimer_tm))

    # Create a pandas dataframe from the dimer pairs and melting temperatures
    df = pd.DataFrame({'Sequence1': [d[0] for d in dimers],
                       'Sequence2': [d[1] for d in dimers],
                       'MeltingTemp': melting_temps})

    return df

def is_valid_dna_sequence(sequence):
    """Checks if a DNA sequence only contains valid characters."""
    valid_chars = set('ACGTN')
    return all(char in valid_chars for char in sequence.upper())

def main():
    i=0
    with open(primer_seq, 'r') as infile, open(primer_seq+'.fasta', 'w') as outfile:
        # Initialize a dictionary to keep track of sequence names
        seen_names = {}
        i+=1
        for line in infile:
            columns = line.strip().split('\t')
            if len(columns) == 1:
                # Remove any leading/trailing whitespace from the sequence
                sequence = line.strip()
                if not (is_valid_dna_sequence(sequence)):
                    print("error in sequence on: " + line)
                    break
                # Generate a unique sequence name
                name = f'sequence{i}'
            else:
                sequence, name = columns
                if not (is_valid_dna_sequence(sequence)):
                    print("error in sequence on: " + name)
                    break
            # Check if the name has already been used
            if name in seen_names:
                # Increment the count for this name
                count = seen_names[name] + 1
                # Generate a new name by appending the count to the original name
                new_name = f'{name}_{count}'
                # Update the dictionary with the new count
                seen_names[name] = count
                # Use the new name for the sequence
                name = new_name
            else:
                # Add the name to the dictionary with a count of 1
                seen_names[name] = 1
            # Write the sequence to the output file in FASTA format
            outfile.write(f'>{name}\n{sequence}\n')


    # primer dimers
    dimers = find_oligo_dimers(primer_seq+'.fasta', annealing_temp, salt_conc)
    # print("Writing output to " + out_file + '_primer_dimears.tsv')
    # dimers.to_csv(out_file + '_primer_dimers.tsv', index=False, sep="\t")


    blast_df = find_binding_positions(primer_seq+'.fasta', ref_fasta_file, annealing_temp, req_five, salt_conc)

    # Read in the FASTA file as a dictionary of SeqRecord objects
    fasta_dict = SeqIO.to_dict(SeqIO.parse(primer_seq+'.fasta', "fasta"))

    # Iterate through the rows of the DataFrame and replace the add primer with sequences
    for i, row in blast_df.iterrows():
        query_id = row["Query ID"]
        if query_id in fasta_dict:
            sequence = str(fasta_dict[query_id].seq)
            blast_df.at[i, "Sequence"] = sequence
        else:
            blast_df.at[i, "Sequence"] = ""

    blast_df = blast_df.rename(columns={'Query ID': 'Query Seq'})

    amp_df = find_compatible_pairs(blast_df, max_len=max_amplicon_len)
    # print("Writing output to " + out_file + '.tsv')
    # amp_df.to_csv(out_file + '.tsv', index=False, sep="\t")



    # get interactions
    amplicons = []
    ref_seq = SeqIO.read(ref_fasta_file, "fasta").seq
    for index, row in amp_df.iterrows():
        ref_region = row['ref_region']
        #add in the primer sequences
        if row['direction1'] == '-':
            pf_seq = Seq(row['qseq1_input'])
            pf_seq = str(pf_seq.reverse_complement())
        else:
            pf_seq = Seq(row['qseq1_input'])
        if row['direction2'] == '-':
            pr_seq = Seq(row['qseq2_input'])
            pr_seq = str(pf_seq.reverse_complement())
        else:
            pr_seq = Seq(row['qseq2_input'])
        ref_start, ref_end = ref_region.split(",")
        ref_end = int(ref_end)
        ref_start = int(ref_start)
        if ref_start > ref_end:
            amplicon_seq = ref_seq[ref_end-1:ref_start]
        else:
            amplicon_seq = ref_seq[ref_start-1:ref_end]
        amplicon_seq = str(pf_seq) + str(amplicon_seq) + str(pr_seq)
        amplicons.append({'Forward_primer':  row['qseq1'], 'Reverse_primer':  row['qseq2'], 'amplicon_seq': amplicon_seq})


    tm_data = []
    for pair in itertools.combinations(amplicons, 2):
        # tm = MeltingTemp.Tm_NN(pair[0]['amplicon_seq'], pair[1]['amplicon_seq'], Na=salt_conc, saltcorr=7)
        try:
            tm = calc_tm(seq=pair[0]['amplicon_seq'], seq2=pair[1]['amplicon_seq'])
            tm_data.append({"amplicon1_PF": pair[0]['Forward_primer'],
                "amplicon1_PR": pair[0]['Reverse_primer'],
                "amplicon2_PF": pair[1]['Forward_primer'],
                "amplicon2_PR": pair[1]['Reverse_primer'], 
                'tm': str(tm)})
        except:
            continue

    tm_df = pd.DataFrame(tm_data)
    # print("Writing output to " + out_file + '_amplicon_interactions.tsv')
    # tm_df.to_csv(out_file + '_amplicon_interactions.tsv', index=False, sep="\t")


    with pd.ExcelWriter(out_file + '.xlsx', engine='xlsxwriter') as writer:
        amp_df.to_excel(writer, sheet_name='amplicons', index=False)
        dimers.to_excel(writer, sheet_name='dimers', index=False)
        tm_df.to_excel(writer, sheet_name='interactions', index=False)

    # cleanup
    os.remove(primer_seq+'.fasta')


if __name__ == '__main__':
    main()


