import csv
import argparse

# 01_Variant_Prep_FKRP.py
# 
# This script processes variant count files from the Gargamel pipeline
# Usage: 1_Variant_Prep_FKRP.py -o fkrp_variantCounts_site.txt normalized_FKRP-1_func normalized_FKRP-2_func normalized_FKRP-0_func

def read_and_process_file(filename, full_length, seq_dict):
    with open(filename, 'r') as file:
        lines = file.readlines()
        for line in lines[1:]:
            line = line.rstrip().split("\t")
            site = int(line[2])
            mut_seq = full_length[:site-1] + line[5] + full_length[site:]
            high = int(line[11])
            low = int(line[12])
            if mut_seq in seq_dict:
                seq_dict[mut_seq].extend([high, low])
            else:
                seq_dict[mut_seq] = [site, high, low]

parser = argparse.ArgumentParser(description='Process variant count files from Gargamel')
parser.add_argument('input_files', metavar='N', type=str, nargs='+', help='input variant count files')
parser.add_argument('-o', '--output', type=str, default='output.txt', help='output file name')
args = parser.parse_args()

# Initialize variables
fkrp = []
full_length = "A" # Start codon
old_site = "1"
seq_dict = {}

# Construct the full WT CDS sequence
with open(args.input_files[0], 'r') as file:
    lines = file.readlines()
    fkrp = [line.rstrip().split("\t") for line in lines]

for i in fkrp[1:]:
    site = i[2]
    if site != old_site:
        old_site = site
        full_length += i[4]

# Process variant count files
for input_file in args.input_files:
    read_and_process_file(input_file, full_length, seq_dict)

# Prepare the output
header = ['nt_seq', 'site', 'fkrp1high', 'fkrp1low', 'fkrp2high', 'fkrp2low', 'fkrp0high', 'fkrp0low']
out = [header]

for nt_seq, values in seq_dict.items():
    out.append([nt_seq] + values)

# Write the output
with open(args.output, 'w', newline='') as file:
    writer = csv.writer(file, delimiter='\t')
    writer.writerows(out)
