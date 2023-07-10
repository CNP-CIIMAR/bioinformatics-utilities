#!/usr/bin/python3
# Authors: Leandro de Mattos Pereira
# Pedro Leao - Team Leader Researcher 
# Get genome IDs from a list of proteins

import subprocess
import sys

# Check if the script has been provided with the necessary arguments
if len(sys.argv) < 3:
    print("Usage: python get_genome_ids.py <input_filename> <output_filename>")
    print("Description:")
    print("This script takes a file containing protein accessions as input and retrieves the corresponding genome IDs.")
    print("The results are saved in a tabular format in the specified output file.")
    sys.exit(1)

# Input and output filenames provided as command-line arguments
input_filename = sys.argv[1]
output_filename = sys.argv[2]

# Open the input file and process each line
with open(input_filename, 'r') as input_file, open(output_filename, 'w') as output_file:
    for line in input_file:
        line = line.strip()
        protein_accessions = line.split()
        for protein_accession in protein_accessions:
            if protein_accession:
                ipg_result = subprocess.check_output(['efetch', '-db', 'protein', '-id', protein_accession, '-format', 'ipg']).decode().strip()
                output_file.write("{}\t{}\n".format(protein_accession, ipg_result))

print(f"Results saved to {output_filename}")
