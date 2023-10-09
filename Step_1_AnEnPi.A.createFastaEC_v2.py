### Autor: Leandro de Mattos Pereira
### Junior Researcher - Pedro Leao Lab - CNP Team
import re
import sys
import os
import logging
from collections import Counter

def parse_ec_numbers(header):
    ec_pattern = re.compile(r"(EC:|EC\s)\s*(\d+\.\d+\.\d+\.\d+|\d+\.\d+\.\d+|\d+\.\d+|\d+)")
    ec_numbers = ec_pattern.findall(header)
    return [ec[1] for ec in ec_numbers]

def split_fasta(input_file, output_directory):
    with open(input_file, 'r') as file:
        fasta_data = file.read()

    fasta_records = fasta_data.split('>')[1:]

    ec_counts = Counter()  # Initialize a Counter to count EC numbers by class

    for record in fasta_records:
        if record:
            lines = record.split('\n')
            header = lines[0]
            sequence = ''.join(lines[1:])
            ec_numbers = parse_ec_numbers(header)

            for ec_number in ec_numbers:
                output_file = output_directory + '/EC_' + ec_number + '.fasta'
                with open(output_file, 'a') as file:
                    file.write('>' + header + '\n' + sequence + '\n')

                ec_class = ec_number.split('.')[0]  # Extract the EC class (first number)
                ec_counts[ec_class] += 1  # Increment the count for the EC class

    print("Fasta files have been created successfully.")

    # Write EC counts by class to a text file
    output_file_path = output_directory + '/ec_counts.txt'
    with open(output_file_path, 'w') as file:
        file.write("EC counts by class:\n")
        for ec_class in ['1', '2', '3', '4', '5', '6', '7']:
            count = ec_counts.get(ec_class, 0)
            file.write(f"EC {ec_class}: {count}\n")

    # Log EC counts by class
    log_filename = output_directory + '/ec_counts.log'
    logging.basicConfig(filename=log_filename, level=logging.INFO, format='%(message)s')
    for ec_class in ['1', '2', '3', '4', '5', '6', '7']:
        count = ec_counts.get(ec_class, 0)
        logging.info(f"EC {ec_class}: {count}")

    print(f"EC counts have been written to {output_file_path}.")
    print(f"EC counts have been logged to {log_filename}.")

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python split_fasta.py input.pep output_directory")
        sys.exit(1)

    input_file = sys.argv[1]
    output_directory = sys.argv[2]

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    split_fasta(input_file, output_directory)
