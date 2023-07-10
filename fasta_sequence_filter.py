#Autor: Leandro de Mattos Pereira. PhD Computational Biology and Systems, Junior Researcher CNP Team - Leao Laboratory.
#Date 05/19/2023
import argparse
from Bio import SeqIO

def filter_sequences(input_file, output_file, ids_file, keep=True):
    # Lê os IDs a serem mantidos ou excluídos a partir do arquivo
    with open(ids_file) as f:
        ids = set(line.strip() for line in f)

    # Filtra as sequências com base nos IDs
    initial_count = 0
    filtered_count = 0
    with open(input_file) as f, open(output_file, 'w') as out_f:
        for record in SeqIO.parse(f, 'fasta'):
            initial_count += 1
            sequence_id = record.id
            if sequence_id.startswith('>'):
                sequence_id = sequence_id[1:]
            if (sequence_id in ids) == keep:
                SeqIO.write(record, out_f, 'fasta')
                filtered_count += 1

    return initial_count, filtered_count

def main():
    parser = argparse.ArgumentParser(description='Filter FASTA file based on sequence IDs.')
    parser.add_argument('input', type=str, help='Input FASTA file')
    parser.add_argument('output', type=str, help='Output FASTA file')
    parser.add_argument('ids', type=str, help='File containing the IDs to filter')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--exclude', action='store_true', help='Exclude the sequences based on the IDs')
    group.add_argument('--keep', action='store_true', help='Keep only the sequences based on the IDs')
    args = parser.parse_args()

    initial_count, filtered_count = filter_sequences(args.input, args.output, args.ids, keep=args.keep)

    if args.exclude:
        print(f"Sequences in input: {initial_count}")
        print(f"Sequences remaining after filtration: {filtered_count}")
    else:
        print(f"Sequences in input: {initial_count}")
        print(f"Sequences generated in output: {filtered_count}")

    print("Filtration completed successfully.")

if __name__ == "__main__":
    main()
    
#Feel free to customize the README file according to your needs. Make sure to update the installation instructions and provide proper credits, license information, and any other relevant details.
#Licence: This project is licensed under the MIT License.
