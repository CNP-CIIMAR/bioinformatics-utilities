import argparse
from Bio import SeqIO

def filter_sequences(input_file, output_file):
    # Store unique sequence IDs
    unique_ids = set()

    # Filter the sequences and write unique IDs to the output file
    with open(input_file) as f, open(output_file, 'w') as out_f:
        for record in SeqIO.parse(f, 'fasta'):
            sequence_id = record.id
            if sequence_id.startswith('>'):
                sequence_id = sequence_id[1:]
            if sequence_id not in unique_ids:
                SeqIO.write(record, out_f, 'fasta')
                unique_ids.add(sequence_id)

    return len(unique_ids)

def main():
    parser = argparse.ArgumentParser(description='Remove duplicate sequence IDs from FASTA file.')
    parser.add_argument('input', type=str, help='Input FASTA file')
    parser.add_argument('output', type=str, help='Output FASTA file')
    args = parser.parse_args()

    unique_count = filter_sequences(args.input, args.output)

    print(f"Total unique sequence IDs: {unique_count}")
    print("Filtration completed successfully.")

if __name__ == "__main__":
    main()
