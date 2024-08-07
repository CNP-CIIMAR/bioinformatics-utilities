import sys

def read_genome_ids(file_path):
    with open(file_path, 'r') as file:
        return [line.strip() for line in file]

def extract_prefix(genome_id):
    parts = genome_id.split('_')
    if len(parts) > 2:
        return '_'.join(parts[:2])
    return genome_id

def write_prefixes_to_file(genome_list, output_file_name):
    with open(output_file_name, 'w') as file:
        for genome_id in genome_list:
            prefix = extract_prefix(genome_id)
            file.write(prefix + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    genome_list = read_genome_ids(input_file)
    write_prefixes_to_file(genome_list, output_file)

    print(f"Prefixes extracted and saved to {output_file}")
