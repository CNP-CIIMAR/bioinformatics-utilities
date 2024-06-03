import sys
import argparse
import logging

def filter_table(file_path):
    """
    Filter the table based on organism accession numbers.
    Args:
    file_path (str): Path to the input file.

    Returns:
    str: Filtered table as a string.
    """
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
    except Exception as e:
        logging.error(f"Error reading file {file_path}: {e}")
        return ""

    # Remove header and sort the remaining lines
    header = lines[0].strip()
    data_lines = sorted(line.strip() for line in lines[1:] if line.strip())
    
    # Dictionary to store organism accessions
    organism_accessions = {}
    
    for line in data_lines:
        columns = line.split('\t')
        if len(columns) > 2:
            accession = columns[0]
            organism_name = columns[2]
            
            base_accession = accession.split('.')[0]
            
            if organism_name not in organism_accessions:
                organism_accessions[organism_name] = {}
                
            if base_accession not in organism_accessions[organism_name]:
                organism_accessions[organism_name][base_accession] = []
                
            organism_accessions[organism_name][base_accession].append(accession)

    filtered_lines = []

    for line in data_lines:
        columns = line.split('\t')
        if len(columns) > 2:
            accession = columns[0]
            organism_name = columns[2]
            base_accession = accession.split('.')[0]
            
            accessions_list = organism_accessions.get(organism_name, {}).get(base_accession, [])
            
            if accession.startswith('GCF_'):
                if 'GCF_' + base_accession in accessions_list:
                    filtered_lines.append(line)
            elif accession.startswith('GCA_'):
                if 'GCF_' + base_accession not in accessions_list:
                    filtered_lines.append(line)

    filtered_table = header + '\n' + '\n'.join(filtered_lines)
    return filtered_table

def print_help():
    """Print help message."""
    print("Usage: python script.py <path_to_table> <path_to_output>")
    print("Filter the table based on organism accession numbers and generate a new output table.")
    print("Ensure to provide the path to the file containing the table as the first argument and the path to the output file as the second argument.")

def main():
    """Main function to parse arguments and filter table."""
    parser = argparse.ArgumentParser(description="Filter a table based on organism accession numbers.")
    parser.add_argument("table_path", help="Path to the input table file")
    parser.add_argument("output_path", help="Path to the output file")

    args = parser.parse_args()

    table_path = args.table_path
    output_path = args.output_path

    logging.basicConfig(level=logging.INFO)
    logging.info("Starting the filtering process.")

    filtered_table = filter_table(table_path)

    try:
        with open(output_path, 'w') as output_file:
            output_file.write(filtered_table)
        logging.info(f"Filtered table has been written to {output_path}")
    except Exception as e:
        logging.error(f"Error writing to file {output_path}: {e}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print_help()
        sys.exit(1)
    main()
