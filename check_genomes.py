import os
import sys
import shutil
import pandas as pd
from tabulate import tabulate
## Autor: Leandro de Mattos Pereira
## 09/05/2024


def check_genomes_in_directory(genomes_file, directory, output_file, new_directory):
    # Reading genomes from the file
    with open(genomes_file, 'r') as file:
        input_genomes = [line.strip() for line in file.readlines()]
    
    # Reading genomes from the directory
    directory_genomes = os.listdir(directory)
    
    # Checking matches
    matched_genomes = [genome for genome in input_genomes if genome in directory_genomes]
    unmatched_input_genomes = [genome for genome in input_genomes if genome not in directory_genomes]
    unmatched_directory_genomes = [genome for genome in directory_genomes if genome not in input_genomes]
    
    # Creating DataFrame for the table
    data = {
        'Genomas_arquivo_input': pd.Series(input_genomes),
        'Genomas_match_dir': pd.Series(matched_genomes),
        'Genomas_only_dir': pd.Series(unmatched_directory_genomes)
    }
    df = pd.DataFrame(data)
    
    # Ensure new directory exists
    os.makedirs(new_directory, exist_ok=True)
    
    # Move matched files to the new directory
    for genome in matched_genomes:
        shutil.move(os.path.join(directory, genome), os.path.join(new_directory, genome))
    
    # Saving the table to a tabular format file
    with open(output_file, 'w') as f:
        f.write(tabulate(df, headers='keys', tablefmt='grid'))
        # Adding summary statistics at the end of the file
        f.write("\n\nSummary Statistics:\n")
        f.write(f"Total Matches: {len(matched_genomes)}\n")
        f.write(f"Genomes Only in Directory: {len(unmatched_directory_genomes)}\n")
    
    print(tabulate(df, headers='keys', tablefmt='grid'))
    print("\nSummary Statistics:")
    print(f"Total Matches: {len(matched_genomes)}")
    print(f"Genomes Only in Directory: {len(unmatched_directory_genomes)}")
    
    return df

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py <path_to_genomes_file> <directory> <output_file> <new_directory>")
        sys.exit(1)
    
    genomes_file = sys.argv[1]
    directory = sys.argv[2]
    output_file = sys.argv[3]
    new_directory = sys.argv[4]
    
    # Running the function
    df = check_genomes_in_directory(genomes_file, directory, output_file, new_directory)
