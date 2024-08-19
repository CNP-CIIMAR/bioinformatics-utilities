## Autor: Leandro de Mattos Pereira
## CNP Laboratory - Leao Pedro.
import pandas as pd
import sys

def filter_sort_and_remove_rows(input_file, e_value_threshold, output_file):
    # Load the data into a pandas DataFrame
    df = pd.read_csv(input_file, sep='\t', engine='python')
    
    # Convert the e-value column to numeric, forcing non-numeric to NaN (Not a Number)
    df['full_sequence_e-value'] = pd.to_numeric(df['full_sequence_e-value'], errors='coerce')
    
    # Remove rows where 'Genome_accession' starts with 'all_fasta_good.quality.faa'
    df = df[~df['Genome_accession'].astype(str).str.startswith('all_fasta_good.quality.faa')]
    
    # Sort the DataFrame based on the e-value column in ascending order
    df_sorted = df.sort_values(by='full_sequence_e-value', ascending=True)
    
    # Filter the DataFrame based on the e-value threshold
    filtered_df = df_sorted[df_sorted['full_sequence_e-value'] <= e_value_threshold]
    
    # Save the filtered DataFrame to the specified output TSV file with tab separator
    filtered_df.to_csv(output_file, sep='\t', index=False)
    print(f'Filtered, sorted, and cleaned data saved to {output_file}')

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python filter_script.py <input_file.tsv> <e_value_threshold> <output_file.tsv>")
        sys.exit(1)
    
    input_filepath = sys.argv[1]
    e_value = float(sys.argv[2])
    output_filepath = sys.argv[3]
    
    filter_sort_and_remove_rows(input_filepath, e_value, output_filepath)
