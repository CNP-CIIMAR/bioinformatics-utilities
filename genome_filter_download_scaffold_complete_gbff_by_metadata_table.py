import pandas as pd
import os
import sys
import subprocess

def filter_table(input_file, output_file):
    # Leia a tabela
    df = pd.read_csv(input_file, delimiter="\t")

    # Filtrar as linhas de acordo com as condições especificadas
    filtered_df = df[(df['Assembly Level'].isin(['Scaffold', 'Complete'])) & 
                     (df['Lineage'].str.contains('Cyanobacteriota|Pseudomonadota|Myxococcota|Actinomycetota|Eukaryota'))]

    # Salvar a tabela filtrada
    filtered_df.to_csv(output_file, sep='\t', index=False)

    return filtered_df

def download_genomes(df, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    datasets_cmd = "/home/mattoslmp/anaconda3/envs/biopython/bin/datasets"

    for assembly_id in df['Assembly']:
        output_path = os.path.join(output_dir, f"{assembly_id}.zip")
        cmd = [
            datasets_cmd, 'download', 'genome', 'accession', assembly_id,
            '--include', 'gbff',
            '--filename', output_path
        ]
        try:
            subprocess.run(cmd, check=True)
            print(f"Downloaded {output_path}")
        except subprocess.CalledProcessError as e:
            print(f"Failed to download {assembly_id}, error: {str(e)}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: script.py <input_file> <output_file> <output_directory>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    output_dir = sys.argv[3]

    filtered_df = filter_table(input_file, output_file)
    download_genomes(filtered_df, output_dir)

