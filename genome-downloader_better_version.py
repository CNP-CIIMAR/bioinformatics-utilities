import os
import sys
import subprocess

def download_genomes(genome_ids, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    ## Não esqueça de alterar a linha abaixo para o caminho do script datasets do ncbi - obtido no ncbi
    datasets_cmd = "/home/mattoslmp/anaconda3/envs/biopython/bin/datasets"

    for assembly_id in genome_ids:
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
    if len(sys.argv) != 3:
        print("Usage: script.py <genome_ids_file> <output_directory>")
        sys.exit(1)

    genome_ids_file = sys.argv[1]
    output_dir = sys.argv[2]

    # Leia os IDs dos genomas do arquivo
    with open(genome_ids_file, 'r') as f:
        genome_ids = [line.strip() for line in f]

    # Faça o download dos genomas
    download_genomes(genome_ids, output_dir)
