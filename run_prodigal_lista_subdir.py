import os
import sys
import subprocess

if len(sys.argv) != 3:
    print("Uso: python script.py <diretorio_data> <arquivo_lista_subdiretorios>")
    sys.exit(1)

data_directory = sys.argv[1]
subdirectories_file = sys.argv[2]

# Lê a lista de subdiretórios do arquivo
with open(subdirectories_file, "r") as f:
    subdirectories = [line.strip() for line in f]

# Percorre os subdiretórios da lista
for subdir in subdirectories:
    subdir_path = os.path.join(data_directory, subdir)

    # Verifica se o subdiretório existe
    if os.path.isdir(subdir_path):
        input_genome_files = [file for file in os.listdir(subdir_path) if file.endswith(".fna")]

        # Verifica se há pelo menos um arquivo .fna no subdiretório
        if input_genome_files:
            # Executa o Prodigal para cada arquivo .fna
            for input_genome_file in input_genome_files:
                input_genome_path = os.path.join(subdir_path, input_genome_file)
                output_gbk_path = os.path.join(subdir_path, "genomic.gbk")

                # Executa o Prodigal para gerar o arquivo GBK
                subprocess.run(["prodigal", "-i", input_genome_path, "-o", output_gbk_path, "-f", "gbk"])

# Imprime mensagem de conclusão
print("Prodigal concluído para os subdiretórios listados.")
