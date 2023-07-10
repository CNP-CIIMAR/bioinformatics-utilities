## Autor: Leandro de Mattos Pereira - Junior Researcher
## CNP team - Leao Lab - Pedro Leao Team Leader Researcher

import os
import subprocess
import sys
import shutil

def process_contigs(directory):
    output_dir = "16SrRNAmultigenomes"
    os.makedirs(output_dir, exist_ok=True)

    for filename in os.listdir(directory):
        if filename.endswith(".fna"):
            input_file = os.path.join(directory, filename)
            output_file = os.path.join(output_dir, f"{os.path.splitext(filename)[0]}.16Srrna")

            # Executa o comando barrnap
            command = f"./barrnap {input_file} > {output_file}"
            try:
                subprocess.run(command, shell=True, check=True)
                print(f"Arquivo processado: {filename}")
            except subprocess.CalledProcessError as e:
                print(f"Erro ao processar o arquivo {filename}: {e}")
                continue

            # Move o arquivo de saída para o diretório de output
            output_path = os.path.join(output_dir, f"{os.path.splitext(filename)[0]}_processado.rrna.fa")
            shutil.move(output_file, output_path)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Uso: python script.py <diretorio>")
        sys.exit(1)

    directory = sys.argv[1]
    process_contigs(directory)

