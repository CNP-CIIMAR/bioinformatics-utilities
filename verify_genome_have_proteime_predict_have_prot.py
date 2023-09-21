import os
import sys

if len(sys.argv) != 2:
    print("Uso: python script.py <diretório>")
    sys.exit(1)

directory = sys.argv[1]

# Lista de subdiretórios com arquivos .faa vazios
subdirectories_with_empty_faa = []

# Percorre todos os subdiretórios dentro do diretório informado
for subdir in os.listdir(directory):
    subdir_path = os.path.join(directory, subdir)
    if os.path.isdir(subdir_path):
        # Verifica se o subdiretório contém arquivos .faa
        faa_files = [file for file in os.listdir(subdir_path) if file.endswith(".faa")]
        if faa_files:
            # Verifica se os arquivos .faa estão vazios
            empty_faa_files = [file for file in faa_files if os.path.getsize(os.path.join(subdir_path, file)) == 0]
            if empty_faa_files:
                subdirectories_with_empty_faa.append(subdir)

# Imprime os subdiretórios que contêm arquivos .faa vazios
if subdirectories_with_empty_faa:
    print("Subdiretórios com arquivos .faa vazios:")
    for subdir in subdirectories_with_empty_faa:
        print(subdir)
else:
    print("Nenhum subdiretório encontrado com arquivos .faa vazios.")
