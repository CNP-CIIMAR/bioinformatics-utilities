import os
import shutil
import sys

def move_matching_files(source_dir, exclude_file, destination_dir):
    # Cria o diretório de destino se ele não existir
    if not os.path.exists(destination_dir):
        os.makedirs(destination_dir)

    # Lê o arquivo exclude_genome.txt e armazena os nomes
    with open(exclude_file, 'r') as file:
        names = [line.strip() for line in file]

    # Procura por arquivos correspondentes e os move
    for file in os.listdir(source_dir):
        for name in names:
            if file.startswith(name):
                shutil.move(os.path.join(source_dir, file), os.path.join(destination_dir, file))
                break

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <source_dir> <exclude_file> <destination_dir>")
        sys.exit(1)

    source_dir = sys.argv[1]
    exclude_file = sys.argv[2]
    destination_dir = sys.argv[3]

    move_matching_files(source_dir, exclude_file, destination_dir)
