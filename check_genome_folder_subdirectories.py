import os
import argparse

def read_genome_list(file_path):
    with open(file_path, 'r') as f:
        return set(line.strip() for line in f.readlines())

def list_subdirectories(dir_path):
    return set(name for name in os.listdir(dir_path) if os.path.isdir(os.path.join(dir_path, name)))

def main():
    parser = argparse.ArgumentParser(description='Verifique quais subdiretórios estão e não estão presentes com base em um arquivo de lista.')
    parser.add_argument('dir_path', type=str, help='Caminho para o diretório contendo os subdiretórios.')
    parser.add_argument('file_path', type=str, help='Caminho para o arquivo contendo a lista de genomas.')
    
    args = parser.parse_args()
    dir_path = args.dir_path
    file_path = args.file_path
    
    genome_list = read_genome_list(file_path)
    subdirs = list_subdirectories(dir_path)
    
    present = genome_list.intersection(subdirs)
    not_present = genome_list.difference(subdirs)
    
    with open('output_info.txt', 'w') as f:
        f.write("Subdiretórios presentes:\n")
        for name in present:
            f.write(f"{name}\n")
        
        f.write("\nSubdiretórios não presentes:\n")
        for name in not_present:
            f.write(f"{name}\n")

if __name__ == '__main__':
    main()
