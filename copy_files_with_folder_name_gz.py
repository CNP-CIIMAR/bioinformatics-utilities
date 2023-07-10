#!/usr/bin/python3
## Autor: Leandro de Matttos Pereira
# PhD Computational Biology and Systems - Junior Investigador - CNP Team - Leao Lab
# -*- coding: utf-8 -*-
import os
import sys
import shutil
import gzip

def copy_files_with_folder_name(directory, output_directory):
    # Verifica se o diretório de saída já existe
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Percorre todos os subdiretórios do diretório fornecido
    for root, dirs, files in os.walk(directory):
        for file in files:
            # Obtém o caminho completo do arquivo original
            source_file = os.path.join(root, file)
            # Obtém o nome da pasta pai (subdiretório)
            folder_name = os.path.basename(root)
            # Adiciona o nome da pasta pai ao início do nome do arquivo
            new_file_name = folder_name + '_' + file
            # Obtém o caminho completo do arquivo de destino
            destination_file = os.path.join(output_directory, new_file_name)

            # Verifica se o arquivo está compactado em formato .gz
            if file.endswith('.gz'):
                # Verifica se o arquivo .gz é válido
                try:
                    with gzip.open(source_file, 'rb') as f_in:
                        f_in.read(1)  # Tenta ler pelo menos um byte
                except (IOError, EOFError, OSError, ValueError) as e:
                    print('Erro ao descompactar o arquivo {}: {}'.format(source_file, str(e)))
                    continue

                # Descompacta o arquivo
                with gzip.open(source_file, 'rb') as f_in:
                    with open(destination_file, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)

                # Remove a extensão .gz do nome do arquivo de destino
                destination_file_no_ext = os.path.splitext(destination_file)[0]
                os.rename(destination_file, destination_file_no_ext)
            else:
                # Copia o arquivo para o diretório de destino
                shutil.copy(source_file, destination_file)

if __name__ == '__main__':
    # Verifica o número de argumentos fornecidos
    if len(sys.argv) == 1:
        print('Este script copia arquivos de subdiretórios para um novo diretório, adicionando o nome do subdiretório ao início do nome do arquivo.')
        print('Uso: python script.py <diretório_origem> <diretório_destino>')
        sys.exit(0)
    elif len(sys.argv) != 3:
        print('Argumentos inválidos. Use python script.py --help para obter mais informações.')
        sys.exit(1)

    # Obtém os argumentos de linha de comando
    if sys.argv[1] == '--help':
        print('Este script copia arquivos de subdiretórios para um novo diretório, adicionando o nome do subdiretório ao início do nome do arquivo.')
        print('Uso: python script.py <diretório_origem> <diretório_destino>')
        sys.exit(0)

    source_directory = sys.argv[1]
    output_directory = sys.argv[2]

    # Executa a função de cópia de arquivos
    copy_files_with_folder_name(source_directory, output_directory)
