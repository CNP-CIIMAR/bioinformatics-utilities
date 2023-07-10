##Autores: Leandro de Mattos Pereira
##### Pedro Leão, Team Leader - CNP Team - CIIMAR
### Script function: Copy files indicated in one file.txt from several subdirectories to other specified directory.

import os
import shutil
import sys

# Verifica se foram passados argumentos suficientes
if len(sys.argv) < 3:
    print("Uso: python script.py <arquivo_texto> <diretorio_saida>")
    sys.exit(1)

# Obtém os argumentos de linha de comando
arquivo_texto = sys.argv[1]
diretorio_saida = sys.argv[2]

# Verifica se o arquivo de texto existe
if not os.path.isfile(arquivo_texto):
    print(f"O arquivo {arquivo_texto} não existe.")
    sys.exit(1)

# Cria o diretório de saída se ele não existir
os.makedirs(diretorio_saida, exist_ok=True)

# Lê os nomes dos arquivos do arquivo de texto
with open(arquivo_texto, 'r') as file:
    nomes_arquivos = file.read().splitlines()

# Percorre todos os subdiretórios e arquivos no diretório atual
for root, dirs, files in os.walk(os.getcwd()):
    for file in files:
        # Verifica se o arquivo possui a extensão .fna e corresponde a um dos nomes dos arquivos
        if file.endswith('.fna') and any(nome_arquivo in file for nome_arquivo in nomes_arquivos):
            # Caminho completo do arquivo de entrada
            caminho_arquivo_entrada = os.path.join(root, file)
            # Caminho completo do arquivo de saída
            caminho_arquivo_saida = os.path.join(diretorio_saida, file)
            # Verifica se o arquivo de entrada e o arquivo de saída são os mesmos
            if os.path.abspath(caminho_arquivo_entrada) != os.path.abspath(caminho_arquivo_saida):
                # Copia o arquivo para o diretório de saída
                shutil.copy(caminho_arquivo_entrada, caminho_arquivo_saida)
