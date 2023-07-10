#!/usr/bin/python3
#autor: Leandro de Mattos Pereira, PhD in Computational Biology and Systems
## 05/16/2023

import sys

def comparar_arquivos(arquivo1, arquivo2):
    with open(arquivo1, 'r') as file1, open(arquivo2, 'r') as file2:
        linhas_arquivo1 = file1.readlines()
        linhas_arquivo2 = file2.readlines()

    diferentes = []

    for linha in linhas_arquivo1:
        if linha not in linhas_arquivo2:
            diferentes.append((linha.strip(), arquivo1))

    for linha in linhas_arquivo2:
        if linha not in linhas_arquivo1:
            diferentes.append((linha.strip(), arquivo2))

    return diferentes

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Por favor, forneça o nome de dois arquivos como argumentos.")
    else:
        arquivo1 = sys.argv[1]
        arquivo2 = sys.argv[2]

        diferentes = comparar_arquivos(arquivo1, arquivo2)

        if len(diferentes) == 0:
            print("Os arquivos têm o mesmo conteúdo.")
        else:
            print("Linhas diferentes encontradas:")
            for linha, arquivo in diferentes:
                print(f"Linha '{linha}' encontrada no arquivo '{arquivo}'.")
