import os
import argparse

def rename_files(diretorio_raiz):
    # Percorre todos os subdiretórios do diretório raiz
    for subdiretorio in os.listdir(diretorio_raiz):
        caminho_subdiretorio = os.path.join(diretorio_raiz, subdiretorio)

        # Verifica se o item no diretório é um subdiretório
        if os.path.isdir(caminho_subdiretorio):
            # Obtém o nome do subdiretório
            nome_subdiretorio = subdiretorio

            # Renomeia os arquivos dentro do subdiretório
            for arquivo in os.listdir(caminho_subdiretorio):
                caminho_arquivo_antigo = os.path.join(caminho_subdiretorio, arquivo)
                nome_arquivo, extensao = os.path.splitext(arquivo)

                # Verifica se o arquivo é o arquivo 'genomic.gbff' ou 'protein.faa'
                if extensao == '.gbff':
                    novo_nome_arquivo = f'{nome_subdiretorio}_genomics.gbff'
                elif extensao == '.faa':
                    novo_nome_arquivo = f'{nome_subdiretorio}_protein.faa'
                else:
                    continue  # Ignora outros tipos de arquivos

                caminho_arquivo_novo = os.path.join(caminho_subdiretorio, novo_nome_arquivo)

                # Renomeia o arquivo
                os.rename(caminho_arquivo_antigo, caminho_arquivo_novo)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Script para renomear arquivos em subdiretórios.')
    parser.add_argument('diretorio_raiz', help='Diretório raiz que contém os subdiretórios a serem processados')

    args = parser.parse_args()
    rename_files(args.diretorio_raiz)
