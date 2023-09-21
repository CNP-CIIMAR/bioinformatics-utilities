import os
import re
import pandas as pd


def verificar_id_no_header(fasta_id, fasta_file):
    with open(fasta_file, 'r') as f:
        header = ''
        for line in f:
            if line.startswith('>'):
                header = line.strip()[1:]
                if fasta_id in header:
                    # Remover a parte após o caractere '#' e remover o próprio caractere '#'
                    header = re.sub(r'\s*#.*$', '', header).replace('#', '')
                    return header
                else:
                    header = ''
    return ''


def processar_arquivos_fasta(diretorio, lista_ids):
    tabela_resultado = []
    for root, dirs, files in os.walk(diretorio):
        for file in files:
            if file.endswith('.faa'):
                arquivo_fasta = os.path.join(root, file)
                for fasta_id in lista_ids:
                    id_encontrado = verificar_id_no_header(fasta_id, arquivo_fasta)
                    if id_encontrado:
                        nome_arquivo = file
                        especie = re.search(r'\[(.*?)\]', id_encontrado)
                        especie = especie.group(1) if especie else ''
                        tabela_resultado.append([fasta_id, id_encontrado, nome_arquivo, especie])

    return tabela_resultado


def salvar_resultado_em_tabela(tabela_resultado, nome_arquivo_saida):
    df = pd.DataFrame(tabela_resultado, columns=["ID Query", "ID found", "Filename", "Specie"])
    df.to_csv(nome_arquivo_saida, sep='\t', index=False)


def main(diretorio, arquivo_ids, nome_arquivo_saida):
    with open(arquivo_ids, 'r') as f:
        lista_ids = [line.strip() for line in f if line.strip()]

    tabela_resultado = processar_arquivos_fasta(diretorio, lista_ids)
    salvar_resultado_em_tabela(tabela_resultado, nome_arquivo_saida)


if __name__ == '__main__':
    import sys

    if len(sys.argv) != 4:
        print("Uso: python script.py <diretorio> <arquivo_ids> <nome_arquivo_saida>")
        sys.exit(1)

    diretorio = sys.argv[1]
    arquivo_ids = sys.argv[2]
    nome_arquivo_saida = sys.argv[3]

    main(diretorio, arquivo_ids, nome_arquivo_saida)
