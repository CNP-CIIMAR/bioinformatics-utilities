import os
import argparse
import re

def ler_enzimas(arquivo_enzimas):
    with open(arquivo_enzimas, 'r', encoding='utf-8') as f:
        enzimas = [linha.strip() for linha in f.readlines()]
    return enzimas

def extrair_nome_specie(nome_arquivo):
    match = re.search(r'Result_(.*?)\.fasta', nome_arquivo)
    if match:
        return match.group(1)
    return "Desconhecido"

def buscar_enzimas(diretorio, enzimas):
    # Lista para armazenar os resultados
    resultados = []
    
    # Percorrer todos os subdiretÃ³rios e arquivos
    for root, _, files in os.walk(diretorio):
        for file in files:
            if file.endswith(".gbk"):
                file_path = os.path.join(root, file)
                especie = extrair_nome_specie(file_path)
                print(f"Lendo arquivo: {file_path}")
                try:
                    with open(file_path, 'r', encoding='utf-8') as f:
                        conteudo = f.read()
                        for enzima in enzimas:
                            comprimento_minimo = len(enzima) // 2
                            encontrado = False
                            for i in range(len(enzima) - comprimento_minimo + 1):
                                substring = enzima[i:i + comprimento_minimo]
                                if substring in conteudo:
                                    encontrado = True
                                    # Encontrar a anotaÃ§Ã£o correspondente
                                    anotacoes = re.findall(r'/description="(.*?)"', conteudo)
                                    for anotacao in anotacoes:
                                        if substring in anotacao:
                                            resultados.append(f"{especie} | {file_path} | {enzima} | {anotacao}")
                                            break
                                    break
                except Exception as e:
                    print(f"Erro ao ler o arquivo {file_path}: {e}")
    
    return resultados

def salvar_resultados(resultados, arquivo_saida):
    with open(arquivo_saida, 'w', encoding='utf-8') as f:
        for resultado in resultados:
            f.write(resultado + '\n')
    print(f"Resultados salvos em: {arquivo_saida}")

def main():
    parser = argparse.ArgumentParser(description='Buscar arquivos .gbk que contÃªm nomes de enzimas.')
    parser.add_argument('diretorio', type=str, help='DiretÃ³rio raiz para buscar arquivos .gbk')
    parser.add_argument('arquivo_enzimas', type=str, help='Arquivo txt contendo a lista de enzimas')
    parser.add_argument('arquivo_saida', type=str, help='Arquivo de saÃ­da para salvar os resultados')
    
    args = parser.parse_args()
    
    diretorio = args.diretorio
    arquivo_enzimas = args.arquivo_enzimas
    arquivo_saida = args.arquivo_saida
    
    enzimas = ler_enzimas(arquivo_enzimas)
    
    print(f"Buscando no diretÃ³rio: {diretorio}")
    print(f"Procurando pelas enzimas listadas em: {arquivo_enzimas}")
    
    resultados = buscar_enzimas(diretorio, enzimas)
    
    if resultados:
        print("Resultados encontrados:")
        for resultado in resultados:
            print(resultado)
        salvar_resultados(resultados, arquivo_saida)
    else:
        print("Nenhum arquivo contendo as enzimas especificadas foi encontrado.")

if __name__ == "__main__":
    main()
