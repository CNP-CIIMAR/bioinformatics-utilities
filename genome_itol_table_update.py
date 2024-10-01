import pandas as pd
import argparse
import sys
import os

def extrair_prefixo(genome_id):
    """
    Extrai o prefixo antes do segundo '_'.
    Exemplo: 'GCF_030382115.1_ASM3038211v1_genomic' -> 'GCF_030382115.1'
    """
    partes = genome_id.split('_', 2)  # Divide no mÃ¡ximo em 3 partes
    if len(partes) >= 2:
        return f"{partes[0]}_{partes[1]}"
    return genome_id  # Retorna o original se nÃ£o houver '_'

def limpar_id(prefixo):
    """
    Remove qualquer sufixo apÃ³s o ponto '.'.
    Exemplo: 'GCF_030382115.1' -> 'GCF_030382115'
    """
    return prefixo.split('.')[0]

def atualizar_tabela(genome_ids_path, tabela_path, tabela_saida_path):
    """
    Atualiza a tabela adicionando as colunas 'Label' e 'Color' com base na correspondÃªncia dos Genome IDs.
    """
    # Verificar se os arquivos de entrada existem
    if not os.path.isfile(genome_ids_path):
        print(f"Erro: O arquivo '{genome_ids_path}' nÃ£o foi encontrado.")
        sys.exit(1)
    if not os.path.isfile(tabela_path):
        print(f"Erro: O arquivo '{tabela_path}' nÃ£o foi encontrado.")
        sys.exit(1)

    print("Carregando a lista de IDs de genomas...")
    try:
        # Carregar a lista de IDs de genomas
        with open(genome_ids_path, 'r') as file:
            genome_ids = set(line.strip() for line in file if line.strip())
    except Exception as e:
        print(f"Erro ao ler '{genome_ids_path}': {e}")
        sys.exit(1)

    print("Carregando a tabela de entrada...")
    try:
        # Carregar a tabela
        tabela = pd.read_csv(tabela_path, sep='\t')  # Ajuste o separador se necessÃ¡rio
    except Exception as e:
        print(f"Erro ao ler '{tabela_path}': {e}")
        sys.exit(1)

    # Verificar se a coluna 'Tree node ID' existe
    if 'Tree node ID' not in tabela.columns:
        print("Erro: A coluna 'Tree node ID' nÃ£o foi encontrada na tabela.")
        sys.exit(1)

    print("Processando a tabela...")
    # Extrair e limpar os prefixos
    tabela['Prefixo'] = tabela['Tree node ID'].apply(extrair_prefixo)
    tabela['Prefixo_Limpo'] = tabela['Prefixo'].apply(limpar_id)

    # Verificar correspondÃªncia
    tabela['Match'] = tabela['Prefixo_Limpo'].isin(genome_ids)

    # Preencher 'Label' e 'Color' usando operaÃ§Ãµes vetorizadas
    tabela.loc[tabela['Match'], 'Label'] = 'AMP-L-CG'
    tabela.loc[tabela['Match'], 'Color'] = '#FFA500'  # Hexadecimal para laranja

    # Opcional: Remover as colunas auxiliares
    tabela.drop(['Prefixo', 'Prefixo_Limpo', 'Match'], axis=1, inplace=True)

    print("Salvando a tabela atualizada...")
    try:
        # Salvar a tabela atualizada
        tabela.to_csv(tabela_saida_path, sep='\t', index=False)
        print(f"Tabela atualizada salva em '{tabela_saida_path}'.")
    except Exception as e:
        print(f"Erro ao salvar a tabela atualizada em '{tabela_saida_path}': {e}")
        sys.exit(1)

def main():
    # Configurar o parser de argumentos
    parser = argparse.ArgumentParser(description="Atualiza uma tabela com base em uma lista de Genome IDs.")
    parser.add_argument('-i', '--input_ids', required=True, help="Caminho para o arquivo de IDs de genomas (e.g., genome_ids.txt).")
    parser.add_argument('-t', '--input_table', required=True, help="Caminho para o arquivo da tabela de entrada (e.g., tabela2.tsv).")
    parser.add_argument('-o', '--output_table', required=True, help="Caminho para o arquivo da tabela de saÃ­da (e.g., tabela2_atualizada.tsv).")

    args = parser.parse_args()

    # Exibir os caminhos dos arquivos para verificaÃ§Ã£o
    print(f"Arquivo de IDs de genomas: {args.input_ids}")
    print(f"Tabela de entrada: {args.input_table}")
    print(f"Tabela de saÃ­da: {args.output_table}")

    # Chamar a funÃ§Ã£o para atualizar a tabela
    atualizar_tabela(args.input_ids, args.input_table, args.output_table)

if __name__ == "__main__":
    main()

