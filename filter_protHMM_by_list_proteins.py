import pandas as pd
import argparse

def clean_text(text):
    """Limpa o texto substituindo ou removendo caracteres problemÃ¡ticos, como quebras de linha."""
    if isinstance(text, str):
        return text.replace('\n', ' ').replace('\r', ' ').strip()
    return text

def filter_protein_ids(protein_id_file, input_table_path, output_table_path, delimiter='\t'):
    # Ler a lista de IDs de proteÃ­nas de um arquivo
    with open(protein_id_file, 'r') as file:
        protein_id_list = [line.strip() for line in file.read().splitlines() if line.strip()]

    print("Primeiros 5 IDs de proteÃ­nas carregados:", protein_id_list[:5])  # Para verificar os IDs carregados

    # Tentar carregar a tabela de dados
    try:
        data = pd.read_csv(input_table_path, delimiter=delimiter)
    except pd.errors.ParserError as e:
        print(f"Erro ao ler o arquivo CSV: {e}")
        return

    # Limpar dados antes de filtrar
    for col in data.columns:
        data[col] = data[col].apply(clean_text)

    # Filtrar a tabela para incluir apenas linhas cujo target_name estÃ¡ na lista de IDs de proteÃ­nas
    filtered_data = data[data['target_name'].isin(protein_id_list)]
    
    # Verificar se hÃ¡ dados apÃ³s filtrar
    print("NÃºmero de linhas apÃ³s filtrar:", len(filtered_data))

    # Salvar a tabela filtrada se houver linhas filtradas
    if len(filtered_data) > 0:
        filtered_data.to_csv(output_table_path, index=False, sep='\t')
        print(f"Tabela filtrada salva como {output_table_path}")
    else:
        print("Nenhuma correspondÃªncia encontrada. Verifique a lista de IDs e o arquivo de entrada.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filtrar IDs de proteÃ­nas em uma tabela.")
    parser.add_argument("protein_id_file", help="Caminho para o arquivo contendo os IDs de proteÃ­nas")
    parser.add_argument("input_table_path", help="Caminho para o arquivo CSV de entrada")
    parser.add_argument("output_table_path", help="Caminho para o arquivo CSV de saÃ­da")
    parser.add_argument("--delimiter", help="Delimitador usado no arquivo CSV", default='\t')

    args = parser.parse_args()
    
    filter_protein_ids(args.protein_id_file, args.input_table_path, args.output_table_path, args.delimiter)

