import pandas as pd
import json
import sys
import argparse
#Author: Leandro de Mattos Pereira
# CNP -Team
#Team Leader. Pedro Leão
def process_text_to_dataframe(text):
    # Divide o texto em linhas
    lines = text.strip().split("\n")
    
    # Divide cada linha em campos usando tabulaÃ§Ã£o como delimitador
    rows = [line.split("\t") for line in lines]
    
    # Define the expected keys in the sequence_data dictionary
    expected_keys = ['value', 'length', 'molWeight', 'crc64', 'md5']
    
    # Processa a coluna Sequence para extrair os campos internos
    headers = rows[0] + expected_keys  # Add headers
    new_rows = [headers]
    
    for row in rows[1:]:  # Ignorando a primeira linha (cabeÃ§alho)
        try:
            sequence_data = json.loads(row[11].replace("'", '"'))
            # Extract sequence data and append it to the row
            new_row = row[:11] + [sequence_data.get(key, "") for key in expected_keys] + row[12:]
            
            # Ensure the row has the correct number of columns by filling in missing columns with empty strings
            new_row += [''] * (len(headers) - len(new_row))
            
            new_rows.append(new_row)
        except (IndexError, json.JSONDecodeError):
            print(f"WARNING: Problem processing row: {row}")
            continue
    
    # Cria um DataFrame a partir dos dados extraÃ­dos
    df = pd.DataFrame(new_rows[1:], columns=new_rows[0])
    
    # Renomeia a coluna "value" para "Sequence"
    df.rename(columns={"value": "Sequence"}, inplace=True)
    
    return df

def main():
    parser = argparse.ArgumentParser(description="Processa um arquivo de texto e gera arquivos Excel e TSV.")
    parser.add_argument("input_file", help="Caminho para o arquivo de entrada.")
    parser.add_argument("--excel", default="output.xlsx", help="Nome do arquivo Excel de saÃ­da. PadrÃ£o: output.xlsx.")
    parser.add_argument("--tsv", default="output.tsv", help="Nome do arquivo TSV de saÃ­da. PadrÃ£o: output.tsv.")
    
    args = parser.parse_args()
    
    with open(args.input_file, 'r') as f:
        text = f.read()
    
    df = process_text_to_dataframe(text)
    
    # Salva o DataFrame em Excel e TSV
    df.to_excel(args.excel, index=False)
    df.to_csv(args.tsv, sep="\t", index=False)
    print(f"Arquivos {args.excel} e {args.tsv} gerados com sucesso!")

if __name__ == "__main__":
    main()
