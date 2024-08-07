import pandas as pd
import sys

def read_proteins_list(file_path):
    with open(file_path, 'r') as file:
        proteins_list = [line.strip() for line in file]
    return proteins_list

def create_output_table(proteins_list, dataframe):
    # Criar um dataframe vazio com as mesmas colunas do dataframe original
    columns = dataframe.columns.tolist()
    output_data = {col: [] for col in columns}

    for protein in proteins_list:
        if protein in dataframe["Protein Accession"].values:
            match = dataframe[dataframe["Protein Accession"] == protein]
            for col in columns:
                output_data[col].append(match[col].values[0])
        else:
            for col in columns:
                output_data[col].append("not genome found")
    
    output_df = pd.DataFrame(output_data)
    return output_df

def main(input_file, table_file, output_file):
    proteins_list = read_proteins_list(input_file)
    dataframe = pd.read_csv(table_file, sep='\t')

    output_df = create_output_table(proteins_list, dataframe)
    output_df.to_csv(output_file, index=False, sep='\t')

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Uso: python script.py <input_file> <table_file> <output_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    table_file = sys.argv[2]
    output_file = sys.argv[3]

    main(input_file, table_file, output_file)
