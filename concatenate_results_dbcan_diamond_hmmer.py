import os
import pandas as pd
import argparse

def concatenate_tsv_files(root_dir):
    all_data = []
    
    for root, dirs, files in os.walk(root_dir):
        for file in files:
            if file == "hmmer.out":
                file_path = os.path.join(root, file)
                try:
                    data = pd.read_csv(file_path, sep="\t", header=None, skiprows=1)
                    
                    # Verificar se a linha contém o cabeçalho e, neste caso, descartá-la
                    data = data[~data[0].str.contains("HMM Profile", na=False)]
                    
                    _, dir_name = os.path.split(root)
                    specie_name = dir_name.replace("output_", "")
                    data['Specie'] = specie_name  # Adicionando a coluna 'Specie'
                    all_data.append(data)
                except Exception as e:
                    print(f"Erro ao ler o arquivo {file_path}: {str(e)}")
    
    if all_data:
        final_data = pd.concat(all_data, ignore_index=True)
        
        final_data.to_csv("final_output_hmmer.tsv", sep="\t", index=False, header=[
            "HMM Profile", "Profile Length", "Gene ID", "Gene Length", "E Value",
            "Profile Start", "Profile End", "Gene Start", "Gene End", "Coverage", "Specie"
        ])
    else:
        print("Nenhum dado foi encontrado para concatenação.")

def main():
    parser = argparse.ArgumentParser(description='Concatenate TSV files.')
    parser.add_argument('root_dir', type=str, help='Root directory containing the subdirectories with hmmer.out files.')
    args = parser.parse_args()
    concatenate_tsv_files(args.root_dir)

if __name__ == "__main__":
    main()

