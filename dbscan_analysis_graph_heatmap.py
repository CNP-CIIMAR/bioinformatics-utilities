import os
import pandas as pd
import argparse
from collections import defaultdict, Counter
import seaborn as sns
import matplotlib.pyplot as plt

def collect_data_from_overviews(directory):
    data_diamond = defaultdict(Counter)
    data_dbCAN_sub = defaultdict(Counter)
    data_HMMER = defaultdict(Counter)

    for subdir, _, files in os.walk(directory):
        for file in files:
            if "overview_output_" in file:
                species_name = file.split("overview_output_")[1].rsplit(".txt", 1)[0]
                file_path = os.path.join(subdir, file)
                df = pd.read_csv(file_path, sep="\t", usecols=["HMMER", "dbCAN_sub", "DIAMOND"])

                for _, row in df.iterrows():
                    if row["DIAMOND"] != "-":
                        data_diamond[species_name][row["DIAMOND"]] += 1

                    if row["dbCAN_sub"] != "-":
                        dbCAN_code = row["dbCAN_sub"].split("_")[0]
                        data_dbCAN_sub[species_name][dbCAN_code] += 1

                    if row["HMMER"] != "-":
                        HMMER_code = row["HMMER"].split("(")[0]
                        data_HMMER[species_name][HMMER_code] += 1

    return data_diamond, data_dbCAN_sub, data_HMMER

def plot_data(data, title):
    df = pd.DataFrame(data).fillna(0).T
    # Ordenar as colunas em ordem crescente
    df = df.reindex(sorted(df.columns), axis=1)
    
    plt.figure(figsize=(20, 15))
    ax = sns.heatmap(df, cmap="coolwarm")

    # Ajusta a posiÃ§Ã£o dos ticks para o meio das cÃ©lulas
    ax.set_xticks([x + 0.5 for x in range(df.shape[1])])
    ax.set_yticks([y + 0.5 for y in range(df.shape[0])])
    
    ax.set_xticklabels(df.columns, rotation=90, fontsize=8)  # Rotacionar os labels e ajustar o tamanho da fonte
    ax.set_yticklabels(df.index, fontsize=6)  # Ajustar o tamanho da fonte dos labels no eixo Y

    plt.title(title)
    plt.ylabel('Species')
    plt.xlabel('Subfamily')
    plt.tight_layout()

    # Remove espaÃ§os e apÃ³strofos do tÃ­tulo para nomear o arquivo
    file_name = title.replace(" ", "_").replace("'", "") + ".png"
    plt.savefig(file_name, dpi=300)
    plt.close()
    
def save_count_tables(data, title):
    all_data = []  # Lista para armazenar os dados de todas as tabelas

    for result_type, species_data in data.items():
        for species_name, subfamily_counts in species_data.items():
            if isinstance(subfamily_counts, int):
                # Se a contagem for um valor inteiro, cria uma Ãºnica entrada na tabela
                all_data.append((species_name, subfamily_counts, result_type))
            else:
                # Caso contrÃ¡rio, converte o dicionÃ¡rio Counter em uma lista de tuplas (subfamily, count)
                subfamily_counts_list = [(result_type, count, species_name) for result_type, count in subfamily_counts.items()]
                all_data.extend(subfamily_counts_list)

    # Criar um DataFrame a partir da lista de tuplas
    df = pd.DataFrame(all_data, columns=["Subfamily", "contagem", "Specie"])

    # Remove espaÃ§os e apÃ³strofos do tÃ­tulo para nomear o arquivo
    table_name = title.replace(" ", "_").replace("'", "") + "_Count_Table.tsv"
    df.to_csv(table_name, sep="\t", index=False)  # Especifique sep="\t" e index=False para salvar em formato TSV sem Ã­ndice

        
def main():
    parser = argparse.ArgumentParser(description="Generate heatmaps based on overview files.")
    parser.add_argument("directory", help="Directory containing subdirectories with overview_output files.")
    args = parser.parse_args()

    data_diamond, data_dbCAN_sub, data_HMMER = collect_data_from_overviews(args.directory)
    plot_data(data_diamond, 'DIAMOND Subfamily Counts for Different Species')
    plot_data(data_dbCAN_sub, 'dbCAN_sub Subfamily Counts for Different Species')
    plot_data(data_HMMER, 'HMMER Subfamily Counts for Different Species')
    save_count_tables(data_diamond, 'DIAMOND_Subfamily_Counts_for_Different_Species')
    save_count_tables(data_dbCAN_sub, 'dbCAN_sub_Subfamily_Counts_for_Different_Species')
    save_count_tables(data_HMMER, 'HMMER_Subfamily_Counts_for_Different_Species')

if __name__ == "__main__":
    main()
