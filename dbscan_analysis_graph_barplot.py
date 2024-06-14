import matplotlib.ticker as plticker
import os
import pandas as pd
import argparse
from collections import defaultdict, Counter
import seaborn as sns
import matplotlib.pyplot as plt
import re

# Author: Leandro de Mattos Pereira
# CNP team - Leao Laboratory
# 11/10/2023

def format_species_name(species_name):
    # Substitui o '_' por ' ' e coloca em itálico o nome das espécies que precedem o sufixo LEGE
    parts = species_name.split('_')
    new_parts = []
    italic = True
    for part in parts:
        if italic:
            if "LEGE" in part:
                italic = False
                new_parts.append(part)
            else:
                new_parts.append(f'$\it{{{part}}}$')
        else:
            new_parts.append(part)
    return ' '.join(new_parts).replace('{ ', '{').replace(' }', '}')

def collect_data_from_overviews(directory):
    data_diamond = defaultdict(Counter)
    data_dbCAN_sub = defaultdict(Counter)
    data_HMMER = defaultdict(Counter)
    species_counts = defaultdict(lambda: {'HMMER': 0, 'dbCAN_sub': 0, 'DIAMOND': 0})

    for subdir, _, files in os.walk(directory):
        for file in files:
            if "overview_output_" in file:
                species_name = file.split("overview_output_")[1].rsplit(".txt", 1)[0]
                formatted_species_name = format_species_name(species_name)
                file_path = os.path.join(subdir, file)
                df = pd.read_csv(file_path, sep="\t", usecols=["HMMER", "dbCAN_sub", "DIAMOND"])

                for _, row in df.iterrows():
                    if row["DIAMOND"] != "-":
                        data_diamond[formatted_species_name][row["DIAMOND"]] += 1
                        species_counts[formatted_species_name]['DIAMOND'] += 1

                    if row["dbCAN_sub"] != "-":
                        dbCAN_code = row["dbCAN_sub"].split("_")[0]
                        data_dbCAN_sub[formatted_species_name][dbCAN_code] += 1
                        species_counts[formatted_species_name]['dbCAN_sub'] += 1

                    if row["HMMER"] != "-":
                        HMMER_code = row["HMMER"].split("(")[0]
                        data_HMMER[formatted_species_name][HMMER_code] += 1
                        species_counts[formatted_species_name]['HMMER'] += 1

    return data_diamond, data_dbCAN_sub, data_HMMER, species_counts

def plot_stacked_bar(data, title, file_prefix, top_n=40):
    df = pd.DataFrame(data).fillna(0)
    df = df.T
    
    top_n_subfamilies = df.sum().sort_values(ascending=False).head(top_n).index
    df = df[top_n_subfamilies]
    
    # Ordenar as colunas alfanumericamente
    df = df[sorted(df.columns)]
    
    # Usar a paleta de cores personalizada
    color_palette = [
        '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
        '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
        '#2E86C1', '#C0392B', '#F39C12', '#16A085', '#8E44AD',
        '#58D68D', '#6C3483', '#B03A2E', '#117A65', '#C5E17A',
        '#7D6608', '#76D7C4', '#6E2C00', '#1B4F72', '#641E16',
        '#DFFF00', '#FFBF00', '#FF7F50', '#DE3163', '#9FE2BF',
        '#40E0D0', '#6495ED', '#CCCCFF', '#DFFF00', '#00FFFF',
        '#5A9', '#B5B', '#A52A2A', '#7FFF00', '#FA8072'
    ]
    
    plt.figure(figsize=(12, 20))
    
    ax = df.plot(kind="barh", stacked=True, color=color_palette)
    ax.set_title(title)
    ax.set_xlabel("Counts", fontsize=10)
    ax.set_ylabel("Species", fontsize=6)
    
    plt.xticks(rotation=0, fontsize=8)
    
    loc = plticker.MultipleLocator(base=50)
    ax.xaxis.set_major_locator(loc)
    
    plt.legend(title="Subfamily", loc="center left", bbox_to_anchor=(1.0, 0.5), ncol=2, fontsize=8)
    
    plt.subplots_adjust(left=0.2, right=0.7)
    plt.yticks(fontsize=6)
    
    plt.savefig(f"{file_prefix}_{title.replace(' ', '_')}.png", dpi=300, bbox_inches="tight")
    plt.savefig(f"{file_prefix}_{title.replace(' ', '_')}.svg", format="svg", bbox_inches="tight")
    plt.close()

def save_species_counts_table(species_counts, output_file):
    df = pd.DataFrame.from_dict(species_counts, orient='index')
    df.to_csv(output_file, index_label='Species', sep='\t')

def save_detailed_counts_table(data, output_file):
    df = pd.DataFrame(data).fillna(0)
    df.to_csv(output_file, sep='\t')

def main():
    parser = argparse.ArgumentParser(description="Generate stacked bar charts based on overview files.")
    parser.add_argument("directory", help="Directory containing subdirectories with overview_output files.")
    args = parser.parse_args()

    data_diamond, data_dbCAN_sub, data_HMMER, species_counts = collect_data_from_overviews(args.directory)
    
    plot_stacked_bar(data_diamond, 'DIAMOND Subfamily Counts for Different Species', 'barplot')
    plot_stacked_bar(data_dbCAN_sub, 'dbCAN_sub Subfamily Counts for Different Species', 'barplot')
    plot_stacked_bar(data_HMMER, 'HMMER Subfamily Counts for Different Species', 'barplot')
    
    save_species_counts_table(species_counts, 'species_counts_table.tsv')
    save_detailed_counts_table(data_diamond, 'detailed_diamond_counts.tsv')
    save_detailed_counts_table(data_dbCAN_sub, 'detailed_dbCAN_sub_counts.tsv')
    save_detailed_counts_table(data_HMMER, 'detailed_HMMER_counts.tsv')

if __name__ == "__main__":
    main()


