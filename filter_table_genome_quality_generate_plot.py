import pandas as pd
import argparse
import matplotlib.pyplot as plt

def filter_genomes(input_file, output_file, completeness_threshold, contamination_threshold):
    # Ler a tabela a partir do arquivo CSV
    df = pd.read_csv(input_file, delimiter='\t')
    
    # Filtrar os genomas com base nos thresholds fornecidos
    filtered_df = df[(df['Completeness'] >= completeness_threshold) & (df['Contamination'] <= contamination_threshold)]
    
    # Salvar a tabela filtrada em um novo arquivo CSV
    filtered_df.to_csv(output_file, sep='\t', index=False)
    
    print(f'Tabela filtrada salva em {output_file}')
    
    return filtered_df

def classify_genomes(df, completeness_threshold, contamination_threshold):
    high_quality = df[(df['Completeness'] > 90) & (df['Contamination'] < 5)]
    medium_quality = df[(df['Completeness'] >= 50) & (df['Completeness'] <= 90) & (df['Contamination'] < 10)]
    low_quality = df[(df['Completeness'] < 50) & (df['Contamination'] < 10)]
    min_quality = df[(df['Completeness'] >= completeness_threshold) & (df['Contamination'] <= contamination_threshold)]
    
    return len(high_quality), len(medium_quality), len(low_quality), len(min_quality)

def plot_genome_quality(high_quality, medium_quality, low_quality, min_quality, base_filename, completeness_threshold, contamination_threshold):
    labels = ['High-quality draft\n(>90% complete, <5% contamination)',
              'Medium-quality draft\n(â‰¥50% complete, <10% contamination)',
              'Low-quality drafts\n(<50% complete, <10% contamination)',
              f'Genomes with\nâ‰¥{completeness_threshold}% complete, â‰¤{contamination_threshold}% contamination\n*criterium used by researcher']
    counts = [high_quality, medium_quality, low_quality, min_quality]
    
    # Ordenar as barras para que a barra com maior nÃºmero de genomas apareÃ§a primeiro
    sorted_indices = sorted(range(len(counts)), key=lambda k: counts[k], reverse=True)
    sorted_labels = [labels[i] for i in sorted_indices]
    sorted_counts = [counts[i] for i in sorted_indices]
    
    plt.figure(figsize=(12, 8))
    plt.bar(sorted_labels, sorted_counts, color=['green', 'orange', 'red', 'blue'])
    plt.xlabel('Genome Quality')
    plt.ylabel('Number of Genomes/MAGs')
    plt.title('Mandatory standard metrics according to Genomic Standards Consortium (GSC)')
    
    # Salvar o grÃ¡fico em diferentes formatos
    plt.savefig(f'{base_filename}.png', format='png', dpi=300)
    plt.savefig(f'{base_filename}.svg', format='svg', dpi=300)
    plt.savefig(f'{base_filename}.jpeg', format='jpeg', dpi=300)
    plt.show()

def main():
    # Configurar o argparse para receber argumentos da linha de comando
    parser = argparse.ArgumentParser(description='Filtrar genomas com base em Completeness e Contamination.')
    parser.add_argument('input_file', type=str, help='Arquivo CSV de entrada contendo a tabela de genomas.')
    parser.add_argument('output_file', type=str, help='Arquivo CSV de saÃ­da para salvar a tabela filtrada.')
    parser.add_argument('base_filename', type=str, help='Nome base do arquivo para salvar os grÃ¡ficos.')
    parser.add_argument('completeness_threshold', type=float, help='Valor mÃ­nimo de Completeness.')
    parser.add_argument('contamination_threshold', type=float, help='Valor mÃ¡ximo de Contamination.')

    args = parser.parse_args()
    
    # Filtrar os genomas
    filtered_genomes = filter_genomes(args.input_file, args.output_file, args.completeness_threshold, args.contamination_threshold)
    
    # Classificar os genomas
    high_quality, medium_quality, low_quality, min_quality = classify_genomes(filtered_genomes, args.completeness_threshold, args.contamination_threshold)
    
    # Plotar o grÃ¡fico
    plot_genome_quality(high_quality, medium_quality, low_quality, min_quality, args.base_filename, args.completeness_threshold, args.contamination_threshold)

if __name__ == "__main__":
    main()
