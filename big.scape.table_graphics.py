## Authors: Leandro de Mattos Pereira
## CNP -team - Dr. Pedro Leao, Team Leader.

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def plot_big_scape_class_by_order(input_file):
    # Carregar o arquivo CSV
    df = pd.read_csv(input_file)

    # Extrair a coluna "Taxonomy" e dividir em diferentes níveis taxonômicos
    taxonomic_levels = df['Taxonomy'].str.split('; ')

    # Obter a Order a partir da sétima posição
    order_values = taxonomic_levels.apply(lambda x: x[6] if len(x) >= 7 else '')

    # Adicionar a coluna "Order" ao DataFrame
    df['Order'] = order_values

    # Contar a quantidade de ocorrências de cada classe de BGC por Order
    class_counts = df.groupby(['Order', 'BiG-SCAPE class']).size().unstack()

    # Configurar a paleta de cores para as Orders
    palette = sns.color_palette('Set3', len(class_counts.columns))

    # Configurar o estilo do gráfico
    sns.set(style='whitegrid')

    # Gerar o gráfico de barras empilhadas
    fig, ax = plt.subplots(figsize=(10, 8))
    class_counts.plot(kind='barh', stacked=True, ax=ax, color=palette)

    # Personalizar o gráfico
    ax.set_xlabel('Quantidade', fontsize=12)
    ax.set_ylabel('Order', fontsize=12)
    ax.set_title('Distribuição de BGC Classes por Order', fontsize=14, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(title='BiG-SCAPE Class')

    # Ajustar a posição da legenda
    ax.legend(title='BiG-SCAPE Class', loc='center left', bbox_to_anchor=(1, 0.5), fontsize='small')

    plt.tight_layout()

    # Exibir o gráfico na tela
    plt.show()

    # Salvar o gráfico em formatos de alta qualidade
    fig.savefig('big_scape_class_distribution_graph.svg', format='svg', dpi=300)
    fig.savefig('big_scape_class_distribution_graph.png', format='png', dpi=300)
    fig.savefig('big_scape_class_distribution_graph.jpeg', format='jpeg', dpi=300)

if __name__ == '__main__':
    # Configurar o argumento de linha de comando
    parser = argparse.ArgumentParser(description='Gerar um gráfico de barras empilhadas da distribuição de BGC Classes por Order.')
    parser.add_argument('input_file', help='Caminho para o arquivo CSV de entrada.')

    # Parse dos argumentos
    args = parser.parse_args()

    # Gerar o gráfico
    plot_big_scape_class_by_order(args.input_file)
