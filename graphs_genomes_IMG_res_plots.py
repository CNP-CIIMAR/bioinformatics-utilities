# Função para gerar os gráficos e salvá-los em alta resolução nos formatos PNG, SVG e JPEG
def generate_plots(df, output_folder):
    numeric_columns = df.select_dtypes(include=[float, int]).columns.drop(["taxon_oid", "Taxon Object ID", "Proposal GOLD ID"])

    for column in numeric_columns:
        plt.figure(figsize=(10, 5))
        plt.title(column)
        ax = plt.gca()

        # Calcular a média e o desvio padrão para cada combinação de lagoas e pontos de coleta
        means = df.groupby(["Lakes sediments", "Sample"])[column].mean()
        stds = df.groupby(["Lakes sediments", "Sample"])[column].std()
        labels = [f"{l} {p}" for l, p in means.index]

        # Plotar os gráficos de boxplot
        df.boxplot(column=column, by=["Lakes sediments", "Sample"], ax=ax)

        # Adicionar as barras de erro representando o desvio padrão
        xticklabels = ax.get_xticklabels()
        xticks = ax.get_xticks()
        for label, tick, mean, std in zip(labels, xticks, means, stds):
            plt.errorbar(tick, mean, yerr=std, fmt='o', color='red', markersize=5, capsize=4)

        plt.suptitle("")  # Remover o título do subplot gerado pelo pandas
        plt.xlabel("Lagoas e pontos de coleta")
        plt.ylabel(column)
        plt.xticks(ticks=xticks, labels=labels, rotation=45)
        plt.tight_layout()

        # Salvar o gráfico em PNG, SVG e JPEG
        output_file = f"{output_folder}/{column}."
        plt.savefig(output_file + "png", format='png', dpi=300)
        plt.savefig(output_file + "svg", format='svg')
        plt.savefig(output_file + "jpeg", format='jpeg', dpi=300)

        plt.close()

# Função para ler os dados da tabela TSV
def read_data(file_path):
    return pd.read_csv(file_path, sep="\t")

# Função principal para processar os argumentos da linha de comando
def main():
    parser = argparse.ArgumentParser(description="Script para gerar gráficos boxplot a partir de uma tabela TSV.")
    parser.add_argument("input_file", help="Caminho para o arquivo TSV contendo os dados.")
    parser.add_argument("output_folder", help="Pasta onde os gráficos serão salvos.")
    args = parser.parse_args()

    # Ler os dados do arquivo TSV
    df = read_data(args.input_file)

    # Gerar e salvar os gráficos em alta resolução
    generate_plots(df, args.output_folder)

if __name__ == "__main__":
    main()
