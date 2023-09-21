import pandas as pd
import sys

def process_tables(file2, output_file, count_file, *files1):
    # Carregar a segunda tabela
    table2 = pd.read_csv(file2, sep='\t')

    # Criar uma cópia da segunda tabela
    new_table = table2.copy()

    # Lista para armazenar os nomes das colunas 'GENE'
    gene_columns = []

    # Processar cada tabela de entrada
    for i, file1 in enumerate(files1, start=1):
        # Carregar a tabela de entrada
        table1 = pd.read_csv(file1, sep='\t')

        # Adicionar a coluna 'GENE' à nova tabela
        gene_column = 'GENE_fad' + chr(ord('A') + i - 1)  # 'GENE_fadA', 'GENE_fadB', etc.
        new_table[gene_column] = ''
        gene_columns.append(gene_column)

        # Iterar sobre cada linha da primeira tabela
        for filename in table1['Filename']:
            # Remover o sufixo "_protein.faa"
            assembly = filename.replace('_protein.faa', '')

            # Procurar por uma correspondência na coluna 'Assembly' da nova tabela
            matches = new_table['Assembly'].str.contains(assembly)

            # Preencher a célula correspondente na coluna 'GENE' com 'fadA' se houver correspondência
            new_table.loc[matches, gene_column] = 'fad' + chr(ord('A') + i - 1)

    # Salvar a nova tabela
    new_table.to_csv(output_file, sep='\t', index=False)

    # Contar o número de correspondências para cada tabela de entrada
    match_counts = new_table[gene_columns].apply(lambda col: col[col != ''].count())

    # Identificar as linhas onde todas as colunas 'GENE' são preenchidas
    common_matches = new_table[gene_columns].dropna().shape[0]

    # Imprimir as contagens em um arquivo de saída
    with open(count_file, 'w') as f:
        for col, count in match_counts.items():
            f.write(f'{col}: {count}\n')
        f.write(f'Common matches: {common_matches}\n')

def print_help():
    print("""
Usage: python script.py TABLE2 OUTPUT COUNTS TABLE1...

This script processes one or more TABLE1-style input tables against a TABLE2-style input table.
For each TABLE1, it adds a column to TABLE2, which is filled with 'fadX' where there is a match
between 'Filename' in TABLE1 and 'Assembly' in TABLE2. 'X' changes for each TABLE1 (fadA for the first, fadB for the second, etc.).

Arguments:
    TABLE2:  The name of the TABLE2-style input file. This table should be a TSV file with a 'Assembly' column.
    OUTPUT:  The name of the output file. The output table will be a TSV file.
    COUNTS:  The name of the counts output file. The counts of matches for each gene column and the number of common matches will be written to this file.
    TABLE1:  One or more names of TABLE1-style input files. These tables should be TSV files with a 'Filename' column.

Example:
    python script.py table2.txt output.txt counts.txt table1a.txt table1b.txt
""")

# Check for help argument
if '-h' in sys.argv or '--help' in sys.argv:
    print_help()
    sys.exit()

# Check for correct number of arguments
if len(sys.argv) < 5:
    print("Error: too few arguments")
    print_help()
    sys.exit()

# Recuperar os nomes dos arquivos a partir dos argumentos do script
file2 = sys.argv[1]
output_file = sys.argv[2]
count_file = sys.argv[3]
files1 = sys.argv[4:]

# Processar as tabelas
process_tables(file2, output_file, count_file, *files1)

