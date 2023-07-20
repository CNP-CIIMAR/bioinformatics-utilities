import sys

def filter_table(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    filtered_lines = []
    organism_names = {}

    for line in lines:
        line = line.strip()
        if line:  # Ignorar linhas vazias
            columns = line.split('\t')
            if len(columns) > 2:  # Verificar se há colunas suficientes
                accession = columns[0]
                organism_name = columns[2]

                if organism_name not in organism_names:
                    organism_names[organism_name] = []

                if accession.startswith('GCF_'):
                    organism_names[organism_name].append(accession)
                elif accession.startswith('GCA_') and len(organism_names[organism_name]) == 0:
                    organism_names[organism_name].append(accession)

    for line in lines:
        line = line.strip()
        if line:  # Ignorar linhas vazias
            columns = line.split('\t')
            if len(columns) > 2:  # Verificar se há colunas suficientes
                accession = columns[0]
                organism_name = columns[2]

                if accession in organism_names.get(organism_name, []):
                    filtered_lines.append(line)

    filtered_table = '\n'.join(filtered_lines)
    return filtered_table

def print_help():
    print("Uso: python script.py <caminho_para_tabela> <caminho_para_output>")
    print("Filtre a tabela com base nos números de acesso de organismos e gere uma nova tabela de saída.")
    print("Certifique-se de fornecer o caminho para o arquivo contendo a tabela como o primeiro argumento e o caminho para o arquivo de saída como o segundo argumento.")

# Verificar se foram fornecidos os caminhos da tabela e do arquivo de saída como argumentos
if len(sys.argv) < 3:
    print_help()
    sys.exit(1)

table_path = sys.argv[1]  # Caminho para o arquivo contendo a tabela
output_path = sys.argv[2]  # Caminho para o arquivo de saída

filtered_table = filter_table(table_path)

# Escrever a tabela filtrada no arquivo de saída
with open(output_path, 'w') as output_file:
    output_file.write(filtered_table)

print("Tabela filtrada foi escrita em", output_path)

