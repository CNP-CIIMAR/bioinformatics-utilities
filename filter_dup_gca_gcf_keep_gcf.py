import sys

def filter_table(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    filtered_lines = []
    organism_names = {}

    for line in lines:
        line = line.strip()
        columns = line.split('\t')
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
        columns = line.split('\t')
        accession = columns[0]
        organism_name = columns[2]

        if accession in organism_names[organism_name]:
            filtered_lines.append(line)

    filtered_table = '\n'.join(filtered_lines)
    return filtered_table

# Exemplo de uso:
table_path = sys.argv[1]  # Caminho para o arquivo contendo a tabela
filtered_table = filter_table(table_path)
print(filtered_table)
