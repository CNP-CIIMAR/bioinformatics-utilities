import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

if len(sys.argv) != 4:
    print("Usage: python script.py <table_file> <output_file> <species_file>")
    sys.exit(1)

# Tabela do InterProScan
table_file = sys.argv[1]
# Nome do arquivo de saída
output_file = sys.argv[2]
# Tabela contendo as colunas Protein.accession e species
species_file = sys.argv[3]

# Leitura da tabela como um DataFrame
df = pd.read_csv(table_file, sep='\t')

# Agrupamento dos dados por Protein.accession
grouped = df.groupby('Protein.accession')

# Dicionário de cores por Signature.description
color_mapping = {}

# Função para obter uma cor única para cada Signature.description
def get_color(signature_description):
    if signature_description not in color_mapping:
        color_mapping[signature_description] = "#%06x" % (hash(signature_description) % 0xFFFFFF)
    return color_mapping[signature_description]

# Dicionário de shapes disponíveis
shapes = {
    0: 'RE',  # Rectangle
    1: 'HH',  # Horizontal Hexagon
    2: 'HV',  # Vertical Hexagon
    3: 'EL',  # Ellipse
    4: 'DI',  # Rhombus (Diamond)
    5: 'TR',  # Right Pointing Triangle
    6: 'TL',  # Left Pointing Triangle
    7: 'PL',  # Left Pointing Pentagram
    8: 'PR',  # Right Pointing Pentagram
    9: 'PU',  # Up Pointing Pentagram
    10: 'PD',  # Down Pointing Pentagram
    11: 'OC',  # Octagon
    12: 'GP',  # Rectangle (Gap; Black filled rectangle with 1/3 normal height)
}

# Função para selecionar a melhor combinação de shapes
def select_shape_combination(lengths):
    sorted_lengths = sorted(lengths)
    selected_shapes = []
    shape_mapping = {}
    shape_counter = 0

    for length in lengths:
        if length not in shape_mapping:
            shape_mapping[length] = shapes[shape_counter % len(shapes)]
            shape_counter += 1
        selected_shapes.append(shape_mapping[length])

    return selected_shapes

# Lista para armazenar as informações das proteínas
protein_info = []
# Dicionário para armazenar a contagem de repetições de cada domínio
domain_count = {}

# Leitura da tabela de espécies
species_df = pd.read_csv(species_file, sep='\t')

# Função para obter o nome da espécie para uma dada Protein.accession
def get_species(protein_accession):
    species = species_df.loc[species_df['Protein.accession'] == protein_accession, 'species']
    if len(species) > 0:
        return species.values[0]
    return ''

# Iteração sobre os grupos
for _, group in grouped:
    accession = group['Protein.accession'].iloc[0]
    sequence_length = group['Sequence.length'].iloc[0]
    lengths = group['Sequence.length'].tolist()
    start_locations = group['Start.location'].tolist()
    stop_locations = group['Stop.location'].tolist()
    signature_descriptions = group['Signature.description'].tolist()

    color_combination = [get_color(desc) for desc in signature_descriptions]
    shape_combination = select_shape_combination(lengths)

    protein_domains = []
    for i in range(len(lengths)):
        domain_info = f"{shape_combination[i]}|{start_locations[i]}|{stop_locations[i]}|{color_combination[i]}"
        protein_domains.append(domain_info)
        # Atualizar a contagem de repetições do domínio
        domain = signature_descriptions[i]
        if domain in domain_count:
            domain_count[domain] += 1
        else:
            domain_count[domain] = 1

    protein_species = get_species(accession)

    protein_info.append((accession, sequence_length, protein_domains, protein_species))

# Ordenar as informações com base na ordem de Protein.accession
#protein_info.sort(key=lambda x: int(x[1]) if x[1] != 'nan' else 0)
#protein_info.sort(key=lambda x: int(float(x[1])) if str(x[1]) != 'nan' else float('inf'))
# Sort the protein_info list based on the order of Protein.accession in the input
protein_info.sort(key=lambda x: df.loc[df['Protein.accession'] == x[0]].index[0])

# Impressão das informações na tela
for info in protein_info:
    accession, sequence_length, domains, species = info
    domains_str = ','.join(domains)
    print(f"{accession},{sequence_length},{domains_str}")

# Escrita das informações no arquivo de saída
#with open(output_file, 'w') as file:
#    file.write('Protein.accession\tSequence.length\tshape\tStart.location\tStop.location\n')
#    for info in protein_info:
#        accession, sequence_length, domains, species = info
#        domains_str = ','.join(domains)
#        file.write(f'{accession},{sequence_length},{domains_str}\n')

# Escrita das informações no arquivo de saída
with open(output_file, 'w') as file:
    file.write('Protein.accession\tSequence.length\tshape\tStart.location\tStop.location\n')
    for info in protein_info:
        accession, sequence_length, domains, species = info
        domains_str = ','.join(domains)
        file.write(f'{accession},{sequence_length},{domains_str}\n')

# Criação do gráfico
