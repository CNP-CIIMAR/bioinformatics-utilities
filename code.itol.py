# Iteração sobre os grupos
for _, group in grouped:
    accession = group['Protein.accession'].iloc[0]
    lengths = group['Sequence.length'].tolist()
    start_locations = group['Start.location'].tolist()
    stop_locations = group['Stop.location'].tolist()

    color_combination = select_color_combination(len(lengths))
    shape_combination = select_shape_combination(lengths)

    protein_domains = []
    for i in range(len(lengths)):
        domain_info = f"{shape_combination[i]}|{start_locations[i]}|{stop_locations[i]}|{color_combination[i]}"
        protein_domains.append(domain_info)

    protein_info.append(f"{accession},{','.join(protein_domains)}")

# Ordenar as informações por comprimento da proteína
protein_info.sort(key=lambda x: int(x.split(',')[1]))

# Impressão das informações na tela
for info in protein_info:
    print(info)

# Escrita das informações no arquivo de saída
with open(output_file, 'w') as file:
    file.write('Protein.accession\tSequence.length\tshape\tStart.location\tStop.location\n')
for info in protein_info:
protein = import sys
import pandas as pd

# Tabela do InterProScan
table_file = sys.argv[1]
# Nome do arquivo de saída
output_file = sys.argv[2]

# Leitura da tabela como um DataFrame
df = pd.read_csv(table_file, sep='\t')

# Dicionário de cores disponíveis
colors = {
    0: '#ff0000',  # Red
    1: '#00ff00',  # Green
    2: '#0000ff',  # Blue
    3: '#ffff00',  # Yellow
    4: '#ff00ff',  # Magenta
    5: '#00ffff',  # Cyan
}

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

# Função para selecionar a melhor combinação de cores
def select_color_combination(n):
    selected_colors = []
    for i in range(n):
        selected_colors.append(colors[i % len(colors)])
    return selected_colors

# Função para selecionar a melhor combinação de shapes
def select_shape_combination(lengths):
    sorted_lengths = sorted(lengths)
    selected_shapes = []
    for length in lengths:
        index = sorted_lengths.index(length)
        selected_shapes.append(shapes[index % len(shapes)])
    return selected_shapes

# Agrupamento dos dados por Protein.accession
grouped = df.groupby('Protein.accession')

# Lista para armazenar as informações das proteínas
protein_info = []

# Iteração sobre os grupos
for _, group in grouped:
    accession = group['Protein.accession'].iloc[0]
    lengths = group['Sequence.length'].tolist()
    start_locations = group['Start.location'].tolist()
    stop_locations = group['Stop.location'].tolist()

    color_combination = select_color_combination(len(lengths))
    shape_combination = select_shape_combination(lengths)

    protein_domains = []
    for i in range(len(lengths)):
        domain_info = f"{shape_combination[i]}|{start_locations[i]}|{stop_locations[i]}|{color_combination[i]}"
        protein_domains.append(domain_info)

    protein_info.append(f"{accession},{','.join(protein_domains)}")

# Ordenar as informações por comprimento da proteína
protein_info.sort(key=lambda x: int(x.split(',')[1]))

# Impressão das informações na tela
for info in protein_info:
    print(info)

# Escrita das informações no arquivo de saída
with open(output_file, 'w') as file:
    file.write('Protein.accession\tSequence.length\tshape\tStart.location\tStop.location\n')
for info in protein_info:
protein = info.split(',')
accession = protein[0]
domains = protein[1:]
for domain in domains:
domain_info = domain.split('|')
shape = domain_info[0]
start_location = domain_info[1]
stop_location = domain_info[2]
color = domain_info[3]
file.write(f"{accession}\t{length}\t{shape}\t{start_location}\t{stop_location}\t{color}\n")
