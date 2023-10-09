import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

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
# Dicionário para armazenar a contagem de repetições de cada domínio em cada proteína
domain_count = defaultdict(lambda: defaultdict(int))

# Leitura da tabela de espécies
species_df = pd.read_csv(species_file, sep='\t')

# Função para obter o nome da espécie para uma dada Protein.accession
def get_species(protein_accession):
    species = species_df.loc[species_df['Protein.accession'] == protein_accession, 'species']
    if len(species) > 0:
        return species.values[0]
    return ''

# Função para obter o nome do filo para uma dada Protein.accession
def get_phylum(protein_accession):
    phylum = species_df.loc[species_df['Protein.accession'] == protein_accession, 'phylum']
    if len(phylum) > 0:
        return phylum.values[0]
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
        domain_count[accession][domain] += 1

    protein_species = get_species(accession)
    protein_phylum = get_phylum(accession)

    protein_info.append((accession, sequence_length, protein_domains, protein_species, protein_phylum))

# Ordenar as informações com base na ordem de Protein.accession
protein_info.sort(key=lambda x: int(float(x[1])) if str(x[1]) != 'nan' else float('inf'))

# Impressão das informações na tela
for info in protein_info:
    accession, sequence_length, domains, species, phylum = info
    domains_str = ','.join(domains)
    print(f"{accession},{sequence_length},{domains_str}")

# Escrita das informações no arquivo de saída
with open(output_file, 'w') as file:
    file.write('Protein.accession\tSequence.length\tshape\tStart.location\tStop.location\n')
    for info in protein_info:
        accession, sequence_length, domains, species, phylum = info
        domains_str = ','.join(domains)
        file.write(f'{accession},{sequence_length},{domains_str}\n')

# Criação da tabela de contagem de domínios por proteína e por espécie
count_table = pd.DataFrame.from_dict(domain_count).fillna(0)

# Impressão da tabela no browser
count_table_html = count_table.to_html()
with open("domain_count.html", "w") as file:
    file.write(count_table_html)

# Impressão da tabela no arquivo de texto
count_table.to_csv("domain_count.txt", sep='\t', index=True)
# Bar plot of domain counts
domain_counts = count_table.sum(axis=1)
sorted_counts = domain_counts.sort_values(ascending=False)
sorted_counts.plot(kind='bar')
plt.xlabel('Protein Accession')
plt.ylabel('Domain Count')
plt.title('Domain Counts per Protein')

# Save the figure in high resolution (300 dpi)
plt.savefig('domain_counts_300dpi.png', dpi=300)

# Save the figure in high resolution (500 dpi)
plt.savefig('domain_counts_500dpi.png', dpi=500)

# Display the plot on the screen
plt.show()

# Pie chart of species distribution
species_counts = pd.Series([info[3] for info in protein_info]).value_counts()

# Select top 40 species
top_species_counts = species_counts[:40]
other_species_count = species_counts[40:].sum()

# Combine other species into a single entry
top_species_counts['Others'] = other_species_count

# Plot the pie chart
plt.pie(top_species_counts, labels=top_species_counts.index, autopct='%1.1f%%')
plt.axis('equal')
plt.title('Protein Distribution by Species')

# Save the figure in high resolution (300 dpi)
plt.savefig('species_distribution_300dpi.png', dpi=300)

# Save the figure in high resolution (500 dpi)
plt.savefig('species_distribution_500dpi.png', dpi=500)

# Display the pie chart on the screen
plt.show()

# Pie chart of phylum distribution
phylum_counts = pd.Series([info[4] for info in protein_info]).value_counts()

# Select top 24 phyla
top_phylum_counts = phylum_counts[:24]
other_phylum_count = phylum_counts[24:].sum()

# Combine other phyla into a single entry
top_phylum_counts['Others'] = other_phylum_count

# Plot the pie chart
plt.pie(top_phylum_counts, labels=top_phylum_counts.index, autopct='%1.1f%%')
plt.axis('equal')
plt.title('Protein Distribution by Phylum')

# Save the figure in high resolution (300 dpi)
plt.savefig('phylum_distribution_300dpi.png', dpi=300)

# Save the figure in high resolution (500 dpi)
plt.savefig('phylum_distribution_500dpi.png', dpi=500)

# Display the pie chart on the screen
plt.show()

############################################################# Parte II
# Função para gerar o desenho esquemático das proteínas com múltiplas assinaturas FAAL
def generate_protein_diagram(protein_info, output_file):
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    ax.set_xlim(0, 100)
    ax.set_ylim(0, len(protein_info) * 10)
    ax.set_axis_off()

    for i, info in enumerate(protein_info):
        accession, sequence_length, domains, species, phylum = info

        # Only consider proteins with more than one FAAL domain
        faal_domains = [domain for domain in domains if 'FAAL' in domain]

        if len(faal_domains) > 1:
            y = i * 10

            for domain in faal_domains:
                shape, start, stop, color = domain.split('|')
                shape_width = float(stop) - float(start)
                shape_height = sequence_length / 1000.0

                rectangle = plt.Rectangle((float(start), y), shape_width, shape_height, facecolor=color, edgecolor='black')
                ax.add_patch(rectangle)

                plt.text(float(start) + shape_width / 2, y + shape_height / 2, shape, ha='center', va='center')

                # Update the domain count
                domain_name = domain.split('|')[3]
                domain_count[accession][domain_name] += 1

    plt.savefig(output_file, dpi=300)
    plt.show()
