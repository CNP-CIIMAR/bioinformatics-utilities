#Authors:
#Leandro de Mattos Pereira, Junior Researcher
# Pedro Leao - team Leader CNP 
import sys
from ete3 import NCBITaxa

def download_taxonomic_rank(species_list):
    ncbi = NCBITaxa()

    # Ler as espécies do arquivo
    with open(species_list, 'r') as file:
        species = [line.strip() for line in file if line.strip()]

    # Mapear os IDs taxonômicos das espécies
    taxids = ncbi.get_name_translator(species)

    # Obter os ranks taxonômicos para cada espécie
    ranks = {}
    for species_name, taxid_list in taxids.items():
        for taxid in taxid_list:
            lineage = ncbi.get_lineage(taxid)
            lineage_ranks = ncbi.get_rank(lineage)
            species_rank = {rank: ncbi.get_taxid_translator([taxid])[taxid] for (taxid, rank) in lineage_ranks.items()}
            ranks[species_name] = species_rank
            break  # Considerar apenas o primeiro ID taxonômico encontrado

    return ranks

# Verificar se o arquivo de lista de espécies foi fornecido como argumento
if len(sys.argv) < 3:
    print("Uso: python script.py arquivo_lista_especies.txt arquivo_saida.tsv")
    sys.exit(1)

# Argumento do arquivo de lista de espécies
species_list_file = sys.argv[1]
# Argumento do arquivo de saída
output_file = sys.argv[2]

# Baixar os ranks taxonômicos
taxonomic_ranks = download_taxonomic_rank(species_list_file)

# Obter uma lista de todos os ranks taxonômicos encontrados
all_ranks = set()
for ranks in taxonomic_ranks.values():
    all_ranks.update(ranks.keys())

# Ordenar os ranks em ordem alfabética
sorted_ranks = sorted(all_ranks)

# Escrever os resultados no arquivo de saída
with open(output_file, 'w') as file:
    # Escrever a linha de cabeçalho com os títulos das colunas
    file.write("Espécie\t")
    for rank in sorted_ranks:
        file.write(rank + "\t")
    file.write("\n")

    # Escrever os dados das espécies e ranks taxonômicos
    for species_name, ranks in taxonomic_ranks.items():
        file.write(species_name + "\t")
        for rank in sorted_ranks:
            if rank in ranks:
                file.write(ranks[rank] + "\t")
            else:
                file.write("\t")  # Deixar a célula vazia se o rank não estiver presente para a espécie
        file.write("\n")

print("Os ranks taxonômicos foram salvos no arquivo", output_file)
