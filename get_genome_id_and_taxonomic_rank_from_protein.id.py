import pandas as pd
from Bio import Entrez
import subprocess
import sys
## Authors: Leandro de Mattos Pereira
## Date: 22 February 2024
## CNP -team - Dr. Pedro Leao, Team Leader.

# Verifique se os argumentos necessÃ¡rios foram fornecidos
if len(sys.argv) < 4:
    print("Usage: python script.py <input_filename> <genome_output_filename> <taxonomic_output_filename>")
    sys.exit(1)

input_filename = sys.argv[1]
genome_output_filename = sys.argv[2]
taxonomic_output_filename = sys.argv[3]

# Define your email here for Entrez
Entrez.email = 'your_email@example.com'

def get_taxonomic_rank(protein_accession):
    try:
        handle = Entrez.efetch(db='protein', id=protein_accession, retmode='xml')
        records = Entrez.read(handle)
        species = records[0]['GBSeq_organism']
        lineage = records[0]['GBSeq_taxonomy']
        handle.close()
        return species, lineage
    except Exception as e:
        print(f"Error retrieving data for Protein Accession: {protein_accession}")
        print(str(e))
        return None, None

# Prepara a tabela de genomas
with open(input_filename, 'r') as input_file, open(genome_output_filename, 'w') as genome_file:
    for line in input_file:
        protein_accession = line.strip()
        if protein_accession:
            # Chamada ao subprocesso para efetch
            try:
                ipg_result = subprocess.check_output(['efetch', '-db', 'protein', '-id', protein_accession, '-format', 'ipg'], universal_newlines=True).strip()
                refseq_info = None
                insdc_info = None
                for line in ipg_result.split('\n'):
                    if 'RefSeq' in line:
                        refseq_info = line.split()[-1]  # Presume-se que o ID desejado esteja na Ãºltima coluna
                    elif 'INSDC' in line:
                        insdc_info = line.split()[-1]
                genome_accession = refseq_info if refseq_info else insdc_info
                if genome_accession:
                    genome_file.write(f"{protein_accession}\t{genome_accession}\n")
            except subprocess.CalledProcessError as e:
                print(f"Failed to fetch data for {protein_accession}: {e}")

# Prepara a tabela taxonÃ´mica usando a saÃ­da do passo anterior
data = []
with open(genome_output_filename, 'r') as genome_file:
    for line in genome_file:
        protein_accession, genome_accession = line.strip().split('\t')
        species, lineage = get_taxonomic_rank(protein_accession)
        if species and lineage:
            data.append([protein_accession, genome_accession, species, lineage])

df_taxonomic = pd.DataFrame(data, columns=['Protein Accession', 'Genome Accession', 'Species', 'Lineage'])
df_taxonomic.to_csv(taxonomic_output_filename, sep='\t', index=False)

print(f"Genome IDs saved to {genome_output_filename}")
print(f"Taxonomic information saved to {taxonomic_output_filename}")
