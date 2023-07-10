#############Autor Leandro de Mattos Pereira, Junior Researcher - CNP Lab.

from Bio import Entrez
import pandas as pd
import sys

def get_taxonomic_rank(protein_accession):
    Entrez.email = 'seu-email@aqui.com'
    handle = Entrez.efetch(db='protein', id=protein_accession, retmode='xml', idtype='acc')
    record = Entrez.read(handle)
    species = record[0]['GBSeq_organism']
    lineage = record[0]['GBSeq_taxonomy']
    handle.close()
    return species, lineage

def main(file_path, output_file_path):
    with open(file_path, 'r') as file:
        protein_accessions = file.read().splitlines()

    data = []

    for protein_accession in protein_accessions:
        try:
            species, lineage = get_taxonomic_rank(protein_accession)
            data.append([protein_accession, species, lineage])
        except Exception as e:
            print(f"Error retrieving data for Protein Accession: {protein_accession}")
            print(f"Protein Accession: {protein_accession} | Assembly Accession not found")
         #   print(str(e))
       # print()

    df = pd.DataFrame(data, columns=['Protein Accession', 'Species', 'Lineage'])
    df['Lineage'] = df['Lineage'].str.replace(';', '; ')

   # df.to_csv(output_file_path, index=False)
    df.to_csv(output_file_path, sep='\t', index=False)
    print(f"Results saved to {output_file_path}")

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python script.py protein_accessions.txt output.csv")
    else:
        file_path = sys.argv[1]
        output_file_path = sys.argv[2]
        main(file_path, output_file_path)

