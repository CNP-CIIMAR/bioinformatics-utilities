##########################
###Leandro de Mattos Pereira - CNP Team - Junior Resarcher.
### Code for get the taxonomy for the genomes annotated with Bigscape. 

#The code put the name of Specie in Organism collum and get the taxonomy. The fasta file need to be named acoording the specie or genus

import pandas as pd
import argparse
from Bio import Entrez
import time
import urllib

def fetch_taxonomy(genus_name):
    Entrez.email = ''  # Provide your email address for Entrez API

    # Search NCBI Taxonomy for the given genus name
    handle = Entrez.esearch(db='taxonomy', term=genus_name, retmode='xml')
    record = Entrez.read(handle)

    if record['IdList']:
        taxonomy_id = record['IdList'][0]
        
        # Intervalo de espera de 1 segundo antes de fazer a requisição
        time.sleep(1)

        # Fetch the taxonomy record for the given taxonomy ID
        handle = Entrez.efetch(db='taxonomy', id=taxonomy_id, retmode='xml')
        record = Entrez.read(handle)

        # Extract the taxonomic lineage
        lineage = record[0]['LineageEx']
        taxonomy = '; '.join([entry['ScientificName'] for entry in lineage])
        
        # Add the genus name after the taxonomic rank
        genus = genus_name.split()[0]
        taxonomy += '; ' + genus
        
        return taxonomy
    else:
        return None


def fill_organism(df):
    # Loop over the rows of the DataFrame
    for index, row in df.iterrows():
        # Check if the Organism cell is empty
        if pd.isnull(row['Organism']):
            # Get the value of the cluster identifier
            cluster_id = row['BGC']
            
            # Extract the name from the cluster identifier
            name = cluster_id.split('.cluster')[0].replace('_', ' ')

            # Fill the empty cell with the extracted name
            df.at[index, 'Organism'] = name

    return df

def fill_taxonomy(df):
    for index, row in df.iterrows():
        if pd.isnull(row['Taxonomy']):
            organism = row['Organism']
            
            if pd.isna(organism):
                continue
            
            genus = organism.split()[0]

            taxonomy = fetch_taxonomy(genus)

            if taxonomy:
                df.at[index, 'Taxonomy'] = taxonomy

    return df

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Preenche células vazias na coluna "Taxonomy" com informações obtidas do gênero presente na coluna "Organism" usando a biblioteca Biopython para consultar a base de dados do NCBI Taxonomy')
    parser.add_argument('input_file', help='Caminho do arquivo CSV de entrada')
    parser.add_argument('output_file', help='Caminho do arquivo CSV de saída')

    args = parser.parse_args()

    df = pd.read_csv(args.input_file, delimiter='\t')

    df_filled_organism = fill_organism(df)
    df_filled_taxonomy = fill_taxonomy(df_filled_organism)

    df_filled_taxonomy.to_csv(args.output_file, index=False)

