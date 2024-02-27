import pandas as pd
import argparse
import pandas as pd
import argparse
## Authors: Leandro de Mattos Pereira
## Date: 22 February 2024
## CNP -team - Dr. Pedro Leao, Team Leader.

def process_table(file_path):
    # Reading data from the CSV file
    df = pd.read_csv(file_path, sep='\t')
    
    # Creating a unique OTU ID for each Protein Accession
    df['OTU ID'] = range(1, len(df) + 1)
    
    # Splitting the Lineage column into multiple columns based on the separator ';'
    lineage_cols = df['Lineage'].str.split(';', expand=True)
    
    # Naming the columns according to the taxonomic ranking
    taxonomic_ranks = ['Superkingdom', 'Phylum', 'Class', 'Order', 'Suborder', 'Family', 'Genus', 'Species']
    lineage_cols.columns = taxonomic_ranks[:lineage_cols.shape[1]]
    
    # Concatenating the new columns to the original DataFrame
    df = pd.concat([df, lineage_cols], axis=1)
    
    # Removing the original Lineage column
    df.drop('Lineage', axis=1, inplace=True)
    
    return df

def aggregate_and_save(df, level, with_protein_accession, taxonomic_ranks):
    # Aggregate data based on the current taxonomic level
    cols = taxonomic_ranks[:level + 1]
    if level == len(taxonomic_ranks) - 1:  # For the Species level
        df_grouped = df.groupby(cols).size().reset_index(name='Counts')
        if with_protein_accession:  # If Protein Accession is to be included
            df_grouped_with_pa = pd.merge(df_grouped, df[['Species', 'Protein Accession', 'OTU ID']], on='Species')
            df_grouped_with_pa.insert(0, 'OTU ID', df_grouped_with_pa.pop('OTU ID'))  # Move OTU ID to the first column
            arquivo_saida_with_pa = f'tabela_{level}_biom_with_protein_accession.tsv'
            df_grouped_with_pa.columns = map(str.lower, df_grouped_with_pa.columns)  # Change headers to lowercase
            df_grouped_with_pa.to_csv(arquivo_saida_with_pa, sep='\t', index=False)
            print(f'Arquivo {arquivo_saida_with_pa} gerado com sucesso.')
        else:  # Without Protein Accession
            df_grouped.insert(0, 'OTU ID', (df_grouped.index + 1).astype(str).str.zfill(4))  # Add OTU ID as first column
            arquivo_saida = f'tabela_{level}_biom_without_protein_accession.tsv'
            df_grouped.columns = map(str.lower, df_grouped.columns)  # Change headers to lowercase
            df_grouped.to_csv(arquivo_saida, sep='\t', index=False)
            print(f'Arquivo {arquivo_saida} gerado com sucesso.')
    else:
        df_grouped = df.groupby(cols).size().reset_index(name='Counts')
        df_grouped.insert(0, 'OTU ID', (df_grouped.index + 1).astype(str).str.zfill(4))  # Add OTU ID as first column
        arquivo_saida = f'tabela_{level}_biom_without_protein_accession.tsv'
        df_grouped.columns = map(str.lower, df_grouped.columns)  # Change headers to lowercase
        df_grouped.to_csv(arquivo_saida, sep='\t', index=False)
        print(f'Arquivo {arquivo_saida} gerado com sucesso.')


if __name__ == "__main__":
    # Setting up argparse to accept a file path via command line
    parser = argparse.ArgumentParser(description="Processa uma tabela e divide a coluna Lineage em mÃºltiplos nÃ­veis taxonÃ´micos.")
    parser.add_argument('file_path', type=str, help='O caminho do arquivo CSV contendo a tabela para processar.')
    
    args = parser.parse_args()
    
    # Processing the table
    df_processed = process_table(args.file_path)
    

    # Defining taxonomic ranks
    taxonomic_ranks = ['Superkingdom', 'Phylum', 'Class', 'Order', 'Suborder', 'Family', 'Genus', 'Species']
    
    # Aggregating and saving for all taxonomic levels
    for level in range(len(taxonomic_ranks)):
        aggregate_and_save(df_processed, level, False, taxonomic_ranks)  # Without Protein Accession
        if level == len(taxonomic_ranks) - 1:  # For the Species level
            aggregate_and_save(df_processed, level, True, taxonomic_ranks)  # With Protein Accession

