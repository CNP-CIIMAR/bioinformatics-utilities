import pandas as pd
import argparse
## Authors: Leandro de Mattos Pereira
## Date: 22 February 2024
## CNP -team - Dr. Pedro Leao, Team Leader.


def process_table(file_path):
    # Lendo os dados do arquivo CSV
    df = pd.read_csv(file_path, sep='\t')
    
    # Criando um OTU ID Ãºnico para cada Protein Accession
    df['OTU ID'] = range(1, len(df) + 1)
    
    # Dividindo a coluna Lineage em mÃºltiplas colunas com base no separador ';'
    lineage_cols = df['Lineage'].str.split(';', expand=True)
    
    # Nomeando as colunas de acordo com o ranking taxonÃ´mico
    taxonomic_ranks = ['Superkingdom', 'Phylum', 'Class', 'Order', 'Suborder', 'Family', 'Genus', 'Species']
    lineage_cols.columns = taxonomic_ranks[:lineage_cols.shape[1]]
    
    # Concatenando as novas colunas ao DataFrame original
    df = pd.concat([df, lineage_cols], axis=1)
    
    # Removendo a coluna Lineage original
    df.drop('Lineage', axis=1, inplace=True)
    
    return df

def aggregate_and_save(df, level, with_protein_accession, taxonomic_ranks):
    # Agregar dados baseados no nÃ­vel taxonÃ´mico atual
    cols = taxonomic_ranks[:level + 1]
    if level == len(taxonomic_ranks) - 1:  # Para o nÃ­vel de Species
        df_grouped = df.groupby(cols).size().reset_index(name='Counts')
        if with_protein_accession:  # Se necessÃ¡rio incluir o Protein Accession
            df_grouped_with_pa = pd.merge(df_grouped, df[['Species', 'Protein Accession']], on='Species')
            df_grouped_with_pa['OTU_ID'] = (df_grouped_with_pa.index + 1).astype(str).str.zfill(4)
            arquivo_saida_with_pa = f'tabela_{level}_biom_with_protein_accession.tsv'
            df_grouped_with_pa.to_csv(arquivo_saida_with_pa, sep='\t', index=False)
            print(f'Arquivo {arquivo_saida_with_pa} gerado com sucesso.')
    else:
        df_grouped = df.groupby(cols).size().reset_index(name='Counts')
        df_grouped['OTU_ID'] = (df_grouped.index + 1).astype(str).str.zfill(4)  # Adicionando OTU_ID
        df_grouped = df_grouped[['OTU_ID'] + cols + ['Counts']]  # Reordenando as colunas para OTU_ID aparecer primeiro
        arquivo_saida = f'tabela_{level}_biom_without_protein_accession.tsv'
        df_grouped.to_csv(arquivo_saida, sep='\t', index=False)
        print(f'Arquivo {arquivo_saida} gerado com sucesso.')


if __name__ == "__main__":
    # Configurando o argparse para aceitar um caminho de arquivo via linha de comando
    parser = argparse.ArgumentParser(description="Processa uma tabela e divide a coluna Lineage em mÃºltiplos nÃ­veis taxonÃ´micos.")
    parser.add_argument('file_path', type=str, help='O caminho do arquivo CSV contendo a tabela para processar.')
    
    args = parser.parse_args()
    
    # Processando a tabela
    df_processed = process_table(args.file_path)
    
    # Definindo os nÃ­veis taxonÃ´micos
    taxonomic_ranks = ['Superkingdom', 'Phylum', 'Class', 'Order', 'Suborder', 'Family', 'Genus', 'Species']
    
    # Agregar e salvar para todos os nÃ­veis taxonÃ´micos
    for level in range(len(taxonomic_ranks)):
        aggregate_and_save(df_processed, level, False, taxonomic_ranks)  # Sem Protein Accession
        if level == len(taxonomic_ranks) - 1:  # Para o nÃ­vel de Species
            aggregate_and_save(df_processed, level, True, taxonomic_ranks)  # Com Protein Accession
