#!/usr/bin/env python3
# organize_big_slice.py

import os
import pandas as pd
import shutil
import argparse
from collections import defaultdict
import re
import logging

def parse_arguments():
    """
    Parsa os argumentos da linha de comando.
    """
    parser = argparse.ArgumentParser(description="Organiza diretórios antiSMASH em datasets para BiG-SLiCE baseados no nível taxonômico Order.")
    parser.add_argument('--bigslice_dir', type=str, required=True,
                        help='Caminho para o diretório de entrada do BiG-SLiCE (onde os datasets serão criados).')
    parser.add_argument('--antismash_dir', type=str, required=True,
                        help='Caminho para o diretório de resultados do antiSMASH.')
    parser.add_argument('--taxonomy_table', type=str, required=True,
                        help='Caminho para a tabela de taxonomia (CSV).')
    parser.add_argument('--assembly_column', type=str, default='Assembly',
                        help='Nome da coluna na tabela de taxonomia que contém o Assembly Accession (default: "Assembly").')
    parser.add_argument('--lineage_column', type=str, default='Lineage',
                        help='Nome da coluna na tabela de taxonomia que contém a Lineage (default: "Lineage").')
    parser.add_argument('--order_index', type=int, required=True,
                        help='Índice (0-based) na lista de níveis da Lineage onde o Order está localizado.')
    parser.add_argument('--use_symlinks', action='store_true',
                        help='Use links simbólicos em vez de copiar arquivos.')
    parser.add_argument('--log_file', type=str, default='organize_big_slice.log',
                        help='Caminho para o arquivo de log (default: "organize_big_slice.log").')
    return parser.parse_args()

def setup_logging(log_file):
    """
    Configura o sistema de logging.
    """
    logging.basicConfig(
        filename=log_file,
        filemode='w',
        level=logging.INFO,
        format='%(asctime)s %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    # Também loga no console
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(levelname)s: %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)

def extract_assembly_accession(dir_name):
    """
    Extrai o Assembly Accession do nome do diretório usando regex.
    Exemplo: 'Result_GCA_000007865.1.fasta_genomic.fna' -> 'GCA_000007865.1'
             'Result_GCA_910592015.1_Diversispora_eburnea_AZ414A_genomic.fna' -> 'GCA_910592015.1'
             'Result_GCA_027600645.1_ASM2760064v1_genomic.fna' -> 'GCA_027600645.1'
    """
    pattern = r'Result_(GCA|GCF)_\d+\.\d+'  # Captura GCA_123456789.1 ou GCF_123456789.1
    match = re.search(pattern, dir_name)
    if match:
        return match.group(0).split('_')[1] + '_' + match.group(0).split('_')[2]
    else:
        return None

def get_order(lineage, order_index):
    """
    Extrai o nível taxonômico 'Order' a partir da coluna 'Lineage'.
    :param lineage: String contendo a Lineage separada por ';'
    :param order_index: Índice (0-based) onde 'Order' está localizado na lista de níveis
    :return: String com o nome do 'Order' ou 'Unknown' se não encontrado
    """
    levels = [lvl.strip() for lvl in lineage.split(';')]
    try:
        if order_index < len(levels):
            order = levels[order_index]
            return order if order else 'Unknown'
        else:
            return 'Unknown'
    except IndexError:
        return 'Unknown'

def main():
    args = parse_arguments()
    setup_logging(args.log_file)
    
    bigslice_dir = args.bigslice_dir
    antismash_dir = args.antismash_dir
    taxonomy_table_path = args.taxonomy_table
    assembly_column = args.assembly_column
    lineage_column = args.lineage_column
    order_index = args.order_index
    use_symlinks = args.use_symlinks
    
    logging.info("Iniciando o processo de organização dos datasets para BiG-SLiCE.")
    logging.info(f"Diretório BiG-SLiCE: {bigslice_dir}")
    logging.info(f"Diretório antiSMASH: {antismash_dir}")
    logging.info(f"Tabela de Taxonomia: {taxonomy_table_path}")
    logging.info(f"Coluna de Assembly: {assembly_column}")
    logging.info(f"Coluna de Lineage: {lineage_column}")
    logging.info(f"Índice do Order na Lineage: {order_index}")
    logging.info(f"Usar links simbólicos: {'Sim' if use_symlinks else 'Não'}")
    
    # Verificar se os caminhos fornecidos existem
    if not os.path.isdir(antismash_dir):
        logging.error(f"O diretório antiSMASH '{antismash_dir}' não existe.")
        return
    if not os.path.isfile(taxonomy_table_path):
        logging.error(f"A tabela de taxonomia '{taxonomy_table_path}' não existe.")
        return
    os.makedirs(bigslice_dir, exist_ok=True)
    
    # Carregar a tabela de taxonomia
    logging.info("Carregando a tabela de taxonomia...")
    try:
        df_tax = pd.read_csv(taxonomy_table_path)
    except Exception as e:
        logging.error(f"Erro ao ler a tabela de taxonomia: {e}")
        return
    
    # Verificar colunas necessárias
    required_columns = [assembly_column, lineage_column]
    for col in required_columns:
        if col not in df_tax.columns:
            logging.error(f"A coluna '{col}' é necessária na tabela de taxonomia.")
            return
    
    # Extrair o nível taxonômico 'Order' para cada Assembly Accession
    logging.info("Extraindo o nível taxonômico 'Order'...")
    df_tax['Order'] = df_tax[lineage_column].apply(lambda x: get_order(x, order_index))
    
    # Criar um mapeamento de Assembly Accession para Order
    accession_to_order = pd.Series(df_tax.Order.values, index=df_tax[assembly_column]).to_dict()
    
    # Obter a lista única de Orders
    unique_orders = sorted(df_tax['Order'].unique())
    
    # Mapear cada Order para um número de dataset
    order_to_dataset = {order: f"dataset_{i+1}" for i, order in enumerate(unique_orders)}
    
    logging.info(f"Total de Orders únicos: {len(unique_orders)}")
    
    # Inicializar um dicionário para armazenar informações de taxonomia por dataset
    taxonomy_data = defaultdict(list)
    
    # Inicializar o arquivo datasets.tsv
    datasets_tsv_path = os.path.join(bigslice_dir, "datasets.tsv")
    try:
        with open(datasets_tsv_path, 'w') as f:
            f.write("#Dataset name\tPath to dataset folder\tPath to taxonomy file\tDescription of the dataset\n")
    except Exception as e:
        logging.error(f"Erro ao criar o arquivo 'datasets.tsv': {e}")
        return
    
    # Criar a pasta taxonomy dentro do bigslice_dir
    taxonomy_folder = os.path.join(bigslice_dir, "taxonomy")
    os.makedirs(taxonomy_folder, exist_ok=True)
    
    # Listar todos os diretórios no antismash_dir
    all_result_dirs = [d for d in os.listdir(antismash_dir) if os.path.isdir(os.path.join(antismash_dir, d))]
    logging.info(f"Total de diretórios de resultado antiSMASH encontrados: {len(all_result_dirs)}")
    
    # Iterar sobre os diretórios do antiSMASH
    for idx, result_dir in enumerate(all_result_dirs, 1):
        if idx % 1000 == 0:
            logging.info(f"Processando {idx}/{len(all_result_dirs)} diretórios...")
        
        accession = extract_assembly_accession(result_dir)
        if not accession:
            logging.warning(f"O nome do diretório '{result_dir}' não corresponde ao padrão esperado. Pulando...")
            continue
        
        # Obter o Order correspondente
        order = accession_to_order.get(accession, 'Unknown')
        
        # Obter o nome do dataset correspondente ao Order
        dataset_name = order_to_dataset.get(order, 'dataset_unknown')
        
        # Caminho completo do diretório de resultado
        result_dir_path = os.path.join(antismash_dir, result_dir)
        
        # Caminho do dataset de destino
        dataset_dir_path = os.path.join(bigslice_dir, dataset_name)
        os.makedirs(dataset_dir_path, exist_ok=True)
        
        # Nome da pasta do genome no dataset
        genome_folder_name = f"genome_{accession}"
        genome_dir_path = os.path.join(dataset_dir_path, genome_folder_name)
        os.makedirs(genome_dir_path, exist_ok=True)
        
        # Copiar os arquivos .gbk para o diretório do genome no dataset
        gbk_files_copied = False
        for file in os.listdir(result_dir_path):
            if file.endswith(".gbk"):
                src_file = os.path.join(result_dir_path, file)
                dst_file = os.path.join(genome_dir_path, file)
                try:
                    if use_symlinks:
                        # Criar link simbólico (verifica se já existe)
                        if not os.path.exists(dst_file):
                            os.symlink(src_file, dst_file)
                            logging.info(f"Link simbólico criado: {dst_file} -> {src_file}")
                    else:
                        # Copiar o arquivo
                        shutil.copy2(src_file, dst_file)
                        logging.info(f"Arquivo copiado: {src_file} -> {dst_file}")
                    gbk_files_copied = True
                except Exception as e:
                    logging.error(f"Erro ao copiar '{src_file}' para '{dst_file}': {e}")
        
        if not gbk_files_copied:
            logging.warning(f"Nenhum arquivo '.gbk' encontrado no diretório '{result_dir_path}'.")
        
        # Adicionar informações de taxonomia ao dicionário
        tax_info = df_tax[df_tax[assembly_column] == accession]
        if not tax_info.empty:
            tax_info = tax_info.iloc[0]
            lineage_levels = [lvl.strip() for lvl in tax_info[lineage_column].split(';')]
            genome_folder_entry = {
                'Genome folder name': f"{genome_folder_name}/",
                'Kingdom/Domain': lineage_levels[2] if len(lineage_levels) > 2 else 'Unknown',
                'Class': lineage_levels[4] if len(lineage_levels) > 4 else 'Unknown',
                'Order': order,
                'Family': lineage_levels[6] if len(lineage_levels) > 6 else 'Unknown',
                'Genus': lineage_levels[7] if len(lineage_levels) > 7 else 'Unknown',
                'Species': lineage_levels[8] if len(lineage_levels) > 8 else 'Unknown',
                'Organism/Strain': lineage_levels[9] if len(lineage_levels) > 9 else 'Unknown'
            }
            taxonomy_data[dataset_name].append(genome_folder_entry)
        else:
            logging.warning(f"Informações de taxonomia para '{accession}' não encontradas.")
    
    # === Gerar os Arquivos taxonomy_dataset_X.tsv ===
    
    logging.info("Gerando arquivos 'taxonomy_dataset_X.tsv'...")
    for dataset, entries in taxonomy_data.items():
        taxonomy_file_path = os.path.join(taxonomy_folder, f"taxonomy_{dataset}.tsv")
        df_taxonomy = pd.DataFrame(entries)
        try:
            df_taxonomy.to_csv(taxonomy_file_path, sep='\t', index=False)
            logging.info(f"Arquivo de taxonomia criado: {taxonomy_file_path}")
        except Exception as e:
            logging.error(f"Erro ao criar o arquivo de taxonomia '{taxonomy_file_path}': {e}")
            continue
        
        # Adicionar entrada no datasets.tsv
        with open(datasets_tsv_path, 'a') as f:
            description = f"Dataset agrupado por Order: {dataset.replace('dataset_', '')}"
            f.write(f"{dataset}\t{dataset}\ttaxonomy/taxonomy_{dataset}.tsv\t{description}\n")
    
    logging.info("Gerando o arquivo 'datasets.tsv'...")
    
    logging.info("Organização de datasets concluída com sucesso!")
    logging.info(f"Arquivo 'datasets.tsv' criado em: {datasets_tsv_path}")

if __name__ == "__main__":
    main()
