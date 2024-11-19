#!/usr/bin/env python3
# organize_big_slice.py

# Author: Leandro de Mattos Pereira -november 2024

import os
import sys
import pandas as pd
import shutil
import argparse
from collections import defaultdict
import re
import logging
from logging.handlers import RotatingFileHandler

def parse_arguments():
    """
    Analisa os argumentos da linha de comando.
    """
    parser = argparse.ArgumentParser(description="Organiza diretórios antiSMASH em datasets para BiG-SLiCE baseados no nível taxonômico Order.")
    parser.add_argument('--bigslice_dir', type=str, required=True,
                        help='Caminho para o diretório de entrada do BiG-SLiCE (onde os datasets serão criados).')
    parser.add_argument('--antismash_dir', type=str, required=True,
                        help='Caminho para o diretório de resultados do antiSMASH.')
    parser.add_argument('--taxonomy_table', type=str, required=True,
                        help='Caminho para a tabela de taxonomia (TSV).')
    parser.add_argument('--assembly_column', type=str, default='Assembly Accession',
                        help='Nome da coluna na tabela de taxonomia que contém o Assembly Accession (default: "Assembly Accession").')
    parser.add_argument('--lineage_column', type=str, default='Lineage',
                        help='Nome da coluna na tabela de taxonomia que contém a Lineage (default: "Lineage").')
    parser.add_argument('--log_file', type=str, default='organize_big_slice.log',
                        help='Caminho para o arquivo de log (default: "organize_big_slice.log").')
    parser.add_argument('--dry_run', action='store_true',
                        help='Realiza uma execução de teste sem fazer alterações.')
    return parser.parse_args()

def setup_logging(log_file):
    """
    Configura o sistema de logging com rotação.
    """
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    # Handler para arquivo com rotação
    file_handler = RotatingFileHandler(log_file, maxBytes=10**6, backupCount=5)
    file_handler.setLevel(logging.INFO)
    file_formatter = logging.Formatter('%(asctime)s %(levelname)s: %(message)s', '%Y-%m-%d %H:%M:%S')
    file_handler.setFormatter(file_formatter)
    logger.addHandler(file_handler)

    # Handler para console
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    console_formatter = logging.Formatter('%(levelname)s: %(message)s')
    console.setFormatter(console_formatter)
    logger.addHandler(console)

def extract_assembly_accession(dir_name):
    """
    Extrai o Assembly Accession do nome do diretório usando regex.
    Exemplos:
        'Result_GCA_000007865.1.fasta_genomic.fna' -> 'GCA_000007865.1'
        'Result_GCA_910592015.1_Diversispora_eburnea_AZ414A_genomic.fna' -> 'GCA_910592015.1'
        'Result_GCA_027600645.1_ASM2760064v1_genomic.fna' -> 'GCA_027600645.1'
    """
    pattern = r'Result_(GCA|GCF)_(\d+\.\d+)'  # Captura GCA_123456789.1 ou GCF_123456789.1
    match = re.search(pattern, dir_name, re.IGNORECASE)
    if match:
        return f"{match.group(1)}_{match.group(2)}"
    else:
        return None

def get_order(lineage):
    """
    Extrai o nível taxonômico 'Order' a partir da coluna 'Lineage'.
    Procura pelo nível que termina com 'ales'.
    
    :param lineage: String contendo a Lineage separada por ';'
    :return: String com o nome do 'Order' ou 'Unknown' se não encontrado
    """
    levels = [lvl.strip() for lvl in lineage.split(';')]
    for level in levels:
        if level.endswith('ales'):
            return level
    return 'Unknown'

def load_taxonomy(taxonomy_table_path, assembly_column, lineage_column):
    """
    Carrega e processa a tabela de taxonomia.
    
    :return: DataFrame com uma coluna adicional 'Order'
    """
    logging.info("Carregando a tabela de taxonomia...")
    try:
        df_tax = pd.read_csv(taxonomy_table_path, sep='\t')
    except FileNotFoundError:
        logging.error(f"A tabela de taxonomia '{taxonomy_table_path}' não foi encontrada.")
        sys.exit(1)
    except pd.errors.ParserError as e:
        logging.error(f"Erro ao analisar a tabela de taxonomia: {e}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Erro inesperado ao ler a tabela de taxonomia: {e}")
        sys.exit(1)
    
    # Verificar nomes das colunas
    if assembly_column not in df_tax.columns or lineage_column not in df_tax.columns:
        logging.error(f"As colunas '{assembly_column}' e/ou '{lineage_column}' não foram encontradas na tabela de taxonomia.")
        logging.error(f"Colunas disponíveis: {df_tax.columns.tolist()}")
        sys.exit(1)
    
    # Extrair o nível 'Order'
    logging.info("Extraindo os níveis taxonômicos 'Order'...")
    df_tax['Order'] = df_tax[lineage_column].apply(get_order)
    return df_tax

def create_dataset_mapping(df_tax, assembly_column):
    """
    Cria um mapeamento de Order para nomes de datasets.
    
    :return: Tupla (dicionário order_to_dataset, dicionário accession_to_order)
    """
    unique_orders = sorted(df_tax['Order'].unique())
    order_to_dataset = {order: f"dataset_{i+1}" for i, order in enumerate(unique_orders)}
    accession_to_order = pd.Series(df_tax.Order.values, index=df_tax[assembly_column]).to_dict()
    
    logging.info(f"Total de Orders únicos: {len(unique_orders)}")
    return order_to_dataset, accession_to_order

def initialize_datasets_tsv(bigslice_dir):
    """
    Inicializa o arquivo datasets.tsv com os cabeçalhos.
    
    :return: Caminho para o arquivo datasets.tsv
    """
    datasets_tsv_path = os.path.join(bigslice_dir, "datasets.tsv")
    try:
        with open(datasets_tsv_path, 'w') as f:
            f.write("#Dataset name\tPath to dataset folder\tPath to taxonomy file\tDescription of the dataset\n")
    except Exception as e:
        logging.error(f"Erro ao criar o arquivo 'datasets.tsv': {e}")
        sys.exit(1)
    return datasets_tsv_path

def process_antismash_directories(antismash_dir, bigslice_dir, dry_run, 
                                 order_to_dataset, accession_to_order, df_tax, assembly_column, lineage_column):
    """
    Processa os diretórios de resultados do antiSMASH e organiza-os em datasets do BiG-SLiCE.
    
    :return: Dicionário taxonomy_data mapeando nomes de datasets para entradas de taxonomia
    """
    taxonomy_data = defaultdict(list)
    taxonomy_folder = os.path.join(bigslice_dir, "taxonomy")
    if not dry_run:
        os.makedirs(taxonomy_folder, exist_ok=True)
    
    all_result_dirs = [d for d in os.listdir(antismash_dir) if os.path.isdir(os.path.join(antismash_dir, d))]
    logging.info(f"Total de diretórios de resultado antiSMASH encontrados: {len(all_result_dirs)}")
    
    for idx, result_dir in enumerate(all_result_dirs, 1):
        if idx % 1000 == 0:
            logging.info(f"Processando {idx}/{len(all_result_dirs)} diretórios...")
        
        accession = extract_assembly_accession(result_dir)
        if not accession:
            logging.warning(f"O nome do diretório '{result_dir}' não corresponde ao padrão esperado. Pulando...")
            continue
        
        order = accession_to_order.get(accession, 'Unknown')
        dataset_name = order_to_dataset.get(order, 'dataset_unknown')
        
        result_dir_path = os.path.join(antismash_dir, result_dir)
        dataset_dir_path = os.path.join(bigslice_dir, dataset_name)
        genome_folder_name = f"genome_{accession}"
        genome_dir_path = os.path.join(dataset_dir_path, genome_folder_name)
        
        if not dry_run:
            os.makedirs(genome_dir_path, exist_ok=True)
        
        gbk_files_copied = False
        for file in os.listdir(result_dir_path):
            if file.endswith(".gbk"):
                src_file = os.path.join(result_dir_path, file)
                dst_file = os.path.join(genome_dir_path, file)
                try:
                    if not dry_run:
                        shutil.copy2(src_file, dst_file)
                    logging.info(f"Arquivo copiado: {src_file} -> {dst_file}")
                    gbk_files_copied = True
                except Exception as e:
                    logging.error(f"Erro ao copiar o arquivo '{src_file}' para '{dst_file}': {e}")
        
        if not gbk_files_copied:
            logging.warning(f"Nenhum arquivo '.gbk' encontrado no diretório '{result_dir_path}'.")
        
        # Recuperar informações de taxonomia
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
    
    return taxonomy_data

def generate_taxonomy_files(taxonomy_data, taxonomy_folder, bigslice_dir, dry_run):
    """
    Gera arquivos TSV de taxonomia para cada dataset e atualiza o arquivo datasets.tsv.
    """
    datasets_tsv_path = os.path.join(bigslice_dir, "datasets.tsv")
    for dataset, entries in taxonomy_data.items():
        taxonomy_file_path = os.path.join(taxonomy_folder, f"taxonomy_{dataset}.tsv")
        df_taxonomy = pd.DataFrame(entries)
        if not dry_run:
            try:
                df_taxonomy.to_csv(taxonomy_file_path, sep='\t', index=False)
                logging.info(f"Arquivo de taxonomia criado: {taxonomy_file_path}")
            except Exception as e:
                logging.error(f"Erro ao criar o arquivo de taxonomia '{taxonomy_file_path}': {e}")
                continue
        
        # Atualizar datasets.tsv
        try:
            with open(datasets_tsv_path, 'a') as f:
                description = f"Dataset agrupado por Order: {dataset.replace('dataset_', '')}"
                f.write(f"{dataset}\t{dataset}\ttaxonomy/taxonomy_{dataset}.tsv\t{description}\n")
        except Exception as e:
            logging.error(f"Erro ao atualizar 'datasets.tsv' com o dataset '{dataset}': {e}")

def main():
    args = parse_arguments()
    setup_logging(args.log_file)
    
    bigslice_dir = args.bigslice_dir
    antismash_dir = args.antismash_dir
    taxonomy_table_path = args.taxonomy_table
    assembly_column = args.assembly_column
    lineage_column = args.lineage_column
    dry_run = args.dry_run
    
    logging.info("Iniciando o processo de organização dos datasets para BiG-SLiCE.")
    logging.info(f"Diretório BiG-SLiCE: {bigslice_dir}")
    logging.info(f"Diretório antiSMASH: {antismash_dir}")
    logging.info(f"Tabela de Taxonomia: {taxonomy_table_path}")
    logging.info(f"Coluna de Assembly: {assembly_column}")
    logging.info(f"Coluna de Lineage: {lineage_column}")
    logging.info(f"Modo Dry Run: {'Ativado' if dry_run else 'Desativado'}")
    
    # Validar diretórios e arquivos de entrada
    if not os.path.isdir(antismash_dir):
        logging.error(f"O diretório antiSMASH '{antismash_dir}' não existe.")
        sys.exit(1)
    if not os.path.isfile(taxonomy_table_path):
        logging.error(f"A tabela de taxonomia '{taxonomy_table_path}' não existe.")
        sys.exit(1)
    if not dry_run:
        os.makedirs(bigslice_dir, exist_ok=True)
    
    # Carregar e processar dados de taxonomia
    df_tax = load_taxonomy(taxonomy_table_path, assembly_column, lineage_column)
    
    # Criar mapeamentos
    order_to_dataset, accession_to_order = create_dataset_mapping(df_tax, assembly_column)
    
    # Inicializar datasets.tsv
    datasets_tsv_path = initialize_datasets_tsv(bigslice_dir)
    
    # Processar diretórios do antiSMASH
    taxonomy_data = process_antismash_directories(
        antismash_dir=antismash_dir,
        bigslice_dir=bigslice_dir,
        dry_run=dry_run,
        order_to_dataset=order_to_dataset,
        accession_to_order=accession_to_order,
        df_tax=df_tax,
        assembly_column=assembly_column,
        lineage_column=lineage_column
    )
    
    # Gerar arquivos de taxonomia e atualizar datasets.tsv
    if not dry_run:
        generate_taxonomy_files(
            taxonomy_data=taxonomy_data,
            taxonomy_folder=os.path.join(bigslice_dir, "taxonomy"),
            bigslice_dir=bigslice_dir,
            dry_run=dry_run
        )
    
    logging.info("Organização de datasets concluída com sucesso!")
    logging.info(f"Arquivo 'datasets.tsv' criado em: {datasets_tsv_path}")

if __name__ == "__main__":
    main()


