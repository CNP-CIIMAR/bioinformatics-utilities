import os
import sys
import subprocess
import logging
import zipfile
import shutil
import argparse
import glob

def extract_assembly_accession(assembly_id):
    """
    Extrai o assembly accession (prefixo até o segundo '_') do assembly ID completo.
    
    Exemplo:
        Input: 'GCA_000523875.1_ASM52387v1_genomic.fna'
        Output: 'GCA_000523875.1'
    """
    # Remove a extensão '.fna', '.faa' ou '.gbff' se presente
    for suffix in ['.fna', '.faa', '.gbff']:
        if assembly_id.endswith(suffix):
            assembly_id = assembly_id[:-len(suffix)]
            break
    
    parts = assembly_id.split('_')
    if len(parts) < 2:
        logging.error(f"Assembly ID '{assembly_id}' não possui pelo menos dois '_' para extração.")
        return None
    
    # Junta as primeiras duas partes para obter o assembly accession
    return '_'.join(parts[:2])

def setup_logging(log_file):
    """
    Configura o sistema de logging para registrar no console e em um arquivo.
    """
    # Garantir que o diretório para o log exista
    log_dir = os.path.dirname(log_file)
    os.makedirs(log_dir, exist_ok=True)
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file, mode='a'),
            logging.StreamHandler(sys.stdout)
        ]
    )

def map_data_types_to_include(data_types):
    """
    Mapeia os tipos de dados especificados para a opção '--include' do datasets.
    
    Parâmetros:
        data_types (list): Lista de tipos de dados ('genome', 'protein', 'cds')
        
    Retorna:
        list: Lista de opções correspondentes para '--include'
    """
    mapping = {
        'genome': 'genome',      # Para arquivos .fna
        'protein': 'protein',    # Para arquivos .faa
        'cds': 'cds'             # Para arquivos CDS
    }
    include_options = []
    for data_type in data_types:
        mapped = mapping.get(data_type.lower())
        if mapped:
            include_options.append(mapped)
        else:
            logging.warning(f"Tipo de dado desconhecido ou não suportado: {data_type}")
    return include_options

def download_genomes(genome_ids, output_dir, datasets_cmd, genomes_obtained_dir, include_options, assembly_level):
    """
    Faz o download dos genomas/proteínas/cds, descompacta os arquivos, move os .fna/.faa/.gbff para o diretório final,
    e mantém registros de sucessos e falhas.
    
    Parâmetros:
        genome_ids (list): Lista de assembly IDs.
        output_dir (str): Diretório de saída para downloads e descompactação.
        datasets_cmd (str): Caminho para o executável 'datasets'.
        genomes_obtained_dir (str): Diretório final para armazenar os arquivos extraídos.
        include_options (list): Lista de opções para '--include'.
        assembly_level (str or None): Níveis de montagem a serem aplicados (apenas se 'genome' estiver incluído).
    """
    # Diretório para armazenar os arquivos extraídos
    os.makedirs(genomes_obtained_dir, exist_ok=True)
    logging.info(f"Diretório 'genomes_obtained' criado ou já existe: {genomes_obtained_dir}")
    
    # Listas para rastrear sucessos e falhas
    success_downloads = []
    failed_downloads = []
    
    for original_id in genome_ids:
        clean_id = extract_assembly_accession(original_id)
        if not clean_id:
            logging.warning(f"Pular download para assembly ID inválido: {original_id}")
            failed_downloads.append(original_id)
            continue

        zip_filename = f"{clean_id}.zip"
        zip_path = os.path.join(output_dir, zip_filename)
        
        cmd = [
            datasets_cmd, 'download', 'genome', 'accession', clean_id,
            '--include', ','.join(include_options),  # Combina múltiplas opções
            '--filename', zip_path
        ]

        # Adicionar a opção --assembly-level se 'genome' estiver incluído e assembly_level estiver especificado
        if 'genome' in include_options and assembly_level:
            cmd.extend(['--assembly-level', assembly_level])
            logging.info(f"Aplicando '--assembly-level {assembly_level}' para o assembly '{clean_id}'")
        
        try:
            logging.info(f"Iniciando download para: {clean_id} com opção '--include {','.join(include_options)}'")
            subprocess.run(cmd, check=True)
            logging.info(f"Baixado com sucesso: {zip_path}")
            
            # Descompactar o arquivo zip
            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                zip_ref.extractall(output_dir)
            logging.info(f"Descompactado: {zip_path}")
            
            # Definir os padrões de busca para os arquivos extraídos com base nas opções incluídas
            extracted_files = {}
            if 'genome' in include_options:
                # Busca por arquivos que terminam com '_genomic.fna'
                genome_pattern = os.path.join(output_dir, 'ncbi_dataset', 'data', clean_id, '*_genomic.fna')
                genome_matches = glob.glob(genome_pattern)
                if genome_matches:
                    # Assume que há apenas um arquivo genome por assembly
                    extracted_files['genome'] = genome_matches[0]
                else:
                    logging.error(f"Arquivo genome não encontrado para {clean_id} no padrão: {genome_pattern}")
                    failed_downloads.append(original_id)
                    continue
            
            if 'protein' in include_options:
                # Busca por 'protein.faa' no diretório específico
                protein_path = os.path.join(output_dir, 'ncbi_dataset', 'data', clean_id, 'protein.faa')
                if os.path.isfile(protein_path):
                    extracted_files['protein'] = protein_path
                else:
                    logging.error(f"Arquivo 'protein.faa' não encontrado para {clean_id} em: {protein_path}")
                    failed_downloads.append(original_id)
                    continue
            
            if 'cds' in include_options:
                # Busca por 'cds_from_genomic.fna' no diretório específico
                cds_path = os.path.join(output_dir, 'ncbi_dataset', 'data', clean_id, 'cds_from_genomic.fna')
                if os.path.isfile(cds_path):
                    extracted_files['cds'] = cds_path
                else:
                    logging.error(f"Arquivo 'cds_from_genomic.fna' não encontrado para {clean_id} em: {cds_path}")
                    failed_downloads.append(original_id)
                    continue
            
            # Renomear e mover os arquivos extraídos para o diretório final
            for data_type, src_path in extracted_files.items():
                if data_type == 'genome':
                    final_filename = f"{clean_id}_genomic.fna"
                elif data_type == 'protein':
                    final_filename = f"{clean_id}.faa"
                elif data_type == 'cds':
                    final_filename = f"{clean_id}_cds.fna"
                else:
                    final_filename = os.path.basename(src_path)  # Caso inesperado
                
                final_path = os.path.join(genomes_obtained_dir, final_filename)
                shutil.move(src_path, final_path)
                logging.info(f"Mover {src_path} para {final_path}")
            
            success_downloads.append(clean_id)
        
        except subprocess.CalledProcessError as e:
            logging.error(f"Falha ao baixar {clean_id}, erro: {str(e)}")
            failed_downloads.append(original_id)
        except zipfile.BadZipFile as e:
            logging.error(f"Falha ao descompactar {zip_path}, erro: {str(e)}")
            failed_downloads.append(original_id)
        except Exception as e:
            logging.error(f"Erro inesperado com {clean_id}, erro: {str(e)}")
            failed_downloads.append(original_id)
    
    # Gerar o relatório
    report_path = os.path.join(output_dir, 'download_report.txt')
    with open(report_path, 'w') as report_file:
        report_file.write("=== Relatório de Download de Genomas/Proteínas/CDS ===\n\n")
        report_file.write(f"Total de itens solicitados: {len(genome_ids)}\n")
        report_file.write(f"Baixados com sucesso: {len(success_downloads)}\n")
        report_file.write(f"Falharam no download: {len(failed_downloads)}\n\n")
        
        report_file.write("**Lista de Itens Baixados com Sucesso:**\n")
        for gid in success_downloads:
            report_file.write(f"{gid}\n")
        
        report_file.write("\n**Lista de Itens que Falharam no Download:**\n")
        for gid in failed_downloads:
            report_file.write(f"{gid}\n")
    
    logging.info(f"Relatório de download gerado em: {report_path}")

def main():
    # Definir internamente o caminho para o comando 'datasets'
    DATASETS_CMD = "/usr/local/bin/datasets"  # Atualize conforme necessário
    
    # Configurar o parser de argumentos
    parser = argparse.ArgumentParser(
        description=(
            "Script para baixar genomas, proteínas ou CDS do NCBI usando o comando 'datasets'.\n\n"
            "Uso:\n"
            "  python3 download_ncbi_datasetsv2.py <genome_ids_file> <output_dir> [--include <types>] [--assembly-level <levels>]\n\n"
            "Exemplos:\n"
            "  - Baixar apenas genomas:\n"
            "      python3 download_ncbi_datasetsv2.py genome_ids.txt res\n\n"
            "  - Baixar genomas e proteínas:\n"
            "      python3 download_ncbi_datasetsv2.py genome_ids.txt res --include genome protein\n\n"
            "  - Baixar apenas proteínas:\n"
            "      python3 download_ncbi_datasetsv2.py genome_ids.txt res --include protein\n\n"
            "  - Baixar genomas com níveis de montagem específicos:\n"
            "      python3 download_ncbi_datasetsv2.py genome_ids.txt res --include genome --assembly-level chromosome,complete"
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument("genome_ids_file", help="Caminho para o arquivo contendo os IDs de genoma.")
    parser.add_argument("output_dir", help="Diretório de saída onde os arquivos serão baixados e processados.")
    
    # Adicionar argumento para especificar os tipos de dados a serem incluídos
    parser.add_argument(
        "--include", 
        type=str, 
        nargs='+', 
        default=["genome"], 
        help="Tipos de dados a serem incluídos: 'genome', 'protein', 'cds'. Exemplo: --include genome protein"
    )
    
    # A opção --assembly-level será incluída apenas se 'genome' estiver nos tipos de dados incluídos
    parser.add_argument(
        "--assembly-level", 
        type=str, 
        default=None, 
        help="Limitar a genomas em um ou mais níveis de montagem (comma-separated): 'chromosome', 'complete', 'contig', 'scaffold'. Exemplo: --assembly-level chromosome,complete (apenas aplicável se 'genome' estiver incluído)"
    )
    
    args = parser.parse_args()
    
    genome_ids_file = args.genome_ids_file
    output_dir = args.output_dir
    assembly_level = args.assembly_level
    include_types = args.include
    
    # Diretório final para armazenar os arquivos extraídos
    genomes_obtained_dir = os.path.join(output_dir, 'genomes_obtained')

    # Configura o logging
    log_file = os.path.join(output_dir, 'download_log.txt')
    setup_logging(log_file)

    # Verifica se o comando 'datasets' existe no caminho fornecido
    if not os.path.isfile(DATASETS_CMD) or not os.access(DATASETS_CMD, os.X_OK):
        logging.error(f"O comando 'datasets' não foi encontrado ou não é executável em: {DATASETS_CMD}")
        sys.exit(1)

    # Leia os IDs dos genomas do arquivo
    try:
        with open(genome_ids_file, 'r') as f:
            genome_ids = [line.strip() for line in f if line.strip()]
        if not genome_ids:
            logging.error("O arquivo de IDs de genoma está vazio.")
            sys.exit(1)
    except Exception as e:
        logging.error(f"Erro ao ler o arquivo de IDs de genoma: {str(e)}")
        sys.exit(1)

    # Mapeia os tipos de dados para as opções '--include'
    include_options = map_data_types_to_include(include_types)
    if not include_options:
        logging.error("Nenhum tipo de dado válido para incluir foi especificado.")
        sys.exit(1)
    
    # Log das opções de inclusão
    logging.info(f"Tipos de dados incluídos para download: {', '.join(include_options)}")

    # Faça o download dos genomas/proteínas/cds
    download_genomes(genome_ids, output_dir, DATASETS_CMD, genomes_obtained_dir, include_options, assembly_level)

if __name__ == "__main__":
    main()
