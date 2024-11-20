import os
import argparse
from pathlib import Path
from Bio import SeqIO
import csv
import shutil

def contains_domain(feature, domains):
    """
    Verifica se uma feature contém um dos domínios especificados em uma variedade de qualificadores.

    Args:
        feature (SeqFeature): Uma feature do registro GenBank.
        domains (list): Lista de domínios a procurar.

    Returns:
        dict: Contagem dos domínios encontrados.
    """
    domain_count = {domain: 0 for domain in domains}
    relevant_qualifiers = ["detection_rule", "detection_rules", "NRPS_PKS", "sec_met_domain", 
                           "domains", "aSDomain", "description", "domain_id", "label"]
    for qual in relevant_qualifiers:
        if qual in feature.qualifiers:
            for value in feature.qualifiers[qual]:
                for domain in domains:
                    if domain in value:
                        domain_count[domain] += 1
    return domain_count

def extract_genome_id(subdir_name):
    """
    Extrai o ID do genoma do nome do subdiretório.
    O ID é definido como o valor entre "Result_" e o segundo "_" após "GCA_" ou "GCF_".

    Args:
        subdir_name (str): Nome do subdiretório.

    Returns:
        str: ID do genoma extraído.
    """
    parts = subdir_name.split("_")
    if len(parts) > 2 and parts[1] in {"GCA", "GCF"}:
        return f"{parts[1]}_{parts[2]}"
    return ""

def format_size(size_mb):
    """
    Formata o tamanho em Megabases e Gigabases.

    Args:
        size_mb (float): Tamanho em Megabases.

    Returns:
        tuple: Tamanho formatado em Megabases e Gigabases.
    """
    size_gb = size_mb / 1000
    return f"{size_mb:.2f} Mb", f"{size_gb:.2f} Gb"

def process_gbk_files(input_dir, log_file, search_amp_binding=False):
    """
    Processa arquivos .gbk em um diretório de entrada, verifica a presença de domínios
    e calcula o tamanho total em megabases e gigabases de cada subdiretório.
    Além disso, copia os arquivos .gbk identificados com 'CAL_domain' para um diretório filtrado
    com prefixo do Genome ID.

    Args:
        input_dir (Path): Diretório contendo subdiretórios com arquivos .gbk.
        log_file (Path): Arquivo de log para salvar informações sobre o processamento.
        search_amp_binding (bool): Se True, também procura por 'AMP-binding' na análise.
    """
    input_dir = Path(input_dir)
    log_file = Path(log_file)
    
    if not input_dir.is_dir():
        print(f"O diretório de entrada '{input_dir}' não existe ou não é um diretório válido.")
        return

    summary_data = []
    total_size_mb = 0.0
    subdirs_accessed = 0
    subdirs_with_cal = 0

    # Diretório principal para armazenar os subdiretórios filtrados com sufixo _CAL
    filtered_main_dir = input_dir.parent / "filtrados_subdir_CAL"
    filtered_main_dir.mkdir(exist_ok=True)

    # Verifica o uso de disco no diretório pai onde filtered_main_dir será criado
    try:
        total, used, free = shutil.disk_usage(input_dir.parent)
        free_mb = free / 1_000_000  # Convertendo bytes para megabases
        free_mb_str, free_gb_str = format_size(free_mb)
    except Exception as e:
        print(f"Erro ao verificar o uso de disco no diretório '{input_dir.parent}': {e}")
        return

    with open(log_file, "w") as log:
        for subdir in input_dir.iterdir():
            if subdir.is_dir() and subdir.name.startswith("Result_"):
                subdirs_accessed += 1
                genome_id = extract_genome_id(subdir.name)
                
                if not genome_id:
                    log.write(f"Subdiretório {subdir.name} ignorado (Genome ID não encontrado)\n")
                    # Mesmo que Genome ID não seja encontrado, ainda queremos incluí-lo no relatório com 0 contagens
                    # Portanto, definimos contagens como 0 e não copiamos nenhum arquivo
                    summary_entry = {
                        "Assembly": "N/A",
                        "CAL_domain": 0,
                        "AMP-binding": 0 if search_amp_binding else "N/A",
                        "Total_size_MB": "0.00 Mb",
                        "Total_size_GB": "0.00 Gb"
                    }
                    summary_data.append(summary_entry)
                    continue  # Pula para o próximo subdiretório

                gbk_files = list(subdir.glob("*.gbk"))
                if len(gbk_files) == 0:
                    log.write(f"Subdiretório {subdir.name} ignorado (nenhum arquivo .gbk encontrado)\n")
                    # Inclui o Genome ID no relatório com 0 contagens
                    summary_entry = {
                        "Assembly": genome_id,
                        "CAL_domain": 0,
                        "AMP-binding": 0 if search_amp_binding else "N/A",
                        "Total_size_MB": "0.00 Mb",
                        "Total_size_GB": "0.00 Gb"
                    }
                    summary_data.append(summary_entry)
                    continue  # Pula para o próximo subdiretório

                # Preferir arquivos com ".region" no nome
                region_gbk_files = list(subdir.glob("*.region*.gbk"))
                if region_gbk_files:
                    gbk_files = region_gbk_files

                cal_count = 0
                amp_count = 0
                subdir_size_mb = 0.0  # em megabases
                files_with_cal = 0

                # Lista para armazenar os arquivos que contêm CAL_domain
                gbk_files_to_copy = []

                for gbk_file in gbk_files:
                    try:
                        file_has_cal = False
                        file_cal_count = 0
                        file_amp_count = 0
                        file_size_bp = 0

                        for record in SeqIO.parse(gbk_file, "genbank"):
                            file_size_bp += len(record.seq)
                            for feature in record.features:
                                domains_to_search = ["CAL_domain"]
                                if search_amp_binding:
                                    domains_to_search.append("AMP-binding")
                                domain_counts = contains_domain(feature, domains_to_search)
                                if domain_counts["CAL_domain"] > 0:
                                    file_cal_count += domain_counts["CAL_domain"]
                                    file_has_cal = True
                                if search_amp_binding and domain_counts.get("AMP-binding", 0) > 0:
                                    file_amp_count += domain_counts["AMP-binding"]

                        if file_has_cal:
                            cal_count += file_cal_count
                            if search_amp_binding:
                                amp_count += file_amp_count
                            subdir_size_mb += file_size_bp / 1_000_000  # Convertendo para megabases
                            files_with_cal += 1
                            size_mb_str, size_gb_str = format_size(file_size_bp / 1_000_000)
                            log.write(f"  Arquivo com 'CAL_domain' encontrado: {gbk_file.name} (Tamanho: {size_mb_str} / {size_gb_str})\n")

                            # Adiciona o arquivo à lista de arquivos a serem copiados
                            gbk_files_to_copy.append(gbk_file)
                    except Exception as e:
                        log.write(f"Erro ao processar o arquivo {gbk_file}: {e}\n")

                # Verifica se há pelo menos um arquivo com CAL_domain
                if cal_count > 0:
                    subdirs_with_cal += 1
                    size_mb_str, size_gb_str = format_size(subdir_size_mb)
                    summary_entry = {
                        "Assembly": genome_id,
                        "CAL_domain": cal_count,
                        "AMP-binding": amp_count if search_amp_binding else "N/A",
                        "Total_size_MB": size_mb_str,
                        "Total_size_GB": size_gb_str
                    }
                    summary_data.append(summary_entry)
                    total_size_mb += subdir_size_mb
                    log.write(f"Subdiretório {subdir.name} processado: CAL_domain={cal_count}, AMP-binding={amp_count if search_amp_binding else 'N/A'}, Total_size_MB={size_mb_str} / Total_size_GB={size_gb_str}\n\n")

                    # Cria o subdiretório filtrado e copia apenas os arquivos identificados com prefixo Genome ID
                    filtered_subdir = filtered_main_dir / f"{subdir.name}_CAL"
                    filtered_subdir.mkdir(exist_ok=True)
                    for gbk_file in gbk_files_to_copy:
                        try:
                            # Novo nome do arquivo com prefixo Genome ID
                            new_filename = f"{genome_id}_{gbk_file.name}"
                            destination_file = filtered_subdir / new_filename
                            shutil.copy2(gbk_file, destination_file)
                        except Exception as e:
                            log.write(f"Erro ao copiar o arquivo {gbk_file} para {filtered_subdir}: {e}\n")
                else:
                    # Inclui o Genome ID no relatório com 0 contagens
                    summary_entry = {
                        "Assembly": genome_id,
                        "CAL_domain": 0,
                        "AMP-binding": amp_count if search_amp_binding else "N/A",
                        "Total_size_MB": "0.00 Mb",
                        "Total_size_GB": "0.00 Gb"
                    }
                    summary_data.append(summary_entry)
                    log.write(f"Subdiretório {subdir.name} não contém arquivos com 'CAL_domain'\n\n")

        # Resumo final no log
        total_size_mb_str, total_size_gb_str = format_size(total_size_mb)
        with open(log_file, "a") as log:
            log.write(f"Total de subdiretórios acessados: {subdirs_accessed}\n")
            log.write(f"Total de subdiretórios com 'CAL_domain': {subdirs_with_cal}\n")
            log.write(f"Tamanho total de todos os arquivos .gbk encontrados: {total_size_mb_str} / {total_size_gb_str}\n")
            log.write(f"Espaço disponível no disco para cópia: {free_mb_str} / {free_gb_str}\n")
            if total_size_mb <= free_mb:
                log.write("Há espaço suficiente no disco para copiar os arquivos.\n")
                print("Há espaço suficiente no disco para copiar os arquivos.")
            else:
                log.write("Não há espaço suficiente no disco para copiar os arquivos.\n")
                print("Não há espaço suficiente no disco para copiar os arquivos.")

    # Escreve a tabela de resumo em um arquivo CSV
    summary_file = input_dir / "summary.csv"
    with open(summary_file, "w", newline="") as csvfile:
        fieldnames = ["Assembly", "CAL_domain", "AMP-binding", "Total_size_MB", "Total_size_GB"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summary_data)

    print(f"Processamento concluído.")
    print(f"Total de subdiretórios acessados: {subdirs_accessed}")
    print(f"Total de subdiretórios com 'CAL_domain': {subdirs_with_cal}")
    print(f"Tamanho total de todos os arquivos .gbk encontrados: {total_size_mb_str} / {total_size_gb_str}")
    print(f"Resumo salvo em: {summary_file}")
    print(f"Log detalhado salvo em: {log_file}")
    print(f"Arquivos filtrados copiados para o diretório: {filtered_main_dir}")

def main():
    parser = argparse.ArgumentParser(description="Analisa subdiretórios com domínios CAL_domain e opcionalmente AMP-binding, calculando o tamanho total dos arquivos .gbk encontrados para verificar espaço em disco disponível. Além disso, copia os arquivos identificados para um diretório filtrado com prefixo Genome ID.")
    parser.add_argument("input_dir", type=str, help="Diretório de entrada contendo subdiretórios com arquivos .gbk.")
    parser.add_argument("log_file", type=str, help="Arquivo de log para salvar informações sobre o processamento.")
    parser.add_argument("--search-amp-binding", action="store_true", help="Se especificado, também procura por 'AMP-binding'.")

    args = parser.parse_args()

    process_gbk_files(args.input_dir, args.log_file, args.search_amp_binding)

if __name__ == "__main__":
    main()

