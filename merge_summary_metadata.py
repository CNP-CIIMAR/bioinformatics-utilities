import pandas as pd
import argparse
import os

# Leandro de Mattos Pereira
# 18 novembro de 2024
#python merge_summary_metadata.py --summary /caminho/para/summary.csv --metadata /caminho/para/Supplementary_table2.tsv --output /caminho/para/merged_summary.csv --search-amp-binding --report /caminho/para/report.txt
## Parâmetros:
## --summary: Caminho completo para o arquivo summary.csv.
## --metadata: Caminho completo para o arquivo Supplementary_table2.tsv.
## --output: Caminho onde o arquivo merged_summary.csv será salvo.
## --search-amp-binding: (Opcional) Inclui a coluna AMP-binding no relatório final.
## --report: (Opcional) Caminho onde o relatório de correspondências será salvo.

def merge_summary_with_metadata(summary_csv, metadata_tsv, output_csv, report_file, search_amp_binding=False):
    """
    Mescla o arquivo summary.csv com os metadados do Supplementary_table2.tsv com base no Assembly Accession.
    Gera um relatório de quantas correspondências possuem informações completas para Bacteria.

    Args:
        summary_csv (str): Caminho para o arquivo summary.csv.
        metadata_tsv (str): Caminho para o arquivo Supplementary_table2.tsv.
        output_csv (str): Caminho para salvar o arquivo CSV mesclado.
        report_file (str or None): Caminho para salvar o relatório de contagem. Se None, não gera relatório.
        search_amp_binding (bool): Flag indicando se a coluna AMP-binding deve ser incluída.
    """
    # Verifica se os arquivos de entrada existem
    if not os.path.isfile(summary_csv):
        print(f"Erro: O arquivo summary.csv '{summary_csv}' não foi encontrado.")
        return
    if not os.path.isfile(metadata_tsv):
        print(f"Erro: O arquivo Supplementary_table2.tsv '{metadata_tsv}' não foi encontrado.")
        return

    # Lê o arquivo summary.csv
    try:
        summary_df = pd.read_csv(summary_csv)
    except Exception as e:
        print(f"Erro ao ler o arquivo summary.csv: {e}")
        return

    # Lê o arquivo Supplementary_table2.tsv
    try:
        metadata_df = pd.read_csv(metadata_tsv, sep='\t')
    except Exception as e:
        print(f"Erro ao ler o arquivo Supplementary_table2.tsv: {e}")
        return

    # Renomeia a coluna 'Assembly Accession' para 'Assembly' para facilitar o merge
    if 'Assembly Accession' not in metadata_df.columns:
        print("Erro: A coluna 'Assembly Accession' não foi encontrada no arquivo de metadados.")
        return
    metadata_df.rename(columns={'Assembly Accession': 'Assembly'}, inplace=True)

    # Mescla os DataFrames com base na coluna 'Assembly'
    # Usar 'right' para manter todas as linhas de metadata_df
    merged_df = pd.merge(summary_df, metadata_df, on='Assembly', how='right')

    # Reordena as colunas para que as colunas de metadados venham após 'Assembly'
    # Supondo que as colunas de summary_df sejam: 'Assembly', 'CAL_domain', 'AMP-binding', 'Total_size_MB', 'Total_size_GB'
    summary_columns = ['CAL_domain', 'AMP-binding', 'Total_size_MB', 'Total_size_GB']
    # Metadados são todas as colunas que não estão em summary_columns mais 'Assembly'
    metadata_columns = [col for col in metadata_df.columns if col != 'Assembly']
    # Define a nova ordem: 'Assembly', metadados, resumo
    new_column_order = ['Assembly'] + metadata_columns + summary_columns
    # Verifica se todas as colunas existem antes de reordenar
    existing_columns = [col for col in new_column_order if col in merged_df.columns]
    merged_df = merged_df[existing_columns]

    # Se search_amp_binding for False, remove a coluna 'AMP-binding'
    if not search_amp_binding:
        if 'AMP-binding' in merged_df.columns:
            merged_df.drop(columns=['AMP-binding'], inplace=True)

    # Salva o DataFrame mesclado em um novo arquivo CSV
    try:
        merged_df.to_csv(output_csv, index=False)
        print(f"Arquivo mesclado salvo em: {output_csv}")
    except Exception as e:
        print(f"Erro ao salvar o arquivo mesclado: {e}")
        return

    # Gera o relatório se o caminho do relatório for fornecido
    if report_file:
        try:
            # Verifica se a coluna 'Lineage' existe
            if 'Lineage' not in merged_df.columns:
                print("Erro: A coluna 'Lineage' não foi encontrada no DataFrame mesclado.")
                return

            # Filtra linhas onde 'Lineage' contém 'Bacteria'
            bacteria_df = merged_df[merged_df['Lineage'].str.contains('Bacteria', case=False, na=False)]

            # Define as colunas requeridas para a contagem
            required_columns = ['Assembly', 'CAL_domain']
            if search_amp_binding and 'AMP-binding' in merged_df.columns:
                required_columns.append('AMP-binding')

            # Verifica se as colunas requeridas existem
            missing_columns = [col for col in required_columns if col not in merged_df.columns]
            if missing_columns:
                print(f"Erro: As seguintes colunas necessárias para o relatório estão faltando: {missing_columns}")
                return

            # Condição para correspondência completa:
            # - 'Assembly' não é nulo
            # - 'CAL_domain' > 0
            # - Se 'AMP-binding' for requerido, 'AMP-binding' > 0
            conditions = (
                ~bacteria_df['Assembly'].isna() &
                (bacteria_df['CAL_domain'] > 0)
            )
            if search_amp_binding:
                conditions &= (bacteria_df['AMP-binding'] > 0)

            # Conta as correspondências completas
            complete_count = bacteria_df[conditions].shape[0]
            total_bacteria = bacteria_df.shape[0]
            incomplete_count = total_bacteria - complete_count

            # Lista de Assembly IDs onde CAL_domain == 0
            zero_cal_domain_df = bacteria_df[bacteria_df['CAL_domain'] == 0]
            zero_cal_domain_assemblies = zero_cal_domain_df['Assembly'].tolist()

            # Cria o conteúdo do relatório
            report_content = (
                f"Relatório de Correspondências para Bacteria\n"
                f"-----------------------------------------\n"
                f"Total de linhas com 'Lineage' contendo 'Bacteria': {total_bacteria}\n"
                f"Quantidade de correspondências com informações completas {tuple(required_columns)}: {complete_count}\n"
                f"Quantidade de correspondências com informações incompletas: {incomplete_count}\n\n"
                f"Lista de Assembly IDs com 'CAL_domain' igual a 0:\n"
            )

            if zero_cal_domain_assemblies:
                for assembly in zero_cal_domain_assemblies:
                    report_content += f"{assembly}\n"
            else:
                report_content += "Nenhuma correspondência com 'CAL_domain' igual a 0.\n"

            # Escreve o relatório no arquivo especificado
            with open(report_file, 'w') as f:
                f.write(report_content)

            print(f"Relatório salvo em: {report_file}")
        except Exception as e:
            print(f"Erro ao gerar o relatório: {e}")
            return

def main():
    parser = argparse.ArgumentParser(description="Mescla o summary.csv com os metadados do Supplementary_table2.tsv com base no Assembly Accession e gera um relatório de correspondências para Bacteria.")
    parser.add_argument("--summary", type=str, required=True, help="Caminho para o arquivo summary.csv.")
    parser.add_argument("--metadata", type=str, required=True, help="Caminho para o arquivo Supplementary_table2.tsv.")
    parser.add_argument("--output", type=str, required=True, help="Caminho para salvar o arquivo CSV mesclado.")
    parser.add_argument("--search-amp-binding", action="store_true", help="Se especificado, inclui a coluna 'AMP-binding' no relatório final.")
    parser.add_argument("--report", type=str, default=None, help="(Opcional) Caminho para salvar o relatório de correspondências.")

    args = parser.parse_args()

    merge_summary_with_metadata(
        summary_csv=args.summary,
        metadata_tsv=args.metadata,
        output_csv=args.output,
        report_file=args.report,
        search_amp_binding=args.search_amp_binding
    )

if __name__ == "__main__":
    main()
