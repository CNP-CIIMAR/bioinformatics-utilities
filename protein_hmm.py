import os
import subprocess
import argparse
import sys
import pandas as pd

def main(models_dir, fastas_dir, output_dir):
    EVALUE = "0.00000000000000000001"

    # Verificar e criar o diretório de saída, se necessário
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Coletar todos os arquivos .faa
    fastas = [os.path.join(fastas_dir, f) for f in os.listdir(fastas_dir) if f.endswith(".faa")]

    # Coletar todos os modelos HMM
    models = [os.path.join(models_dir, f) for f in os.listdir(models_dir) if f.lower().endswith(".hmm")]

    # Loop em cada modelo HMM
    for model in models:
        model_name = os.path.basename(model).replace(".HMM", "")
        all_results = []

        # Loop em cada arquivo FASTA
        for fasta in fastas:
            # Extrair o número de acesso do genoma do nome do arquivo FASTA
            genome_accession = os.path.basename(fasta).split("_", 2)[:2]
            genome_accession = "_".join(genome_accession)

            # Definir o nome do arquivo de saída com base no modelo e no arquivo FASTA
            output_file = os.path.join(output_dir, f"{model_name}_{genome_accession}_tblout.txt")

            # Incluir o número de acesso do genoma como um comentário no início do arquivo
            with open(output_file, "w") as out_f:
                out_f.write(f"# Genome_accession: {genome_accession}\n")

            cmd = f"hmmsearch --domE {EVALUE} --tblout {output_file} {model} {fasta}"

            try:
                subprocess.check_call(cmd, shell=True)
                df = read_hmmer(output_file)
                
                # Adicionar o número de acesso do genoma como uma coluna na tabela
                df['Genome_accession'] = genome_accession
                
                all_results.append(df)
            except subprocess.CalledProcessError as e:
                print(f"Erro ao executar: {cmd}")
                print(f"Mensagem de erro: {e.output}")

        # Após o loop de fasta, mas ainda dentro do loop de modelo:
        if all_results:
            # Concatenar todos os resultados em um único DataFrame
            final_df = pd.concat(all_results, ignore_index=True)
            
            # Definir o nome do arquivo de saída para a tabela final concatenada
            final_output_file = os.path.join(output_dir, f"{model_name}_all_results.tsv")
            
            # Escrever a tabela final concatenada com o número de acesso do genoma
            write_table_to_tsv(final_df, final_output_file)

    print("Busca e conversão concluídas!")


def read_hmmer(path, program="hmmsearch", format="tblout", add_header_as_index=False, verbose=True):
    assert_acceptable_arguments(program, {"hmmsearch", "hmmscan"})
    assert_acceptable_arguments(format, {"tblout", "domtblout"})

    # Criar um dataframe vazio para 'tblout'
    empty_df = pd.DataFrame(columns=pd.MultiIndex.from_tuples(
        [
            ("identifier", "target_name"),
            ("identifier", "target_accession"),
            ("identifier", "query_name"),
            ("identifier", "query_accession"),
            ("full_sequence", "e-value"),
            ("full_sequence", "score"),
            ("full_sequence", "bias"),
            ("best_domain", "e-value"),
            ("best_domain", "score"),
            ("best_domain", "bias"),
            ("domain_number_estimation", "exp"),
            ("domain_number_estimation", "reg"),
            ("domain_number_estimation", "clu"),
            ("domain_number_estimation", "ov"),
            ("domain_number_estimation", "env"),
            ("domain_number_estimation", "dom"),
            ("domain_number_estimation", "rep"),
            ("domain_number_estimation", "inc"),
            ("identifier", "query_description")
        ]
    ))

    data = []
    header = []
    with open(path, "r") as f:
        for line in f.readlines():
            if line.startswith("#"):
                header.append(line)
            else:
                row = list(filter(bool, line.strip().split(" ")))
                row = row[:18] + [" ".join(row[18:])]
                data.append(row)

    if not data:
        return empty_df

    df = pd.DataFrame(data, columns=empty_df.columns)

    if not df.empty:
        if format == "tblout":
            columns = list(
                map(
                    lambda field: ("identifier", field),
                    [
                        "target_name",
                        "target_accession",
                        "query_name",
                        "query_accession",
                    ],
                )
            ) + list(
                map(
                    lambda field: ("full_sequence", field),
                    ["e-value", "score", "bias"],
                )
            ) + list(
                map(
                    lambda field: ("best_domain", field),
                    ["e-value", "score", "bias"],
                )
            ) + list(
                map(
                    lambda field: ("domain_number_estimation", field),
                    ["exp", "reg", "clu", "ov", "env", "dom", "rep", "inc"],
                )
            ) + list(
                map(
                    lambda field: ("identifier", field),
                    ["query_description"],
                )
            )

            df.columns = pd.MultiIndex.from_tuples(columns)
    else:
        raise ValueError(f"No data found in the file for format '{format}'.")
                

    if add_header_as_index:
        df.set_index(header, append=True, inplace=True)
        df.index.names = ["header", "index"]

    # Extrair a espécie de cada linha da coluna "query_description" usando expressões regulares
    df['Specie'] = df[('identifier', 'query_description')].str.extract(r'\[(.*?)\]')

    # Verifique se a extração foi bem-sucedida e preencha os valores ausentes
    if df['Specie'].isna().any():
        print("Warning: Some entries did not have species information extracted successfully.")
        df['Specie'].fillna("Unknown", inplace=True)

    return df


def assert_acceptable_arguments(argument, acceptable_values):
    if argument not in acceptable_values:
        raise ValueError(
            f"Invalid argument '{argument}'. Acceptable values are: {acceptable_values}"
        )


def write_table_to_tsv(df, output_path):
    df.to_csv(output_path, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Execute hmmsearch em arquivos .faa usando modelos do diretório fornecido.")
    parser.add_argument("models_dir", help="Diretório contendo os modelos HMM (.hmm).")
    parser.add_argument("fastas_dir", help="Diretório contendo os arquivos FASTA (.faa).")
    parser.add_argument("output_dir", help="Diretório para salvar os resultados.")
    args = parser.parse_args()

    main(args.models_dir, args.fastas_dir, args.output_dir)




