import os
import subprocess
import argparse
import sys
import pandas as pd


def main(models_dir):
    EVALUE = "0.00000000000000000001"

    # Coletar todos os arquivos .faa
    fastas = [f for f in os.listdir() if f.endswith(".faa")]

    # Coletar todos os modelos HMM
    models = [os.path.join(models_dir, f) for f in os.listdir(models_dir) if f.endswith(".hmm")]

    # Loop em cada modelo HMM
    for model in models:
        model_name = os.path.basename(model).replace(".hmm", "")

        # Loop em cada arquivo FASTA
        for fasta in fastas:
            # Extrair o número de acesso do nome do arquivo FASTA
            accession = "_".join(fasta.split("_")[:2])

            # Definir o nome do arquivo de saída com base no modelo e no arquivo FASTA
            output_file = f"{model_name}_{accession}_tblout.txt"

            cmd = f"hmmsearch --domE {EVALUE} --tblout {output_file} {model} {fasta}"

            try:
                subprocess.check_call(cmd, shell=True)
                
                # Converta o arquivo tblout para tsv
                tsv_output_file = output_file.replace("_tblout.txt", ".tsv")
                convert_hmmer_table(output_file, tsv_output_file)
            except subprocess.CalledProcessError as e:
                print(f"Erro ao executar: {cmd}")
                print(f"Mensagem de erro: {e.output}")

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


def convert_hmmer_table(input_path, output_path):
    df = read_hmmer(input_path)
    write_table_to_tsv(df, output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Execute hmmsearch em arquivos .faa usando modelos do diretório fornecido.")
    parser.add_argument("models_dir", help="Diretório contendo os modelos HMM (.hmm).")
    args = parser.parse_args()

    main(args.models_dir)

