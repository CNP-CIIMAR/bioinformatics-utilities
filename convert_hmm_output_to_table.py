import sys
import pandas as pd


def read_hmmer(path, program="hmmsearch", format="tblout", add_header_as_index=False, verbose=True):
    """
    Developed using HMMER v3.2.1
    # =========
    # tblout
    # =========
    (1) target name: The name of the target sequence or profile.
    (2) accession: The accession of the target sequence or profile, or â€™-â€™ if none.
    ...
    ...
    (19) description of target: The remainder of the line is the targetâ€™s description line, as free text.
    """
    assert_acceptable_arguments(program, {"hmmsearch", "hmmscan"})
    assert_acceptable_arguments(format, {"tblout", "domtblout"})

    if format in {"tblout", "domtblout"}:
        cut_index = {"tblout": 18, "domtblout": 22}[format]
        data = []
        header = []
        with open(path, "r") as f:
            for line in f.readlines():
                if line.startswith("#"):
                    header.append(line)
                else:
                    row = list(filter(bool, line.strip().split(" ")))
                    row = row[:cut_index] + [" ".join(row[cut_index:])]
                    data.append(row)

        df = pd.DataFrame(data)
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
            if format == "domtblout":
                columns = list(
                    map(
                        lambda field: ("identifier", field),
                        [
                            "target_name",
                            "target_accession",
                            "target_length",
                            "query_name",
                            "query_accession",
                            "query_length",
                        ],
                    )
                ) + list(
                    map(
                        lambda field: ("full_sequence", field),
                        ["e-value", "score", "bias"],
                    )
                ) + list(
                    map(
                        lambda field: ("this_domain", field),
                        [
                            "domain_number",
                            "total_domains",
                            "e-value",
                            "i-value",
                            "score",
                            "bias",
                        ],
                    )
                ) + list(
                    map(
                        lambda field: ("hmm_coord", field),
                        ["from", "to"],
                    )
                ) + list(
                    map(
                        lambda field: ("ali_coord", field),
                        ["from", "to"],
                    )
                ) + list(
                    map(
                        lambda field: ("env_coord", field),
                        ["from", "to"],
                    )
                ) + list(
                    map(
                        lambda field: ("identifier", field),
                        ["acc", "query_description"],
                    )
                )

            df.columns = pd.MultiIndex.from_tuples(columns)
        else:
            df = pd.DataFrame(columns=pd.MultiIndex.from_tuples(columns))
    else:
        raise NotImplementedError(f"Format '{format}' not implemented for reading.")

    if add_header_as_index:
        df.set_index(header, append=True, inplace=True)
        df.index.names = ["header", "index"]

    # Extrair a espécie de cada linha da coluna "query_description" usando expressões regulares
    df['species'] = df[('identifier', 'query_description')].str.extract(r'\[(.*?)\]')

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
    if len(sys.argv) != 3:
        print("Uso: python script.py <input_path> <output_path>")
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = sys.argv[2]

    convert_hmmer_table(input_path, output_path)
