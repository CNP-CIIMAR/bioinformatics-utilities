## Date: 05/07/2023
## Authors: 
## Leandro de Mattos Pereira, Junior Researcher
### Pedro Leao - Team Leader Researcher - CNP Team - Leao Lab.
import os
import subprocess

def Cluster(file_path, threshold, statistics=None):
    print(f"++++++++++> work on {file_path}")
    ref_hash = initial_matrix(file_path)
    print("Do blast...")
    format_db_cmd = f"formatdb -i {file_path} -o T"
    subprocess.check_call(format_db_cmd, shell=True)

    blast_cmd = f"blastall -b 10000000 -p blastp -e 0.1 -d {file_path} -i {file_path} -o {file_path}.tmp"
    subprocess.check_call(blast_cmd, shell=True)
    print(f"Did blast of {file_path}")

    print("Start the parsing...")
    ref_hash = parse_score(ref_hash, f"{file_path}.tmp")
    print("Parsing done")

    ref_hash = balance_hash(ref_hash)
    print("Do interpretation...")
    ref_group = interpretate_hash(ref_hash)
    print("Print Result...")
    write_result(ref_group, file_path)
    print(f"Done for {file_path}")
    os.unlink(f"{file_path}.tmp")

def initial_matrix(file_path):
    ref_hash = {}
    sequence_amount, ref_names = get_amount_seq(file_path)
    amount = len(ref_names)
    for i in range(amount):
        ref_hash[ref_names[i]] = {}
        for j in range(amount):
            ref_hash[ref_names[i]][ref_names[j]] = 0
    return ref_hash

def balance_hash(ref_hash):
    for k1 in ref_hash:
        for k2 in ref_hash[k1]:
            if (ref_hash[k1][k2] is not None and ref_hash[k2][k1] is not None and
                    ref_hash[k1][k2] != ref_hash[k2][k1]):
                ref_hash[k1][k2] = 1
                ref_hash[k2][k1] = 1
                print("not symmetric...")
            elif ((ref_hash[k1][k2] is None and ref_hash[k2][k1] is not None) or
                  (ref_hash[k1][k2] is not None and ref_hash[k2][k1] is None)):
                ref_hash[k1][k2] = 1
                ref_hash[k2][k1] = 1
                print("not symmetric...")

    print("Balanced Matrix")
    return ref_hash

def parse_blast(ref_hash, file_name):
    with open(file_name, 'r') as file:
        for line in file:
            line = line.rstrip('\n')
            qName, sName, qlength, slength, alignmentlength, ident, sim, Evalue, ScoreBit = line.split('\t')
            if float(ScoreBit) >= Threshold:
                ref_hash[qName][sName] = 1

    print("Did the parsing on", file_name)
    return ref_hash


def parse_score(ref_hash, file_name):
    with open(file_name, 'r') as file:
        query = None
        subject = None
        for line in file:
            line = line.rstrip('\n')
            if line.startswith('Query='):
                query = line.split('=')[1].strip()
            elif line.startswith('>'):
                subject = line.split('>')[1].strip()
            elif line.startswith('  Score ='):
                score = float(line.split('=')[1].split('bits')[0].strip())
                if score >= Threshold:
                    ref_hash[query][subject] = 1

    print("Parsing done")
    return ref_hash

from Bio import SearchIO

def parse_blast_old(ref_hash, file_name):
    file = SearchIO.parse(file_name, 'blast')

    for result in file:
        for hit in result.hits:
            for hsp in hit.hsps:
                evalue = hsp.evalue
                score_bit = hit.bits
                if score_bit >= Threshold:
                    ref_hash[result.query_name][hit.id] = 1

    print("Did the parsing on", file_name)
    return ref_hash

def interpretate_hash(ref_hash):
    element_in_cluster = {}  # Hash com as sequências e seus respectivos clusters
    no_cluster = 0  # Número do cluster

    for key in list(ref_hash):
        if key in ref_hash[key]:  # Uma linha será inicializada devido à própria sequência
            no_cluster += 1
            ref_list = get_list([], ref_hash[key])
            del ref_hash[key]
            element = ref_list.pop()

            while element:
                element_in_cluster[element] = no_cluster  # Atribuir o elemento ao cluster correspondente
                if element in ref_hash and element in ref_hash[element]:
                    ref_list = get_list(ref_list, ref_hash[element])
                    del ref_hash[element]
                    element = ref_list.pop()
                else:
                    element = None

    return element_in_cluster

def write_result(element_in_cluster, file_name):
    amount_in_group = [0] * (max(element_in_cluster.values()) + 1)  # Para estatísticas
    cluster_file = {}  # Irá conter as sequências de cada cluster

    for key, cluster in element_in_cluster.items():
        cluster_file.setdefault(cluster, "")
        cluster_file[cluster] += get_sequence_of_fasta_file(file_name, key)
        amount_in_group[cluster] += 1

    # Salvar o resultado em arquivos
    for key, sequences in cluster_file.items():
        file_name = file_name.replace("Sequences", "")
        with open(f"Results/{file_name}_{key}", "w") as file:
            file.write(sequences)

    # Imprimir informações sobre os clusters na saída padrão
    for i in range(1, len(amount_in_group)):
        print(f"Group {i}: {amount_in_group[i]} elements.")

    # Se 'Statistics' estiver definido, um arquivo será gerado com a quantidade de clusters para cada arquivo (mais)
    if statistics is not None:
        amount = len(amount_in_group) - 1
        with open(statistics, "a") as file:
            file.write(f"{file_name} {amount}\n")

def get_list(ref_hash):
    ref_list_ar = []

    for key, value in ref_hash.items():
        if value:
            ref_list_ar.append(key)

    return ref_list_ar


def occur_error(sequence_name, work_dir, error):
    with open(work_dir + sequence_name + ".Error", "w") as file:
        file.write(f"Following Error occurred: {error}")


def get_sequence_of_fasta_file(file_name, name):
    with open(file_name, "r") as file:
        result = ""
        found = False

        for line in file:
            if line.startswith(">"):
                if found:
                    break
                if name == line[1:].strip():
                    found = True
                    result = line
            elif found:
                result += line

    if not result:
        raise ValueError(f"Sequence {name} was not found in file {file_name}")

    return result


def get_amount_seq(file_name):
    with open(file_name, "r") as file:
        counter = 0
        seq_names = []

        for line in file:
            if line.startswith(">"):
                seq_names.append(line[1:].strip())
                counter += 1

    print(f"The file {file_name} has {counter} sequences...")
    return counter, seq_names
