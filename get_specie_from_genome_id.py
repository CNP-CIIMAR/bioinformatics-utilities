import subprocess
import sys
import os
from datetime import datetime

def read_genome_ids(file_path):
    """
    Lê os IDs dos genomas de um arquivo e os retorna como uma lista.
    """
    with open(file_path, 'r') as file:
        return [line.strip() for line in file]

def extract_prefix(genome_id):
    """
    Extrai o prefixo que precede o segundo `_` de um ID de genoma.
    """
    parts = genome_id.split('_')
    if len(parts) > 2:
        return '_'.join(parts[:2])
    return genome_id

def generate_intermediate_list(genome_list, intermediate_file_name):
    """
    Gera uma lista intermediária de IDs extraindo o prefixo até o segundo `_`.
    """
    with open(intermediate_file_name, 'w') as file:
        for genome_id in genome_list:
            prefix = extract_prefix(genome_id)
            file.write(prefix + '\n')

def generate_tsv(genome_list, temp_file_name):
    # Escreve o cabeçalho antes do loop para garantir que seja adicionado apenas uma vez.
    with open(temp_file_name, 'w') as temp_file:
        temp_file.write("Assembly Accession\tOrganism Common Name\tOrganism Name\n")
        for genome_id in genome_list:
            # Comando base para obter informações de cada genoma
            ### É necessario baixar os dois programas: datasets and dataformat do ncbi em https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/
            base_cmd = f"./datasets summary genome accession {genome_id} --as-json-lines | /path/to/dataformat tsv genome --fields accession,organism-common-name,organism-name"
            
            # Executa o comando e captura a saída
            process = subprocess.Popen(base_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()
            
            if process.returncode != 0:
                print(f"Erro ao executar o comando para o ID {genome_id}: {stderr.decode('utf-8')}")
                continue  # Pula para o próximo ID em caso de erro

            # Processa a saída para cada genoma
            output = stdout.decode('utf-8').strip()
            if output:
                temp_file.write(output + '\n')

def finalize_output(temp_file_name, output_file_name):
    # Copia a primeira linha do arquivo temporário para o arquivo final
    with open(temp_file_name, 'r') as temp_file, open(output_file_name, 'w') as output_file:
        output_file.write(temp_file.readline())

    # Usa o comando grep para filtrar as linhas que começam com 'GC' e anexa ao arquivo final
    grep_cmd = f"grep '^GC' {temp_file_name} >> {output_file_name}"
    subprocess.run(grep_cmd, shell=True)

    # Remove o arquivo temporário
    os.remove(temp_file_name)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Uso: script.py <genome_list_file> <output_file_name>")
        sys.exit(1)
    
    genome_list_file = sys.argv[1]
    output_file_name = sys.argv[2]
    temp_file_name = f"temp_output_{datetime.now().strftime('%Y%m%d_%H%M%S')}.tsv"
    intermediate_file_name = f"intermediate_list_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"

    genome_list = read_genome_ids(genome_list_file)
    generate_intermediate_list(genome_list, intermediate_file_name)
    intermediate_genome_list = read_genome_ids(intermediate_file_name)
    generate_tsv(intermediate_genome_list, temp_file_name)
    finalize_output(temp_file_name, output_file_name)

    # Remove o arquivo intermediário
    os.remove(intermediate_file_name)

