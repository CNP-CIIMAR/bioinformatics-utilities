import subprocess
import sys
import os
from datetime import datetime

def read_genome_ids(file_path):
    """
    LÃª os IDs dos genomas de um arquivo e os retorna como uma lista.
    """
    with open(file_path, 'r') as file:
        return [line.strip() for line in file]

def generate_tsv(genome_list, temp_file_name):
    # Escreve o cabeÃ§alho antes do loop para garantir que seja adicionado apenas uma vez.
    with open(temp_file_name, 'w') as temp_file:
        temp_file.write("Assembly Accession\tOrganism Common Name\tOrganism Name\n")
        for genome_id in genome_list:
            # Comando base para obter informaÃ§Ãµes de cada genoma
            base_cmd = f"datasets summary genome accession {genome_id} --as-json-lines | dataformat tsv genome --fields accession,organism-common-name,organism-name"
            
            # Executa o comando e captura a saÃ­da
            process = subprocess.Popen(base_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()
            
            if process.returncode != 0:
                print(f"Erro ao executar o comando para o ID {genome_id}: {stderr.decode('utf-8')}")
                continue  # Pula para o prÃ³ximo ID em caso de erro

            # Processa a saÃ­da para cada genoma
            output = stdout.decode('utf-8').strip()
            if output:
                temp_file.write(output + '\n')

def finalize_output(temp_file_name, output_file_name):
    # Copia a primeira linha do arquivo temporÃ¡rio para o arquivo final
    with open(temp_file_name, 'r') as temp_file, open(output_file_name, 'w') as output_file:
        output_file.write(temp_file.readline())

    # Usa o comando grep para filtrar as linhas que comeÃ§am com 'GC' e anexa ao arquivo final
    grep_cmd = f"grep '^GC' {temp_file_name} >> {output_file_name}"
    subprocess.run(grep_cmd, shell=True)

    # Remove o arquivo temporÃ¡rio
    os.remove(temp_file_name)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Uso: script.py <genome_list_file> <output_file_name>")
        sys.exit(1)
    
    genome_list_file = sys.argv[1]
    output_file_name = sys.argv[2]
    temp_file_name = f"temp_output_{datetime.now().strftime('%Y%m%d_%H%M%S')}.tsv"

    genome_list = read_genome_ids(genome_list_file)
    generate_tsv(genome_list, temp_file_name)
    finalize_output(temp_file_name, output_file_name)
