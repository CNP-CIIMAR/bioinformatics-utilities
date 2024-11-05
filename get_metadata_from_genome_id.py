#Autor: Leandro de Mattos Pereira
# CNP Team - Leam Lab 
import subprocess
import sys
from ete3 import NCBITaxa

# Inicializa o NCBITaxa
ncbi = NCBITaxa()

if len(sys.argv) < 3:
    print("Uso: script.py <input_file_with_assembly_ids> <output_file>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

# Função para obter a linhagem taxonômica completa com base no Organism Tax ID
def get_lineage(tax_id):
    try:
        lineage = ncbi.get_lineage(tax_id)
        names = ncbi.get_taxid_translator(lineage)
        lineage_str = "; ".join([names[t] for t in lineage])
        return lineage_str
    except Exception as e:
        print(f"Erro ao obter linhagem para Tax ID {tax_id}: {e}")
        return "Desconhecido"

# Escreve o cabeçalho manualmente no arquivo de saída, incluindo a coluna Lineage
with open(output_file, 'w') as out_file:
    out_file.write(
        "Assembly Accession\tOrganism Name\tOrganism Common Name\tOrganism Tax ID\tLineage\t"
        "Assembly Level\tBioProject Accession\tBioSample Accession\tGC Percent\t"
        "Total Sequence Length\tSequencing Technology\tRelease Date\tCollection Date\t"
        "BioSample Description\n"
    )

# Lê o arquivo de entrada e processa cada Assembly ID
with open(input_file, 'r') as in_file:
    # Pula o cabeçalho do arquivo de entrada
    next(in_file)

    for line in in_file:
        # Extrai o Assembly ID (presumindo que seja o primeiro campo)
        assembly_id = line.strip().split('\t')[0]

        # Constrói o comando para cada Assembly ID
        command = (
            f"./datasets summary genome accession {assembly_id} --as-json-lines | "
            f"./dataformat tsv genome --fields accession,organism-name,organism-common-name,"
            f"organism-tax-id,assminfo-level,assminfo-bioproject,assminfo-biosample-accession,"
            f"assmstats-gc-percent,assmstats-total-sequence-len,assminfo-sequencing-tech,"
            f"assminfo-release-date,assminfo-biosample-collection-date,"
            f"assminfo-biosample-description-title"
        )

        # Executa o comando e captura a saída
        result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Verifica se o comando foi executado com sucesso
        if result.returncode == 0:
            # Decodifica a saída e processa linha por linha
            output_lines = result.stdout.decode('utf-8').splitlines()

            # Ignora o cabeçalho do `dataformat` na primeira linha, se ele estiver presente
            if output_lines and output_lines[0].startswith("Assembly Accession"):
                output_lines = output_lines[1:]
            
            # Processa cada linha de saída do `dataformat`
            for output_line in output_lines:
                fields = output_line.split('\t')
                
                # Obtém o Organism Tax ID e calcula a linhagem
                tax_id = fields[3] if len(fields) > 3 else ""
                lineage = get_lineage(int(tax_id)) if tax_id.isdigit() else "Desconhecido"
                
                # Insere a linhagem após o Organism Tax ID e escreve no arquivo
                enriched_line = "\t".join(fields[:4] + [lineage] + fields[4:])
                with open(output_file, 'a') as out_file:
                    out_file.write(enriched_line + "\n")
        else:
            # Em caso de erro, imprime a mensagem no console
            print(f"Erro ao processar {assembly_id}: {result.stderr.decode('utf-8')}")
