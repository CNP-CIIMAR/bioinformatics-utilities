## Autor Leandro de Mattos Pereira
## CNP - Team 17/11/2023
## Portugal - CIIMAR - Matosinhos
import sys

def read_tsv(tsv_file):
    """
    Lê um arquivo TSV e extrai os IDs da primeira coluna.
    :param tsv_file: Caminho do arquivo TSV.
    :return: Lista de IDs.
    """
    with open(tsv_file, 'r') as file:
        # Pula o cabeçalho e extrai os IDs
        next(file)
        ids = [line.split('\t')[0] for line in file]
    return ids

def extract_fasta_sequences(ids, fasta_file):
    """
    Extrai sequências do arquivo FASTA com base nos IDs fornecidos.
    :param ids: Lista de IDs para extrair as sequências.
    :param fasta_file: Caminho do arquivo FASTA.
    :return: Dicionário de sequências, com IDs como chaves.
    """
    sequences = {}
    with open(fasta_file, 'r') as file:
        sequence_id = ''
        for line in file:
            if line.startswith('>'):
                sequence_id = line.strip().split()[0][1:]
                if sequence_id in ids:
                    sequences[sequence_id] = line
            elif sequence_id in ids:
                sequences[sequence_id] += line
    return sequences

def write_fasta(sequences, output_file):
    """
    Escreve todas as sequências encontradas em um único arquivo FASTA.
    :param sequences: Dicionário de sequências a serem escritas.
    :param output_file: Caminho do arquivo FASTA de saída.
    """
    with open(output_file, 'w') as file:
        for sequence in sequences.values():
            file.write(sequence)

def help():
    """
    Imprime informações de ajuda sobre como usar o script.
    """
    print("Uso: python extract_fasta.py [TSV_FILE] [FASTA_FILE] [OUTPUT_FILE]")
    print("\nTSV_FILE: Caminho para o arquivo TSV contendo os IDs na primeira coluna.")
    print("FASTA_FILE: Caminho para o arquivo FASTA de onde as sequências serão extraídas.")
    print("OUTPUT_FILE: Caminho para o arquivo FASTA de saída que conterá as sequências extraídas.")

def main(tsv_file, fasta_file, output_file):
    """
    Função principal para processar o arquivo TSV e extrair as sequências do arquivo FASTA.
    :param tsv_file: Caminho do arquivo TSV.
    :param fasta_file: Caminho do arquivo FASTA.
    :param output_file: Caminho do arquivo de saída.
    """
    ids = read_tsv(tsv_file)
    sequences = extract_fasta_sequences(ids, fasta_file)
    write_fasta(sequences, output_file)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        help()
    else:
        tsv_file = sys.argv[1]
        fasta_file = sys.argv[2]
        output_file = sys.argv[3]
        main(tsv_file, fasta_file, output_file)
