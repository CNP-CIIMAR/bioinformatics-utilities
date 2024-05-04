from Bio import SeqIO
import os
import subprocess
import sys
from statistics import mean, stdev
import pandas as pd
from statistics import mean, stdev
import matplotlib.pyplot as plt
from collections import Counter
import numpy as np
from Bio import SeqIO
import pandas as pd
import os
from collections import Counter

#A tabela de códons consolidada, no contexto do seu script, é uma estrutura de dados (no caso, um dicionário Python) que acumula e soma a contagem de cada códon em todos os genomas analisados. Ela é construída agregando as contagens de códons de múltiplos arquivos de sequência de DNA (genomas) para fornecer uma visão geral do uso de códons em todo o conjunto de dados, ao invés de em um único genoma.

def determine_format(filename):
    if filename.endswith('.fna') or filename.endswith('.fasta'):
        return 'fasta'
    elif filename.endswith('.txt'):
        return 'text'  # ou simplesmente 'ignore'
    else:
        raise ValueError(f"Unknown file extension for file {filename}")


def update_codon_counts(sequence, codon_table):
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i+3]
        if codon in codon_table:
            codon_table[codon] += 1
        else:
            codon_table[codon] = 1

def calculate_codon_usage(sequence, total_codons):
    codon_usage = {}
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i+3]
        codon_usage[codon] = codon_usage.get(codon, 0) + 1
    for codon, count in codon_usage.items():
        codon_usage[codon] = count / total_codons
    return codon_usage

def generate_codon_usage_csv(codon_table, codon_to_aa, csv_file_path):
    data = []
    total_codons = sum(codon_table.values())
    
    for codon, count in codon_table.items():
        aminoacid = codon_to_aa.get(codon, '?')
        frequency = count / total_codons
        data.append({"codon": codon, "count": count, "frequency": frequency, "aminoacid": aminoacid})
    
    df = pd.DataFrame(data)
    df.to_csv(csv_file_path, index=False)
    

# Dicionário de Mapeamento de Códons para Aminoácidos
codon_to_aa = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}

def extract_cds_with_prodigal(directory_path):
    # Iterar sobre todos os arquivos no diretório fornecido
    for file in os.listdir(directory_path):
        if file.endswith('.fna') or file.endswith('.fasta'):
            file_path = os.path.join(directory_path, file)
            output_cds_path = os.path.join(directory_path, f"{file}_cds.fna")
            
            # Chamar o Prodigal para extrair CDSs
            subprocess.run(['prodigal', '-i', file_path, '-d', output_cds_path])
def extract_and_save_codon_usage(directory_path):
    consolidated_codon_table = {}
    for file in os.listdir(directory_path):
        if file.endswith("_cds.fna"):  # Assegura que apenas arquivos CDS .fna sejam processados
            file_path = os.path.join(directory_path, file)
            # Não é mais necessário chamar determine_format aqui, pois já sabemos que o arquivo é um .fna
            individual_codon_table = {}
            for record in SeqIO.parse(file_path, "fasta"):  # Podemos usar "fasta" diretamente
                sequence = str(record.seq)
                update_codon_counts(sequence, individual_codon_table)
                update_codon_counts(sequence, consolidated_codon_table)

            total_codons = sum(individual_codon_table.values())
            with open(os.path.join(directory_path, f"{file}_codon_usage.txt"), "w") as outfile:
                for codon, count in individual_codon_table.items():
                    frequency = count / total_codons
                    aminoacid = codon_to_aa.get(codon, '?')
                    outfile.write(f"{codon}: {count}, Frequency: {frequency:.4f}, Amino Acid: {aminoacid}\n")

    total_codons = sum(consolidated_codon_table.values())
    with open(os.path.join(directory_path, "consolidated_codon_usage.txt"), "w") as outfile:
        for codon, count in consolidated_codon_table.items():
            frequency = count / total_codons
            aminoacid = codon_to_aa.get(codon, '?')
            outfile.write(f"{codon}: {count}, Frequency: {frequency:.4f}, Amino Acid: {aminoacid}\n")

    return consolidated_codon_table



def analyze_codon_usage(directory_path):
    codon_counts = {}
    for file in os.listdir(directory_path):
        if "_codon_usage.txt" in file:
            file_path = os.path.join(directory_path, file)
            with open(file_path, 'r') as f:
                for line in f:
                    parts = line.strip().split(", ")
                    codon = parts[0].split(": ")[0]
                    count = int(parts[0].split(": ")[1])
                    
                    if codon not in codon_counts:
                        codon_counts[codon] = []
                    codon_counts[codon].append(count)
    return codon_counts

def plot_top_20_codons(codon_counts):
    codon_means = {codon: mean(counts) for codon, counts in codon_counts.items()}
    top_20_codons = sorted(codon_means.items(), key=lambda x: x[1], reverse=True)[:20]
    top_20_codon_keys = [codon for codon, _ in top_20_codons]
    top_20_data = {codon: codon_counts[codon] for codon in top_20_codon_keys}
    
    plt.figure(figsize=(12, 8))
    plt.boxplot(top_20_data.values(), labels=top_20_data.keys())
    plt.xticks(rotation=45)
    plt.ylabel('Count')
    plt.title('Top 20 Most Abundant Codons')

    # Salvar o gráfico em alta resolução em diferentes formatos
    for format in ['png', 'svg', 'jpeg']:
        plt.savefig(f"top_20_codons.{format}", format=format, dpi=300)
    plt.close()

from Bio import SeqIO
import os
import pandas as pd
from collections import Counter

def count_codons(sequence):
    """Conta os códons em uma sequência de DNA."""
    return Counter([sequence[i:i+3].upper() for i in range(0, len(sequence)-2, 3) if len(sequence[i:i+3]) == 3])

def gc_content(sequence):
    """Calcula a porcentagem de conteúdo GC em uma sequência de DNA."""
    gc_count = sum(1 for base in sequence.upper() if base in ['G', 'C'])
    return (gc_count / len(sequence)) * 100 if sequence else 0
from collections import Counter
import os
from Bio import SeqIO
import pandas as pd
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from Bio import SeqIO

def analyze_genomes_for_hgt(directory_path):
    all_data = []

    for filename in os.listdir(directory_path):
        # Modificação para considerar tanto arquivos de genoma quanto de CDS
        if '_genomic.fna' in filename or '_genomic.fna_cds.fna' in filename:
            file_path = os.path.join(directory_path, filename)
            # Ajuste para identificar corretamente o genome_id com novas extensões
            if '_genomic.fna_cds.fna' in filename:
                genome_id = filename.split('_genomic.fna_cds.fna')[0]
            else:
                genome_id = filename.split('_genomic.fna')[0]
            
            print(f"Opening file: {filename}")  # Mensagem de debug ao abrir cada arquivo

            all_codons_counts = Counter()
            total_at = total_cg = total_bases = 0
            for record in SeqIO.parse(file_path, "fasta"):
                seq = str(record.seq)
                all_codons_counts += count_codons(seq)
                gc = gc_content(seq) * len(seq) / 100
                at = len(seq) - gc
                total_at += at
                total_cg += gc
                total_bases += len(seq)
            
            # Evitando divisão por zero com uma verificação de segurança
            if total_bases > 0:
                genome_average_at = total_at / total_bases
                genome_average_cg = total_cg / total_bases
            else:
                print(f"No bases found in genome: {genome_id}. Skipping genome.")
                continue  # Pula para o próximo arquivo se não houver bases

            total_codons = sum(all_codons_counts.values())
            if total_codons > 0:  # Evitando divisão por zero na média dos codons
                genome_average_counts = {codon: float(count) / total_codons for codon, count in all_codons_counts.items()}
            else:
                print(f"No codons found in genome: {genome_id}. Skipping genome.")
                continue  # Pula para o próximo arquivo se não houver codons

            significant_cds = []
            for record in SeqIO.parse(file_path, "fasta"):
                cds_sequence = str(record.seq)
                cds_codon_counts = count_codons(cds_sequence)
                cds_gc = gc_content(cds_sequence) * len(cds_sequence) / 100
                cds_at = len(cds_sequence) - cds_gc

                criterion_met = []

                # Cálculo do Z-score para Codon Usage
                codon_z_scores = {codon: (cds_codon_counts[codon] - genome_average_counts.get(codon, 0)) / np.sqrt(genome_average_counts.get(codon, 0)) for codon in cds_codon_counts}
                if any(z_score > 1.96 for z_score in codon_z_scores.values()):  # Teste de Z-score com alfa = 0.05 (intervalo de confiança de 95%)
                    criterion_met.append('Codon Usage')

                # Cálculo do Z-score para AT Content
                at_z_score = (cds_at / len(cds_sequence) - genome_average_at) / np.sqrt(genome_average_at)
                if at_z_score > 1.96:
                    criterion_met.append('AT Content')

                # Cálculo do Z-score para CG Content
                gc_z_score = (cds_gc / len(cds_sequence) - genome_average_cg) / np.sqrt(genome_average_cg)
                if gc_z_score > 1.96:
                    criterion_met.append('CG Content')

                if criterion_met:
                    significant_cds.append({
                        'Genome ID': genome_id,
                        'CDS ID': record.id,
                        'CDS Sequence': cds_sequence,
                        'Criterion Met': ', '.join(criterion_met)
                    })

                     # Boxplot para Codon Usage
                 #   plt.figure(figsize=(10, 5))
                 #   plt.boxplot([list(cds_codon_counts.values()), list(genome_average_counts.values())])
                 #   plt.xticks([1, 2], ['CDS Codon Usage', 'Genome Codon Usage'])
                 #   plt.title(f'Codon Usage Distribution in CDS and Genome for {genome_id} - {record.id}')
                 #   plt.ylabel('Codon Count')
                 #   plt.savefig(f"{genome_id}_{record.id}_codon_usage_boxplot.png")
                 #   plt.close()

                    # Boxplot para GC Content
                 #   plt.figure(figsize=(10, 5))
                 #   plt.boxplot([cds_gc, genome_average_cg])
                 #   plt.xticks([1, 2], ['CDS GC Content', 'Genome GC Content'])
                 #   plt.title(f'GC Content Distribution in CDS and Genome for {genome_id} - {record.id}')
                 #   plt.ylabel('GC Content (%)')
                 #   plt.savefig(f"{genome_id}_{record.id}_gc_content_boxplot.png")
                 #   plt.close()

                    # Boxplot para AT Content
                 #   plt.figure(figsize=(10, 5))
                 #   plt.boxplot([cds_at, genome_average_at])
                 #   plt.xticks([1, 2], ['CDS AT Content', 'Genome AT Content'])
                 #   plt.title(f'AT Content Distribution in CDS and Genome for {genome_id} - {record.id}')
                 #   plt.ylabel('AT Content (%)')
                 #   plt.savefig(f"{genome_id}_{record.id}_at_content_boxplot.png")
                 #   plt.close()

            if significant_cds:
                pd.DataFrame(significant_cds).to_csv(f"{directory_path}/{genome_id}_significant_cds.csv", index=False)
                all_data.extend(significant_cds)

    if all_data:
        pd.DataFrame(all_data).to_csv(f"{directory_path}/concatenated_significant_cds.csv", index=False)

# Chamada da função com o diretório contendo os arquivos de genoma e CDS
# Replace 'your_directory_path_here' with the path to your directory containing the '_cds.fna' files.
try:
    import google.colab
    IN_COLAB = True
except:
    IN_COLAB = False

if IN_COLAB:
    from google.colab import drive
    drive.mount('/content/drive')
    print("ForneÃ§a o caminho do diretÃ³rio dentro do seu Google Drive. Exemplo: /content/drive/My Drive/genomes_directory")
    directory_path = input().strip()
else:
    if len(sys.argv) > 1:
        directory_path = sys.argv[1]
    else:
        print("Erro: Por favor, forneÃ§a o caminho completo do diretÃ³rio como argumento na linha de comando.")
        sys.exit()
    # Salvando a tabela consolidada de uso de códons

extract_cds_with_prodigal(directory_path)
consolidated_codon_table = extract_and_save_codon_usage(directory_path) # Garanta que esta função retorne consolidated_codon_table
consolidated_csv_file_path = os.path.join(directory_path, "consolidated_codon_usage.csv")
generate_codon_usage_csv(consolidated_codon_table, codon_to_aa, consolidated_csv_file_path)
print(f"Consolidated codon usage CSV saved to {consolidated_csv_file_path}")
# Removido o analyze_codon_usage(directory_path) redundante
codon_counts = analyze_codon_usage(directory_path) # Garanta que esta função retorne algo útil para plot_top_20_codons
plot_top_20_codons(codon_counts)
analyze_genomes_for_hgt(directory_path)
# Correção: Usando o nome correto da variável.
print("Consolidated Codon Usage:")
for codon, usage in consolidated_codon_table.items():
    print(f"{codon}: {usage}")

