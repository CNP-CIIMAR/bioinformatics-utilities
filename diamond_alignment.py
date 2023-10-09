#Autor: Leandro de mattos Pereira
import subprocess
import csv
from Bio import SeqIO
import re
import argparse

def verificar_ambiente_existente():
    return os.path.exists("diamond_alignment_env") and os.path.exists("uniprot_sprot.fasta")

def instalar_pacotes():
    try:
        if not verificar_ambiente_existente():
            # Criar ambiente virtual
            subprocess.run(["python3", "-m", "venv", "diamond_alignment_env"])
            
            # Ativar ambiente virtual e instalar pacotes
            # Note que a ativaÃ§Ã£o do ambiente virtual dentro do script Ã© um desafio.
            # VocÃª pode preferir ativar o ambiente manualmente antes de executar o script.
            subprocess.run(["source", "diamond_alignment_env/bin/activate"])
            subprocess.run(["pip", "install", "biopython"])
    except Exception as e:
        print(f"Erro durante a instalaÃ§Ã£o dos pacotes: {e}")

def instalar_diamond():
    try:
        if not os.path.exists("/usr/bin/diamond"):
            subprocess.run(["sudo", "apt-get", "update"])
            subprocess.run(["sudo", "apt-get", "install", "-y", "diamond-aligner"])
    except Exception as e:
        print(f"Erro durante a instalaÃ§Ã£o do DIAMOND: {e}")

def baixar_uniprotkb():
    try:
        if not os.path.exists("uniprot_sprot.fasta"):
            subprocess.run(["wget", "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"])
            subprocess.run(["gunzip", "uniprot_sprot.fasta.gz"])
    except Exception as e:
        print(f"Erro durante o download do UniProtKB: {e}")
        
def criar_banco_diamond():
    try:
        subprocess.run(["diamond", "makedb", "--in", "uniprot_sprot.fasta", "-d", "uniprotkb"])
    except Exception as e:
        print(f"Erro durante a criaÃ§Ã£o do banco de dados do DIAMOND: {e}")

def executar_diamond(meu_multifasta):
    subprocess.run(["diamond", "blastp", "-d", "uniprotkb.dmnd", "-q", meu_multifasta, "-o", "resultados.tab", "-f", "6", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"])

def ler_resultados(arquivo_resultados):
    resultados = []
    with open(arquivo_resultados, "r") as f:
        leitor = csv.reader(f, delimiter='\t')
        for linha in leitor:
            resultados.append({
                "qseqid": linha[0],
                "sseqid": linha[1],
                "pident": float(linha[2]),
                "length": int(linha[3]),
                "qstart": int(linha[6]),
                "qend": int(linha[7]),
                "sstart": int(linha[8]),
                "send": int(linha[9]),
            })
    return resultados

def extrair_especie(descricao):
    especie_match = re.search(r'\[([^\]]+)\]', descricao)
    if especie_match:
        return especie_match.group(1)
    return None

def filtrar_resultados(resultados, meu_multifasta, uniprot_fasta):
    resultados_filtrados = []
    query_seqs = SeqIO.to_dict(SeqIO.parse(meu_multifasta, "fasta"))
    uniprot_seqs = SeqIO.to_dict(SeqIO.parse(uniprot_fasta, "fasta"))

    for r in resultados:
        if r["pident"] == 100.0:
            query_seq = query_seqs.get(r["qseqid"])
            uniprot_seq = uniprot_seqs.get(r["sseqid"])
            if query_seq and uniprot_seq:
                query_length = len(query_seq)
                if r["length"] == query_length:
                    query_especie = extrair_especie(query_seq.description)
                    uniprot_especie = extrair_especie(uniprot_seq.description)
                    if query_especie == uniprot_especie:
                        r["query_especie"] = query_especie
                        r["uniprot_especie"] = uniprot_especie
                        r["sequencia_homologa"] = str(uniprot_seq.seq)  # Adiciona a sequÃªncia homÃ³loga
                        resultados_filtrados.append(r)
    return resultados_filtrados


def gerar_multifasta(resultados_filtrados, arquivo_saida, uniprot_fasta):
    sequencias = SeqIO.to_dict(SeqIO.parse(uniprot_fasta, "fasta"))
    with open(arquivo_saida, "w") as f:
        for resultado in resultados_filtrados:
            sseqid = resultado["sseqid"]
            if sseqid in sequencias:
                SeqIO.write(sequencias[sseqid], f, "fasta")

def gerar_tabela(resultados_filtrados, arquivo_saida):
    with open(arquivo_saida, "w", newline='') as f:
        escritor = csv.writer(f)
        escritor.writerow(["qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart", "send", "query_especie", "uniprot_especie", "sequencia_homologa"])
        for resultado in resultados_filtrados:
            escritor.writerow([resultado["qseqid"], resultado["sseqid"], resultado["pident"],
                               resultado["length"], resultado["qstart"], resultado["qend"],
                               resultado["sstart"], resultado["send"],
                               resultado["query_especie"], resultado["uniprot_especie"],
                               resultado["sequencia_homologa"]])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Executa o alinhamento de sequÃªncias usando DIAMOND e filtra os resultados.')
    parser.add_argument('--input', required=True, help='Arquivo multifasta de entrada.')
    args = parser.parse_args()

    meu_multifasta = args.input

    instalar_pacotes()
    instalar_diamond()
    baixar_uniprotkb()
    criar_banco_diamond()
    executar_diamond(meu_multifasta)
    resultados = ler_resultados("resultados.tab")
    resultados_filtrados = filtrar_resultados(resultados, meu_multifasta, "uniprot_sprot.fasta")
    gerar_multifasta(resultados_filtrados, "homologos_uniprotkb.fasta", "uniprot_sprot.fasta")
    gerar_tabela(resultados_filtrados, "tabela_alinhamento.csv")
