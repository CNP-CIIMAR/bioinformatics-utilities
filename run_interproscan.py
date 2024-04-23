import subprocess
import sys
import os

## Author: Leandro de Mattos Pereira
## CNP Laboratory - 23/04/2024
## Dr. Pedro Leão, Team Leader.
## The code automatize the identification of FAAL proteins using the signature results obtained from interproscan

def setup_environment():
    print("Criando ambiente virtual...")
    subprocess.run(['python3', '-m', 'venv', 'interproscan'], check=True)
    subprocess.run(['source', 'interproscan/bin/activate'], shell=True, executable='/bin/bash')

    print("Instalando Java JDK, Perl e dependências Python...")
    subprocess.run(['sudo', 'apt', 'update'], check=True)
    subprocess.run(['sudo', 'apt', 'install', '-y', 'default-jdk', 'perl', 'python3-pip'], check=True)
    
    print("Instalando bedtools...")
    subprocess.run(['sudo', 'apt', 'install', '-y', 'bedtools'], check=True)

    print("Verificando a existência do InterProScan...")
    if os.path.exists('./interproscan.sh'):
        print("InterProScan já está instalado.")
    else:
        print("InterProScan não encontrado. Baixando e descompactando...")
        subprocess.run(['wget', 'http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.67-99.0/interproscan-5.67-99.0-64-bit.tar.gz'], check=True)
        subprocess.run(['tar', '-xzf', 'interproscan-5.67-99.0-64-bit.tar.gz'], check=True)
        print("InterProScan baixado e descompactado com sucesso.")

def main():
    if len(sys.argv) != 2:
        print("Uso: python script.py <caminho_completo_do_arquivo.fasta>")
        sys.exit(1)

    input_fasta = sys.argv[1]
    if not input_fasta.endswith('.fasta'):
        print("O arquivo deve ter a extensão .fasta")
        sys.exit(1)

    setup_environment()  # Configura o ambiente antes de executar o script

    prefixo = os.path.basename(input_fasta)[:-6]

    try:
        subprocess.run(['./interproscan.sh', '-i', input_fasta, '-f', 'tsv', '-dp', '-appl', 'CDD'], check=True)
    except subprocess.CalledProcessError:
        print("Falha ao executar interproscan.")
        sys.exit(1)

    tsv_output = f"{prefixo}.fasta.tsv"

    if not os.path.exists(tsv_output):
        print(f"Erro: arquivo {tsv_output} não encontrado.")
        sys.exit(1)

    try:
        with open(f"{prefixo}.FAAL.tsv", 'w') as out_file:
            subprocess.run(['grep', 'FAAL', tsv_output], stdout=out_file, check=True)
    except subprocess.CalledProcessError:
        print("Falha ao filtrar linhas com 'FAAL'.")
        sys.exit(1)

    try:
        with open(f"{prefixo}.FAAL.bed", 'w') as out_file:
            subprocess.run(['cut', '-f1,7,8', f"{prefixo}.FAAL.tsv"], stdout=out_file, check=True)
    except subprocess.CalledProcessError:
        print("Falha ao usar cut para extrair colunas.")
        sys.exit(1)

    try:
        subprocess.run(['bedtools', 'getfasta', '-fi', input_fasta, '-bed', f"{prefixo}.FAAL.bed", '-fo', f"{prefixo}_output.fasta"], check=True)
    except subprocess.CalledProcessError:
        print("Falha ao extrair sequências fasta usando bedtools.")
        sys.exit(1)

    print("Operações concluídas com sucesso.")

if __name__ == "__main__":
    main()

