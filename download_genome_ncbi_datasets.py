import subprocess
import sys
import logging
import os

def configurar_logging():
    logging.basicConfig(filename='processo_datasets.log', level=logging.INFO, 
                        format='%(asctime)s - %(levelname)s - %(message)s', 
                        datefmt='%Y-%m-%d %H:%M:%S')

def executar_comando_em_background(input_file, output_file):
    comando = [
        "datasets",
        "download",
        "genome",
        "accession",
        "--inputfile",
        input_file,
        "--include",
        "genome",
        "--filename",
        output_file
    ]

    # Configurar logging
    configurar_logging()

    # Abrir um arquivo para redirecionar a saída (similar ao nohup)
    with open('nohup_output.txt', 'w') as f:
        # Iniciar o processo em background, redirecionando saída e erro
        process = subprocess.Popen(comando, stdout=f, stderr=f)

    # Registrar a inicialização do processo
    logging.info(f"Processo iniciado em background com PID {process.pid}")

    return process.pid

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Uso: python script.py <arquivo_input_txt> <arquivo_output_zip>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    pid = executar_comando_em_background(input_file, output_file)
    print(f"Comando sendo executado em background com PID {pid}")

