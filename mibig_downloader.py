import os
import re
import argparse
import requests
from urllib.parse import urljoin
from time import sleep
import zipfile
import shutil

def extrair_ids_mibig(arquivo_txt):
    """
    Extrai os IDs que seguem o padrão MIBIG (BGC seguido de 7 dígitos) do arquivo de texto.
    """
    padrão = re.compile(r'\bBGC\d{7}\b')
    ids = set()
    with open(arquivo_txt, 'r', encoding='utf-8') as f:
        for linha in f:
            correspondências = padrão.findall(linha)
            ids.update(correspondências)
    return list(ids)

def baixar_arquivo(url, caminho_destino, tentativas=3, delay=5):
    """
    Baixa um arquivo a partir de uma URL e salva no caminho de destino.
    Tenta novamente em caso de falha, até o número de tentativas especificado.
    """
    for tentativa in range(1, tentativas + 1):
        try:
            resposta = requests.get(url, stream=True, timeout=10)
            resposta.raise_for_status()
            with open(caminho_destino, 'wb') as f:
                for chunk in resposta.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
            print(f"Baixado: {caminho_destino}")
            return True
        except requests.HTTPError as http_err:
            if resposta.status_code == 404:
                print(f"Arquivo não encontrado (404): {url}")
                return False
            else:
                print(f"Erro HTTP ao baixar {url}: {http_err} (Tentativa {tentativa}/{tentativas})")
        except requests.ConnectionError as conn_err:
            print(f"Erro de conexão ao baixar {url}: {conn_err} (Tentativa {tentativa}/{tentativas})")
        except requests.Timeout as timeout_err:
            print(f"Timeout ao baixar {url}: {timeout_err} (Tentativa {tentativa}/{tentativas})")
        except Exception as err:
            print(f"Erro ao baixar {url}: {err} (Tentativa {tentativa}/{tentativas})")
        
        if tentativa < tentativas:
            print(f"Aguardando {delay} segundos antes da próxima tentativa...")
            sleep(delay)
    
    print(f"Falha ao baixar após {tentativas} tentativas: {url}")
    return False

def descompactar_arquivo(caminho_zip, diretorio_destino):
    """
    Descompacta um arquivo .zip para o diretório de destino especificado.
    """
    try:
        with zipfile.ZipFile(caminho_zip, 'r') as zip_ref:
            # Cria um subdiretório para cada arquivo zip
            nome_pasta = os.path.splitext(os.path.basename(caminho_zip))[0]
            pasta_destino = os.path.join(diretorio_destino, nome_pasta)
            os.makedirs(pasta_destino, exist_ok=True)
            zip_ref.extractall(pasta_destino)
        print(f"Descompactado: {caminho_zip} para {pasta_destino}")
        return pasta_destino
    except zipfile.BadZipFile:
        print(f"Erro: O arquivo {caminho_zip} não é um arquivo zip válido.")
    except Exception as err:
        print(f"Erro ao descompactar {caminho_zip}: {err}")
    return None

def mover_gbk_para_diretorio(diretorio_origem, diretorio_destino):
    """
    Move todos os arquivos .gbk do diretório de origem para o diretório de destino.
    """
    for item in os.listdir(diretorio_origem):
        if item.endswith('.gbk'):
            origem = os.path.join(diretorio_origem, item)
            destino = os.path.join(diretorio_destino, item)
            shutil.move(origem, destino)
            print(f"Mover: {origem} para {destino}")

def main():
    parser = argparse.ArgumentParser(description='Baixa, descompacta e organiza arquivos MIBIG a partir de um arquivo de IDs.')
    parser.add_argument('arquivo_entrada', help='Caminho para o arquivo de entrada (txt) contendo os IDs.')
    parser.add_argument('diretorio_saida', help='Diretório onde os arquivos .zip serão salvos e descompactados.')
    parser.add_argument('--url_base', default='https://mibig.secondarymetabolites.org/repository/', help='URL base do repositório MIBIG.')
    parser.add_argument('--log_falhas', default='falhas_download.txt', help='Arquivo para registrar IDs que falharam no download.')
    args = parser.parse_args()

    arquivo_entrada = args.arquivo_entrada
    diretorio_saida = args.diretorio_saida
    url_base = args.url_base
    log_falhas = args.log_falhas

    # Cria o diretório de saída se não existir
    os.makedirs(diretorio_saida, exist_ok=True)

    # Extrai os IDs MIBIG do arquivo de entrada
    ids_mibig = extrair_ids_mibig(arquivo_entrada)
    print(f"IDs encontrados: {len(ids_mibig)}")

    falhas = []

    for id_mibig in ids_mibig:
        url = urljoin(url_base, f"{id_mibig}/{id_mibig}.zip")
        caminho_destino = os.path.join(diretorio_saida, f"{id_mibig}.zip")
        
        # Verifica se o arquivo já existe para evitar downloads duplicados
        if os.path.exists(caminho_destino):
            print(f"Arquivo já existe, pulando download: {caminho_destino}")
        else:
            sucesso = baixar_arquivo(url, caminho_destino)
            if not sucesso:
                falhas.append(id_mibig)
                continue  # Pula para o próximo ID se o download falhar

        # Descompacta o arquivo .zip
        pasta_destino = descompactar_arquivo(caminho_destino, diretorio_saida)
        if not pasta_destino:
            falhas.append(id_mibig)
            continue  # Pula para o próximo ID se a descompactação falhar

        # Move os arquivos .gbk para o diretório principal
        mover_gbk_para_diretorio(pasta_destino, diretorio_saida)

    # Registra as falhas, se houver
    if falhas:
        with open(os.path.join(diretorio_saida, log_falhas), 'w', encoding='utf-8') as f:
            for id_falha in falhas:
                f.write(f"{id_falha}\n")
        print(f"\nProcesso concluído com falhas. Consulte '{log_falhas}' para ver os IDs que não puderam ser baixados ou descompactados.")
    else:
        print("\nTodos os arquivos foram baixados, descompactados e organizados com sucesso.")

if __name__ == "__main__":
    main()
