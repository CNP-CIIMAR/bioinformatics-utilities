# 27 de august de 2024.
# Autor: Leandro de Mattos Pereira - bbf4 project
# Pedro Leão - team Leader, Cyanobacterial Natureal Products - team
import pandas as pd
import requests
import time
import argparse

# Função para obter o UID correspondente ao Genoma ID (GCF/GCA)
def get_uid_from_genome_id(genome_id):
    search_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=assembly&term={genome_id}&retmode=json"
    try:
        response = requests.get(search_url)
        if response.status_code == 200:
            data = response.json()
            if 'esearchresult' in data and 'idlist' in data['esearchresult'] and data['esearchresult']['idlist']:
                return data['esearchresult']['idlist'][0]  # Retorna o primeiro UID encontrado
            else:
                print(f"UID não encontrado para {genome_id}")
                return None
        elif response.status_code == 429:
            print(f"Erro 429: Excesso de requisições. Aguardando 60 segundos...")
            time.sleep(60)
            return get_uid_from_genome_id(genome_id)  # Tenta novamente após o delay
        else:
            print(f"Erro na requisição ao NCBI para {genome_id}: {response.status_code}")
            return None
    except Exception as e:
        print(f"Erro ao buscar UID para {genome_id}: {e}")
        return None

# Função para obter informações de coleta a partir do UID do NCBI
def get_metadata_from_ncbi(uid):
    ncbi_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=assembly&id={uid}&retmode=json"
    try:
        response = requests.get(ncbi_url)
        if response.status_code == 200:
            try:
                data = response.json()
            except requests.exceptions.JSONDecodeError:
                print(f"Erro ao decodificar JSON para UID {uid}")
                return None
            
            # Verifique se o UID está presente nos resultados
            if uid in data['result']:
                summary = data['result'][uid]
                
                # Extração de informações relevantes
                metadata = {
                    'Collection Date': summary.get('collection_date', None),
                    'Latitude': summary.get('latitude', None),
                    'Longitude': summary.get('longitude', None),
                    'Environment': summary.get('environment', None),
                }
                return metadata
            else:
                print(f"UID {uid} não encontrado na resposta da API.")
                return None
        elif response.status_code == 429:
            print(f"Erro 429: Excesso de requisições. Aguardando 60 segundos...")
            time.sleep(60)
            return get_metadata_from_ncbi(uid)  # Tenta novamente após o delay
        else:
            print(f"Erro na requisição à API do NCBI para UID {uid}: {response.status_code}")
            return None
    except Exception as e:
        print(f"Erro ao buscar metadados para UID {uid}: {e}")
        return None

def main(input_file, output_file):
    # Tentar carregar os dados a partir do arquivo TSV
    try:
        df = pd.read_csv(input_file, sep='\t', on_bad_lines='skip')  # Ignora linhas com problemas
    except Exception as e:
        print(f"Erro ao ler o arquivo {input_file}: {e}")
        return
    
    # Verifique os nomes das colunas carregadas
    print("Nomes das colunas no arquivo:", df.columns)

    # Se a coluna "Assembly" tiver outro nome, ajuste aqui
    if 'Assembly' not in df.columns:
        print("Coluna 'Assembly' não encontrada. Verifique o nome correto da coluna.")
        return
    
    # Adicionar novas colunas para as informações de coleta
    df['Collection Date'] = None
    df['Latitude'] = None
    df['Longitude'] = None
    df['Environment'] = None

    # Iterar sobre as linhas da tabela e buscar as informações de coleta
    for index, row in df.iterrows():
        genome_id = row['Assembly']
        
        # Verifique se o genome_id é uma string válida, não é NaN, e começa com "GCF" ou "GCA"
        if pd.isna(genome_id) or not isinstance(genome_id, str) or not genome_id.startswith(('GCF', 'GCA')):
            print(f"Pular linha {index}: Genoma ID inválido ({genome_id})")
            continue
        
        try:
            uid = get_uid_from_genome_id(genome_id)
            if uid:
                metadata = get_metadata_from_ncbi(uid)
                if metadata:
                    df.at[index, 'Collection Date'] = metadata['Collection Date']
                    df.at[index, 'Latitude'] = metadata['Latitude']
                    df.at[index, 'Longitude'] = metadata['Longitude']
                    df.at[index, 'Environment'] = metadata['Environment']
        except Exception as e:
            print(f"Erro ao processar a linha {index} com Assembly {genome_id}: {e}")

    # Salvar a tabela atualizada com as novas informações
    try:
        df.to_csv(output_file, sep='\t', index=False)
        print(f"Tabela atualizada salva em {output_file}")
    except Exception as e:
        print(f"Erro ao salvar o arquivo {output_file}: {e}")

if __name__ == "__main__":
    # Configuração de argumentos de linha de comando
    parser = argparse.ArgumentParser(description="Atualiza a tabela de genomas com informações adicionais de coleta.")
    parser.add_argument("input_file", help="Caminho para o arquivo TSV de entrada.")
    parser.add_argument("output_file", help="Caminho para o arquivo TSV de saída.")

    # Parse dos argumentos
    args = parser.parse_args()

    # Executar a função principal com os argumentos fornecidos
    main(args.input_file, args.output_file)
