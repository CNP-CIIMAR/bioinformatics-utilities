import sys

def filtrar_tabela(lista_arquivo, tabela_arquivo, arquivo_saida):
    # Ler os IDs da lista
    with open(lista_arquivo, 'r') as lista:
        ids = set(line.strip() for line in lista)

    # Armazenar linhas já escritas para cada ID
    linhas_por_id = {}

    # Ler a tabela e filtrar as linhas correspondentes aos IDs
    with open(tabela_arquivo, 'r') as tabela:
        # Iterar sobre as linhas da tabela
        for linha in tabela:
            # Verificar se a linha não está vazia
            if linha.strip():
                elementos = linha.split()
                # Verificar se a linha tem pelo menos um elemento
                if elementos:
                    # Extrair o ID da primeira coluna
                    id_na_tabela = elementos[0]
                    # Verificar se o ID está na lista
                    if id_na_tabela in ids:
                        # Adicionar a linha ao dicionário, sobrescrevendo se já existe uma linha para este ID
                        linhas_por_id[id_na_tabela] = linha.strip()

    # Escrever as linhas filtradas no arquivo de saída
    with open(arquivo_saida, 'w') as saida:
        # Escrever cabeçalho
        with open(tabela_arquivo, 'r') as tabela:
            saida.write(next(tabela))
        # Escrever linhas filtradas
        for linha in linhas_por_id.values():
            saida.write(linha + '\n')

if __name__ == "__main__":
    # Verificar se foram fornecidos os argumentos corretos
    if len(sys.argv) != 4:
        print("Uso: python script.py lista.txt tabela.txt tabela_filtrada.txt")
        sys.exit(1)

    # Obter nomes dos arquivos dos argumentos de linha de comando
    lista_arquivo = sys.argv[1]
    tabela_arquivo = sys.argv[2]
    arquivo_saida = sys.argv[3]

    # Chamar a função para filtrar a tabela
    filtrar_tabela(lista_arquivo, tabela_arquivo, arquivo_saida)
