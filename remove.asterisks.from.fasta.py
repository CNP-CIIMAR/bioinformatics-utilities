import sys

# Verificar se o número correto de argumentos foi fornecido
if len(sys.argv) != 3:
    print("Uso: python seu_script.py arquivo_entrada.fasta arquivo_saida.fasta")
    sys.exit(1)

# Obter os caminhos dos arquivos de entrada e saída dos argumentos
entrada_path = sys.argv[1]
saida_path = sys.argv[2]

# Abrir o arquivo de entrada em modo de leitura
with open(entrada_path, 'r') as entrada:
    # Abrir o arquivo de saída em modo de escrita
    with open(saida_path, 'w') as saida:
        # Variáveis para armazenar o cabeçalho e a sequência atual
        cabecalho = ''
        sequencia = ''

        # Percorrer as linhas do arquivo de entrada
        for linha in entrada:
            linha = linha.strip()

            # Verificar se a linha é um cabeçalho
            if linha.startswith('>'):
                # Se for um novo cabeçalho, escrever a sequência anterior (se houver) no arquivo de saída
                if cabecalho:
                    saida.write(cabecalho + '\n')
                    saida.write(sequencia.replace('*', '') + '\n')

                # Armazenar o novo cabeçalho e redefinir a sequência
                cabecalho = linha
                sequencia = ''
            else:
                # Adicionar a linha à sequência
                sequencia += linha

        # Escrever o último cabeçalho e sequência no arquivo de saída
        if cabecalho:
            saida.write(cabecalho + '\n')
            saida.write(sequencia.replace('*', '') + '\n')
