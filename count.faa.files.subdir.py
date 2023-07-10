#################### Count arquivos .faa em um subdirectory
import os

# Pasta onde estão os subdiretórios
pasta_principal = './data'

# Variável para armazenar a contagem de arquivos .faa
contagem_arquivos_faa = 0

# Percorrendo a árvore de diretórios
for diretorio_raiz, subdiretorios, arquivos in os.walk(pasta_principal):
    for subdiretorio in subdiretorios:
        diretorio_completo = os.path.join(diretorio_raiz, subdiretorio)

        # Procurando arquivos .faa
        for arquivo in os.listdir(diretorio_completo):
            if arquivo.endswith('.faa'):
                contagem_arquivos_faa += 1

print(f"Total de arquivos .faa encontrados: {contagem_arquivos_faa}")
