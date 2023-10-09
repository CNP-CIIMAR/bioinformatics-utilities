import pandas as pd
import sys

# Verifica se o número correto de argumentos foi fornecido
if len(sys.argv) < 5:
    print("Por favor, forneça os nomes dos arquivos das tabelas fadD, fadB, fadA e fadE como argumentos.")
    sys.exit(1)

# Função para obter a lista de genomas com um gene específico
def obter_genomas(tabela, gene):
    genomas = tabela.loc[tabela.iloc[:, 1:].notnull().sum(axis=1) > 0, 'Genomas'].tolist()
    return genomas

# Lê os arquivos das tabelas
tabela_fadD = pd.read_csv(sys.argv[1], delimiter="\t")
tabela_fadB = pd.read_csv(sys.argv[2], delimiter="\t")
tabela_fadA = pd.read_csv(sys.argv[3], delimiter="\t")
tabela_fadE = pd.read_csv(sys.argv[4], delimiter="\t")

# Obtém os genomas com os genes fadD, fadB, fadA e fadE
genomas_fadD = set(obter_genomas(tabela_fadD, 'fadD'))
genomas_fadB = set(obter_genomas(tabela_fadB, 'fadB'))
genomas_fadA = set(obter_genomas(tabela_fadA, 'fadA'))
genomas_fadE = set(obter_genomas(tabela_fadE, 'fadE'))

# Calcula a interseção dos genomas com os genes fadD, fadB, fadA e fadE
genomas_intersecao_fadD_fadB_fadA_fadE = genomas_fadD.intersection(genomas_fadB, genomas_fadA, genomas_fadE)

# Imprime a quantidade de genomas na interseção de fadD, fadB, fadA e fadE
quantidade_intersecao_fadD_fadB_fadA_fadE = len(genomas_intersecao_fadD_fadB_fadA_fadE)
print("Quantidade de genomas com fadD, fadB, fadA e fadE:", quantidade_intersecao_fadD_fadB_fadA_fadE)
print("Genomas com fadD, fadB, fadA e fadE:")
for genoma in genomas_intersecao_fadD_fadB_fadA_fadE:
    print(genoma)

# Calcula a interseção dos genomas com os genes fadB, fadA e fadE
genomas_intersecao_fadB_fadA_fadE = genomas_fadB.intersection(genomas_fadA, genomas_fadE)

# Imprime a quantidade de genomas na interseção de fadB, fadA e fadE
quantidade_intersecao_fadB_fadA_fadE = len(genomas_intersecao_fadB_fadA_fadE)
print("Quantidade de genomas com fadB, fadA e fadE:", quantidade_intersecao_fadB_fadA_fadE)
print("Genomas com fadB, fadA e fadE:")
for genoma in genomas_intersecao_fadB_fadA_fadE:
    print(genoma)

# Calcula a interseção dos genomas com os genes fadD, fadB e fadA
genomas_intersecao_fadD_fadB_fadA = genomas_fadD.intersection(genomas_fadB, genomas_fadA)

# Imprime a quantidade de genomas na interseção de fadD, fadB e fadA
quantidade_intersecao_fadD_fadB_fadA = len(genomas_intersecao_fadD_fadB_fadA)
print("Quantidade de genomas com fadD, fadB e fadA:", quantidade_intersecao_fadD_fadB_fadA)
print("Genomas com fadD, fadB e fadA:")
for genoma in genomas_intersecao_fadD_fadB_fadA:
    print(genoma)
