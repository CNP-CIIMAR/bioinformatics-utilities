import requests

def download_fasta(pdb_id):
    url = f'https://www.rcsb.org/fasta/entry/{pdb_id}'
    response = requests.get(url)
    if response.ok:
        fasta_data = response.text
        with open(f'{pdb_id}.fasta', 'w') as fasta_file:
            fasta_file.write(fasta_data)

# Obtém a lista de todos os IDs do PDB disponíveis com uma determinada anotação no cabeçalho
def get_pdb_ids(annotation):
    url = f'https://www.rcsb.org/pdb/rest/search'
    query = f'<orgPdbQuery><queryType>org.pdb.query.simple.AdvancedKeywordQuery</queryType><keywords>{annotation}</keywords></orgPdbQuery>'
    response = requests.post(url, data=query)
    if response.ok:
        pdb_ids = response.text.split()
        return pdb_ids

# Itera sobre todos os IDs do PDB com a anotação desejada e baixa os arquivos FASTA
annotation = 'sua_anotacao_aqui'
pdb_ids = get_pdb_ids(annotation)
for pdb_id in pdb_ids:
    download_fasta(pdb_id)
