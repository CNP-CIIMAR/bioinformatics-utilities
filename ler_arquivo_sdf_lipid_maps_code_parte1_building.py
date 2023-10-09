#!/usr/bin/env python
# coding: utf-8

# In[1]:


from rdkit import Chem
from rdkit.Chem.Draw import MolsToGridImage
import pandas as pd

from rdkit import Chem
import pandas as pd

input_file = 'C:/Users/matto/Desktop/find_fatty_acids_NP/structures.sdf'
output_file = 'C:/Users/matto/Desktop/find_fatty_acids_NP/structuresv2.sdf'

# Read the input SDF file
suppl = Chem.SDMolSupplier(input_file)

# Create a new SDF writer for the output file
writer = Chem.SDWriter(output_file)

# Loop over each molecule in the input file
for mol in suppl:
    # Check if the molecule is valid
    if mol is not None:
        try:
            Chem.SanitizeMol(mol)
            # Write the valid molecule to the output file
            writer.write(mol)
        except ValueError:
            # Ignore the molecule if it has an explicit valence error
            print(f"Ignoring molecule with explicit valence error: {mol.GetProp('_Name')}")
            continue

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw

# Carrega o arquivo SDF
suppl = Chem.SDMolSupplier('structures.sdf')

# Cria uma lista para armazenar as informações dos ácidos graxos
acidos_graxos = []

# Define o padrão SMARTS para ácidos graxos
padrao = Chem.MolFromSmarts('[#6](-[#6])[-;!R;!$([#6]=*)](-[#6])[-;!R;!$([#6]=*)](-[#6])[-;!R;!$([#6]=*)](-[#6])[-;!R;!$([#6]=*)](-[#6])[-;!R;!$([#6]=*)]')

# Itera sobre as moléculas contidas no arquivo SDF
for mol in suppl:
    if mol is not None:
        # Verifica se a molécula contém um ácido graxo
        if mol.HasSubstructMatch(padrao):
            # Extrai as informações do ácido graxo
            nome = mol.GetProp('_Name')
            formula_estrutural = Chem.MolToSmiles(mol)
            formula_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
            tamanho_cadeia = max([len(cadeia) for cadeia in Chem.rdmolops.GetMolFrags(mol) if all(atom.GetSymbol() == 'C' for atom in cadeia)]) + 1
            
            # Adiciona as informações do ácido graxo à lista
            acidos_graxos.append({'Nome': nome, 'Fórmula Estrutural': formula_estrutural, 'Fórmula SMILES': formula_smiles, 'Tamanho da Cadeia': tamanho_cadeia})

# Cria um DataFrame com as informações dos ácidos graxos
df = pd.DataFrame(acidos_graxos)

# Imprime o DataFrame na tela
print(df)

# Desenha as moléculas dos ácidos graxos e salva em um arquivo .png
for i, mol in enumerate(suppl):
    if mol is not None:
        # Verifica se a molécula contém um ácido graxo
        if mol.HasSubstructMatch(padrao):
            # Extrai o nome do ácido graxo
            nome = mol.GetProp('_Name')
            
            # Define o tamanho da imagem
            Draw.MolToFile(mol, f'{nome} (Tamanho da Cadeia = {df.iloc[i]["Tamanho da Cadeia"]}).png', size=(500, 500), legend=nome)
            
            # Imprime a molécula na tela
            Draw.MolToImage(mol, size=(500, 500)).show()


# In[ ]:




