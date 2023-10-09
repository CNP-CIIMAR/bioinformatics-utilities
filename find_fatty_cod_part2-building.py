#!/usr/bin/env python
# coding: utf-8

# In[1]:


pip install rdkit
pip install --upgrade rdkit


# In[4]:


#!/usr/bin/env python
import pikachu
import pandas as pd
from rdkit import Chem
import rdkit
from rdkit.Chem import Descriptors
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from PIL import Image
from IPython.display import display
from IPython.display import SVG
from rdkit.Chem import Draw
import os

os.chdir(r'/home/mattoslmp/Chemio_FAAL/')


# Ler a tabela de produtos naturais
df = pd.read_csv('/home/mattoslmp/Chemio_FAAL/NPAtlas_download.tsv',sep='\t')

# Concatenação das colunas genus e origin_species
df["genus_origin_species"] = df["genus"] + " " + df["origin_species"]


# Armazenar o nome e a estrutura química de cada produto natural
produtos = []
for i in range(len(df)):
    species = df.loc[i, 'genus_origin_species']
    nome = df.loc[i, 'compound_names']
    estrutura = Chem.MolFromSmiles(df.loc[i, 'compound_smiles'])
    molfile = Chem.MolToMolBlock(estrutura)
    formula = Chem.rdMolDescriptors.CalcMolFormula(estrutura)
    produtos.append((nome, estrutura, formula, molfile, species))
            
for produto in produtos:
    nome = produto[0]
    estrutura = produto[1]
    formula = produto[2]
    molfile = produto[3]
    species = produto[4]

# Definir lista de ácidos graxos
fatty_acids = [{'name': 'acético', 'smiles': 'CCCCC(=O)O'},{'name': 'propiônico', 'smiles': 'CCCCCCC(=O)O'},
               {'name': 'butírico', 'smiles': 'CCCCCCCCC(=O)O'},{'name': 'valérico', 'smiles': 'CCCCCCCCCCC(=O)O'},
               {'name': 'capróico', 'smiles': 'CCCCCCCCCCC(=O)O'}, {'name': 'enântico', 'smiles': 'CCCCCCCCCCCCC(=O)O'},
               {'name': 'caprílico', 'smiles': 'CCCCCCCCCCCCC(=O)O'}, {'name': 'pelargônico', 'smiles': 'CCCCCCCCCCCCCC(=O)O'},
               {'name': 'láurico', 'smiles': 'CCCCCCCCCCCCCCC(=O)O'}, {'name': 'mirístico', 'smiles': 'CCCCCCCCCCCCCCCC(=O)O'},
               {'name': 'palmítico', 'smiles': 'CCCCCCCCCCCCCCCCC(=O)O'},{'name': 'heptadecanóico', 'smiles': 'CCCCCCCCCCCCCCCCCC(=O)O'},
               {'name': 'esteárico', 'smiles': 'CCCCCCCCCCCCCCCCCCC(=O)O'}, {'name': 'araquídico', 'smiles': 'CCCCCCCCCCCCCCCCCCCC(=O)O'},
               {'name': 'behenico', 'smiles': 'CCCCCCCCCCCCCCCCCCCCC(=O)O'}, 
               {'name': 'lignocérico', 'smiles': 'CCCCCCCCCCCCCCCCCCCCCCC(=O)O'}, 
               {'name': 'cerótico', 'smiles': 'CCCCCCCCCCCCCCCCCCCCCCCC(=O)O'},
               {'name': 'montânico', 'smiles': 'CCCCCCCCCCCCCCCCCCCCCCCCC(=O)O'},
               {'name': 'melíssico', 'smiles': 'CCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O'}]


def search_substructures(fatty_acids, produtos, search_mode='single'):
    """
    Searches for substructures in the structure of natural compounds using RDKit.

    Parameters:
        fatty_acids (list): A list of SMILES string representations of the fatty acid substructures to search for.
        produtos (list): A list of tuples representing the natural products, where each tuple contains the
            product name, product structure (as a RDKit molecule object), product formula, product molfile,
            and the species of origin.
        search_mode (str): The search mode to use. Can be 'single' to find only the first occurrence,
            or 'all' to find all occurrences.

    Returns:
        None.
    """

    for produto in produtos:
        estrutura = produto[1]

        # Perform the substructure search.
        if search_mode == 'single':
            for fatty_acid in fatty_acids:
                fatty_acid_mol = Chem.MolFromSmiles(fatty_acid['smiles'])
                match = estrutura.GetSubstructMatch(fatty_acid_mol)
                if match:
                    Draw.MolToFile(estrutura, produto[0] + '_highlighted.png', highlightAtoms=match)
                    print(fatty_acid['name'], "is a substructure of", produto[0])
                    break
        elif search_mode == 'all':
            for fatty_acid in fatty_acids:
                fatty_acid_mol = Chem.MolFromSmiles(fatty_acid['smiles'])
                matches = estrutura.GetSubstructMatches(fatty_acid_mol)
                if matches:
                    for match in matches:
                        Draw.MolToFile(estrutura, produto[0] + '_highlighted.png', highlightAtoms=match)
                    print(fatty_acid['name'], "is a substructure of", produto[0])

def count_carbons(fatty_acid_mol):
    mol = Chem.MolFromSmiles(smiles)
    chains = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            neighbors = atom.GetNeighbors()
            if len(neighbors) == 2 and all(n.GetSymbol() == 'C' for n in neighbors):
                chain = Chem.PathToSubmol(mol, mol.GetSubstructMatch(Chem.MolFromSmarts('*-C(-[*])=O')))
                chains.append(chain)
    carbon_counts = {}
    for chain in chains:
        carbon_count = chain.GetNumAtoms() - 1
        carbon_counts[carbon_count] = carbon_counts.get(carbon_count, 0) + 1
    return carbon_counts


# In[5]:


print (carbon_counts)


# In[ ]:




