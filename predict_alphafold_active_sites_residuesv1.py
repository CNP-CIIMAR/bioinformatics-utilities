import subprocess
import argparse
import json
import re
import tempfile
import sys
import os
import networkx as nx
import matplotlib.pyplot as plt
import os
import pandas as pd  # Importando o pandas para manipulaÃ§Ã£o de dataframes
import matplotlib.cm as cm
import requests
import pandas as pd
## Developer: Leandro de Mattos Pereira
# CNP Junior Researcher - CNP - team 
# Article authors: Leandro de Mattos Pereira, Anne Liong, Pedro LeÃ£o - Team Leader

required_packages = [
    "absl-py==2.0.0",
    "astunparse==1.6.3",
    "beautifulsoup4==4.12.2",
    "bio==1.5.9",
    "biopython==1.81",
    "biothings-client==0.3.0",
    "bs4==0.0.1",
    "cachetools==5.3.1",
    "certifi==2023.7.22",
    "charset-normalizer==3.2.0",
    "contourpy==1.1.1",
    "cycler==0.11.0",
    "flatbuffers==23.5.26",
    "fonttools==4.42.1",
    "gast==0.4.0",
    "google-auth==2.23.0",
    "google-auth-oauthlib==1.0.0",
    "google-pasta==0.2.0",
    "gprofiler-official==1.0.0",
    "grpcio==1.58.0",
    "h5py==3.9.0",
    "idna==3.4",
    "keras==2.13.1",
    "kiwisolver==1.4.5",
    "libclang==16.0.6",
    "Markdown==3.4.4",
    "MarkupSafe==2.1.3",
    "matplotlib==3.8.0",
    "mygene==3.2.2",
    "numpy==1.24.3",
    "oauthlib==3.2.2",
    "opt-einsum==3.3.0",
    "packaging==23.1",
    "pandas==2.1.1",
    "Pillow==10.0.1",
    "pip==22.2.2",
    "platformdirs==3.10.0",
    "pooch==1.7.0",
    "protobuf==4.24.3",
    "pyasn1==0.5.0",
    "pyasn1-modules==0.3.0",
    "pyparsing==3.1.1",
    "python-dateutil==2.8.2",
    "pytz==2023.3.post1",
    "requests==2.31.0",
    "requests-oauthlib==1.3.1",
    "rsa==4.9",
    "setuptools==63.2.0",
    "six==1.16.0",
    "soupsieve==2.5",
    "tensorboard==2.13.0",
    "tensorboard-data-server==0.7.1",
    "tensorflow==2.13.0",
    "tensorflow-estimator==2.13.0",
    "tensorflow-io-gcs-filesystem==0.34.0",
    "termcolor==2.3.0",
    "tqdm==4.66.1",
    "typing_extensions==4.5.0",
    "tzdata==2023.3",
    "urllib3==1.26.16",
    "Werkzeug==2.3.7",
    "wheel==0.41.2",
    "wrapt==1.15.0"
]

# Verifique se o ambiente Python 'alpha_folder_env' existe
env_name = 'alpha_folder_env'
env_path = os.path.expanduser(f'~/.virtualenvs/{env_name}')
if not os.path.exists(env_path):
    print(f"O ambiente Python '{env_name}' nÃ£o foi encontrado.")
    create_env = input("Gostaria de criar o ambiente e instalar os pacotes? (y/n): ")
    if create_env.lower() == 'y':
        # Crie o ambiente Python e instale os pacotes
        subprocess.run(f"virtualenv {env_path}", shell=True)
        activate_script = os.path.join(env_path, 'bin', 'activate')
        subprocess.run(f"source {activate_script}", shell=True)
        for package in required_packages:
            subprocess.run(f"pip install {package}", shell=True)
        print(f"Ambiente '{env_name}' criado e pacotes instalados.")
    else:
        print("Continuando sem criar o ambiente. Certifique-se de que o ambiente 'alpha_folder_env' esteja ativado.")
else:
    # O ambiente jÃ¡ existe, entÃ£o apenas o ative
    activate_script = os.path.join(env_path, 'bin', 'activate')
    subprocess.run(f"source {activate_script}", shell=True)



# Verifique se o PyMOL existe no diretÃ³rio especificado
pymol_path = '/home/mattoslmp/pymol-open-source/bin/pymol'
if not os.path.exists(pymol_path):
    print("PyMOL nÃ£o foi encontrado no diretÃ³rio especificado.")
    install_pymol = input("Gostaria de instalar o PyMOL? (y/n): ")
    if install_pymol.lower() == 'y':
        # Instale o PyMOL aqui (ajuste o comando conforme suas necessidades)
        subprocess.run("conda install -c schrodinger pymol", shell=True)
    else:
        sys.exit("Encerrando o programa.")
            
            
# DicionÃ¡rio com resÃ­duos de interesse em FAAL32
FAAL32_active_site_residues = {
    "adenosine_binding": {
        "residues": ["Pro-315", "Tyr-342", "Ile-479", "Ser-313", "Glu-314", "Pro-315", "Val-316"],
        "interactions": ["hydrophobic_interactions", "backbone_atoms_interactions", "hydrogen_bonds"]
    },
    "ribose_phosphate_binding": {
        "residues": ["Asp-468", "Arg-482", "Lys-600", "Asp-231", "His-230"],
        "interactions": ["polar_interactions"]
    },
    "hydrophobic_tunnel_for_fatty_acid": {
        "residues": ["Val-210", "Leu-214", "Met-232", "Ile-235", "Thr-236", "Leu-239", "Phe-247", "Phe-277", "Ser-278",
                     "Ala-279", "Leu-310", "Asn-311", "Gly-312", "Ser-313", "Ser-341", "Tyr-342", "Gly-343", "Leu-344",
                     "Ala-345", "Leu-349", "Phe-350"],
        "interactions": ["hydrophobic_interactions"]
    }
}

PYMOL_PATH = os.getenv('PYMOL_PATH', '/home/mattoslmp/pymol-open-source/bin/pymol')
COLORS = ['red', 'blue', 'green', 'yellow', 'magenta', 'cyan', 'orange', 'salmon', 'lime', 'pink']

TMALIGN_PATH = "/home/mattoslmp/alphafold_models/TMalign"

def align_and_get_tmscore(target, reference):
    cmd = f"{TMALIGN_PATH} {target} {reference} > tmp_output.txt"  # Atualizado para usar a variÃ¡vel TMALIGN_PATH
    os.system(cmd)
    
    tmscore = None
    with open("tmp_output.txt", "r") as f:
        for line in f:
            if "TM-score=" in line:
                tmscore = float(line.split()[1])
                break
    os.remove("tmp_output.txt")
    return tmscore

def create_protein_graph(tmscore_df, threshold=0.5):
    G = nx.Graph()
    
    for protein1, tmscore_values in tmscore_df.items():
        for protein2, tmscore in tmscore_values.items():
            if tmscore > threshold:
                G.add_edge(protein1, protein2, weight=tmscore)
    
    return G


def plot_graph(G):
    # ConfiguraÃ§Ã£o da figura
    fig, ax = plt.subplots(figsize=(12, 12))
    
    # Layout
    pos = nx.spring_layout(G, seed=42, iterations=50)
    
    # Cores e tamanhos para os nÃ³s
    node_colors = [G.degree(v) for v in G]
    node_sizes = [G.degree(v) * 200 for v in G]  # Aumentar o tamanho dos nÃ³s
    
    # Desenhar nÃ³s
    nodes = nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors, cmap=plt.cm.viridis, ax=ax)
    
    # Cores para as arestas
    edge_weights = [G[u][v]['weight'] for u, v in G.edges()]
    edge_colors = [cm.viridis(w) for w in edge_weights]
    
    # Desenhar arestas
    edges = nx.draw_networkx_edges(G, pos, ax=ax, edge_color=edge_colors, width=2)
    
    # RÃ³tulos de nÃ³
    labels = {}
    for node in G.nodes():
        labels[node] = node
    nx.draw_networkx_labels(G, pos, labels, font_size=16, ax=ax)  # Aumentar o tamanho da fonte
    
    # RÃ³tulos de aresta
    edge_labels = nx.get_edge_attributes(G, 'weight')
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=12, ax=ax)  # Aumentar o tamanho da fonte
    
    # ConfiguraÃ§Ã£o das colorbars
    sm_nodes = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=plt.Normalize(min(node_colors), max(node_colors)))
    sm_nodes._A = []
    cbar_nodes = plt.colorbar(sm_nodes, ax=ax, orientation="horizontal", label="Node Degree", shrink=0.5)
    
    sm_edges = plt.cm.ScalarMappable(cmap=cm.viridis, norm=plt.Normalize(min(edge_weights), max(edge_weights)))
    sm_edges._A = []
    cbar_edges = plt.colorbar(sm_edges, ax=ax, orientation="horizontal", label="Edge Weight", shrink=0.5)
    
    # FinalizaÃ§Ã£o
    plt.title('Protein Similarity Graph', fontsize=20)  # Adicionar tÃ­tulo
    plt.axis("off")
    plt.show()


# Supondo que G Ã© o seu grafo
# plot_graph(G)

def main_extended(csv_file_path):
    df = pd.read_csv(csv_file_path)  # Carregando o CSV para um dataframe
    
    tmscore_df = {}
    protein_ids = df['Protein.accession'].tolist()
    
    for i, protein1 in enumerate(protein_ids):
        for j, protein2 in enumerate(protein_ids):
            if i >= j:
                continue
            
            tmscore_value = align_and_get_tmscore(f"alphafold_models/{protein1}.pdb", f"alphafold_models/{protein2}.pdb")
            
            if tmscore_value is not None:
                if protein1 not in tmscore_df:
                    tmscore_df[protein1] = {}
                tmscore_df[protein1][protein2] = tmscore_value
    
    # Create the graph
    G = create_protein_graph(tmscore_df)
    
    # Plot the graph
    plot_graph(G)
    
    # Save the TM-score DataFrame
    df_tmscore = pd.DataFrame.from_dict(tmscore_df, orient='index')
    df_tmscore.to_csv('all_vs_all_tmscore_values.csv')

# Chamada da funÃ§Ã£o main_extended
main_extended('/home/mattoslmp/alphafold_models/teste.id.csv')


def run_subprocess(cmd):
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        print(f"Command failed with error code {process.returncode}: {stderr.decode('utf-8')}")
        return None
    return stdout.decode('utf-8')

def download_file(url, filename):
    cmd = f"curl -o {filename} {url}"
    run_subprocess(cmd)

def fetch_alphafold_data(protein_id):
    cmd = f"curl -X 'GET' 'https://alphafold.ebi.ac.uk/api/prediction/{protein_id}' -H 'accept: application/json'"
    response_data = run_subprocess(cmd)

    try:
        data = json.loads(response_data)[0]
        os.makedirs("alphafold_models", exist_ok=True)
        with open(f"alphafold_models/{protein_id}.json", 'w') as json_file:
            json.dump(data, json_file)
        download_file(data["cifUrl"], f"alphafold_models/{protein_id}.cif")
        download_file(data["pdbUrl"], f"alphafold_models/{protein_id}.pdb")
    except json.JSONDecodeError:
        print(f"Failed to fetch data for {protein_id}. Response: {response_data}")

def align_and_color_by_rmsd(target, reference, aligned_filename, alignment_image, color_index=0):
    ref_helix_color = COLORS[color_index % len(COLORS)]
    ref_sheet_color = COLORS[(color_index + 1) % len(COLORS)]
    target_helix_color = COLORS[(color_index + 2) % len(COLORS)]
    target_sheet_color = COLORS[(color_index + 3) % len(COLORS)]

    cmd = f"{PYMOL_PATH} -c -d 'run colorbyrmsd.py; load {reference}, ref; load {target}, target; "
    cmd += f"super target and name CA, ref and name CA; "
    cmd += f"colorbyrmsd target, ref; "
    cmd += f"color {ref_helix_color}, ref and ss h; "
    cmd += f"color {target_helix_color}, target and ss h; "
    cmd += f"color {ref_sheet_color}, ref and ss s; "
    cmd += f"color {target_sheet_color}, target and ss s; "
    cmd += f"ray; png {alignment_image}, dpi=300; save {aligned_filename}, target;'"

    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()

    if process.returncode != 0:
        print(f"PyMOL command failed with error code {process.returncode}: {stderr.decode('utf-8')}")
        return None

    lines = stdout.decode('utf-8').split('\n')
    for line in lines:
        if "RMSD" in line:
            match = re.search(r'RMSD\s*=\s*([\d\.]+)', line)
            if match:
                rmsd_value = float(match.group(1))
                return rmsd_value


def render_structure(pdb_file, output_file, residues_to_highlight, color_index=0, zoomed_output_file=None):
    COLORS = ['red', 'blue', 'green', 'yellow', 'magenta', 'cyan', 'orange', 'salmon', 'lime', 'pink']
    helix_color = COLORS[color_index % len(COLORS)]
    sheet_color = COLORS[(color_index + 1) % len(COLORS)]
    
    residue_string = " or ".join([f"resi {residue.split('-')[1]} and resn {residue.split('-')[0][0:3]}" for residue in residues_to_highlight])
    
    with tempfile.NamedTemporaryFile(mode="w", suffix=".pml", delete=False) as temp_script:
        temp_script.write(f"load {pdb_file}, protein;\n")
        temp_script.write("hide everything;\n")
        temp_script.write("show cartoon;\n")
        temp_script.write(f"color {helix_color}, ss h;\n")
        temp_script.write(f"color {sheet_color}, ss s;\n")
        temp_script.write(f"select activesite, ({residue_string}) and name CA;\n")
        temp_script.write("show sticks, activesite;\n")
        
        # Rotulamento dos resÃ­duos. Utiliza apenas o nome do resÃ­duo e o nÃºmero do resÃ­duo.
     #   temp_script.write('label activesite, (resn, resi)\n')~
        temp_script.write('label activesite, "resn: " + resn + " resi: " + resi\n')

        
        temp_script.write("ray;\n")
        temp_script.write(f"png {output_file}, dpi=300;\n")
        
        if zoomed_output_file:
            temp_script.write("zoom activesite, 5;\n")
            temp_script.write(f"png {zoomed_output_file}, dpi=300;\n")
            
    ret = os.system(f"{PYMOL_PATH} -c -r {temp_script.name}")
    os.remove(temp_script.name)
    
    if ret != 0:
        print(f"PyMOL command failed with error code {ret}")

# Função para verificar se um ID é um ID do UniProtKB
def is_uniprot_id(protein_id):
    return re.match(r'^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$', protein_id) is not None

# Suponha que df seja o seu DataFrame e que ele tenha uma coluna 'Protein.accession'
# df = pd.read_csv("seu_arquivo.csv")

# Crie uma nova coluna 'Protein.accession_Nr-database' que é uma cópia de 'Protein.accession'
df['Protein.accession_Nr-database'] = df['Protein.accession']

# Obtenha uma lista única de IDs em 'Protein.accession' que NÃO são IDs do UniProtKB
unique_ids = [pid for pid in df['Protein.accession'].unique() if not is_uniprot_id(pid)]

def map_ids(id_list, from_db='P_REFSEQ_AC', to_db='ID'):
    url = 'https://www.uniprot.org/uploadlists/'
    params = {
        'from': from_db,
        'to': to_db,
        'format': 'tab',
        'query': ' '.join(id_list)
    }

    response = requests.get(url, params=params)
    if response.ok:
        lines = response.text.split('\n')
        mapped_ids = {}
        for line in lines[1:]:  # Skip the header
            if not line:
                continue
            src_id, target_id = line.split('\t')
            mapped_ids[src_id] = target_id
        return mapped_ids
    else:
        response.raise_for_status()

# Mapeie esses IDs para seus correspondentes em UniProtKB
mapped_ids = map_ids(unique_ids)

# Atualize a coluna 'Protein.accession' com os IDs mapeados
def update_id(protein_id):
    if is_uniprot_id(protein_id):
        return protein_id
    return mapped_ids.get(protein_id, protein_id)

df['Protein.accession'] = df['Protein.accession'].apply(update_id)

def main(csv_file_path):
    from pandas import read_csv  # Importing here to avoid unused import if the function is not used
    
    df = pd.read_csv(csv_file_path)

    # Cria uma nova coluna que é uma cópia da coluna 'Protein.accession'
    df['Protein.accession_Nr-database'] = df['Protein.accession']

    # Obter uma lista única de IDs da coluna 'Protein.accession'
    unique_ids = df['Protein.accession'].unique().tolist()

    # Use a função map_ids para mapear esses IDs para seus correspondentes no UniProtKB
    mapped_ids = map_ids(unique_ids)

    # Atualiza a coluna 'Protein.accession' com os IDs mapeados
    df['Protein.accession'] = df['Protein.accession'].apply(lambda x: mapped_ids.get(x, x))
             
    reference_protein = "O53580"
    fetch_alphafold_data(reference_protein)

    all_residues = []
    for site, info in FAAL32_active_site_residues.items():
        all_residues.extend(info.get("residues", []))

    rmsd_df = {}

    # Render the reference protein
    render_structure(f"alphafold_models/{reference_protein}.pdb", f"alphafold_models/{reference_protein}.png", all_residues, zoomed_output_file=f"alphafold_models/{reference_protein}_zoom.png")

    for protein_id in df['Protein.accession']:
        if protein_id == reference_protein:
            continue

        fetch_alphafold_data(protein_id)

        rmsd_value = align_and_color_by_rmsd(f"alphafold_models/{protein_id}.pdb", f"alphafold_models/{reference_protein}.pdb",
                                             f"alphafold_models/{protein_id}_aligned.pdb", f"alphafold_models/{protein_id}_alignment.png")

        if rmsd_value is not None:
            rmsd_df[protein_id] = {reference_protein: rmsd_value}

        # Render the aligned protein with and without zoom
        render_structure(f"alphafold_models/{protein_id}_aligned.pdb", f"alphafold_models/{protein_id}.png", all_residues, zoomed_output_file=f"alphafold_models/{protein_id}_zoom.png")

    df_rmsd = pd.DataFrame.from_dict(rmsd_df, orient='index')
    df_rmsd.to_csv('rmsd_values.csv')
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch AlphaFold data for protein IDs in a CSV file.")
    parser.add_argument('csv_file_path', type=str, help="Path to the CSV file containing a 'Protein.accession' column with protein IDs.")
    args = parser.parse_args()
    main(args.csv_file_path)
