import pandas as pd
import time
from Bio import Entrez
import xml.etree.ElementTree as ElementTree
from ipyleaflet import Map, CircleMarker, Popup, basemaps, basemap_to_tiles
from ipywidgets import HTML, Output, VBox, HBox
from IPython.display import display
import plotly.express as px
from plotly.io import to_html

# Configuração do email para NCBI Entrez (substitua pelo seu email real)
Entrez.email = 'seu.email@dominio.com'

# Função para buscar metadados do NCBI BioSample
def fetch_biosample_metadata(assembly_accession):
    try:
        search_handle = Entrez.esearch(db='assembly', term=f'{assembly_accession}[Assembly Accession]')
        search_record = Entrez.read(search_handle)
        search_handle.close()
        time.sleep(0.35)

        if search_record['IdList']:
            assembly_uid = search_record['IdList'][0]
            link_handle = Entrez.elink(dbfrom='assembly', id=assembly_uid, db='biosample')
            link_records = Entrez.read(link_handle)
            link_handle.close()
            time.sleep(0.35)

            if link_records[0]['LinkSetDb']:
                biosample_uid = link_records[0]['LinkSetDb'][0]['Link'][0]['Id']
                fetch_handle = Entrez.efetch(db='biosample', id=biosample_uid, rettype='xml')
                xml_data = fetch_handle.read()
                fetch_handle.close()
                time.sleep(0.35)

                root = ElementTree.fromstring(xml_data)
                latitude_longitude = None
                environmental_sample_flag = None
                isolation_source = None
                geographic_location_name = None

                attributes = root.findall('.//Attribute')
                for attribute in attributes:
                    attribute_name = attribute.get('attribute_name', '').lower()
                    if attribute_name == 'lat_lon':
                        latitude_longitude = attribute.text
                    elif attribute_name == 'environmental_sample':
                        environmental_sample_flag = attribute.text
                    elif attribute_name == 'isolation_source':
                        isolation_source = attribute.text
                    elif attribute_name == 'geo_loc_name':
                        geographic_location_name = attribute.text

                if environmental_sample_flag and environmental_sample_flag.lower() == 'true':
                    environmental_sample_flag = 'Environmental sample'

                biome_description = isolation_source or environmental_sample_flag
                latitude, longitude = parse_latitude_longitude(latitude_longitude)

                return {
                    'GenomeID': assembly_accession,
                    'Location': geographic_location_name,
                    'BiomeDistribution': categorize_biome(biome_description),
                    'Latitude': latitude,
                    'Longitude': longitude
                }
    except Exception as e:
        print(f"Erro ao buscar dados para o assembly {assembly_accession}: {e}")
    return None

# Função para analisar coordenadas de latitude e longitude
def parse_latitude_longitude(lat_lon_str):
    if lat_lon_str and isinstance(lat_lon_str, str):
        try:
            parts = lat_lon_str.strip().replace(',', '.').split()
            if len(parts) >= 4:
                latitude_value = float(parts[0])
                latitude_direction = parts[1].upper()
                longitude_value = float(parts[2])
                longitude_direction = parts[3].upper()

                if latitude_direction == 'S':
                    latitude_value = -latitude_value
                if longitude_direction == 'W':
                    longitude_value = -longitude_value

                return latitude_value, longitude_value
        except Exception as e:
            print(f"Não foi possível analisar LatitudeLongitude '{lat_lon_str}': {e}")
    return None, None

# Função para categorizar o tipo de bioma
def categorize_biome(biome_description):
    if biome_description is None:
        return 'Unknown'
    biome_description = biome_description.lower()
    if 'soil' in biome_description:
        return 'Soil'
    elif 'marine' in biome_description or 'sea' in biome_description or 'ocean' in biome_description:
        return 'Marine'
    elif 'lake' in biome_description or 'freshwater' in biome_description:
        return 'Lake'
    elif 'waste' in biome_description or 'wastewater' in biome_description or 'sewage' in biome_description:
        return 'Wastewater'
    elif 'host' in biome_description or 'root' in biome_description or 'nodule' in biome_description:
        return 'Host-Associated'
    elif 'hypersaline' in biome_description:
        return 'Hypersaline'
    elif 'reef' in biome_description:
        return 'Reef'
    elif 'environmental sample' in biome_description:
        return 'Environmental sample'
    else:
        return 'Other'

# Carregar o arquivo de dados diretamente
data_file_path = 'table_example.tsv'  # Substitua pelo seu arquivo TSV
metadata_df = pd.read_csv(data_file_path, sep='\t')

# Verificar se as colunas 'Latitude' e 'Longitude' estão no DataFrame, caso contrário, adicioná-las
if 'Latitude' not in metadata_df.columns or 'Longitude' not in metadata_df.columns:
    metadata_df['Latitude'] = None
    metadata_df['Longitude'] = None

# Preencher os metadados para cada GenomeID
for idx, row in metadata_df.iterrows():
    metadata = fetch_biosample_metadata(row['Assembly'])
    if metadata:
        metadata_df.at[idx, 'Latitude'] = metadata['Latitude']
        metadata_df.at[idx, 'Longitude'] = metadata['Longitude']
        metadata_df.at[idx, 'Location'] = metadata['Location']
        metadata_df.at[idx, 'BiomeDistribution'] = metadata['BiomeDistribution']

# Função para criar a tabela interativa com links para GenomeID
def create_interactive_table():
    html = "<table><tr><th>GenomeID</th><th>Location</th><th>Biome Distribution</th></tr>"
    for idx, row in metadata_df.iterrows():
        genome_id = row['Assembly']
        html += f"<tr><td><a href='#' onclick='highlight_marker(\"{genome_id}\")'>{genome_id}</a></td>"
        html += f"<td>{row['Location']}</td><td>{row['BiomeDistribution']}</td></tr>"
    html += "</table>"
    return HTML(html)

# Criar o mapa com uma visualização global inicial
m = Map(center=(0, 0), zoom=2, basemap=basemap_to_tiles(basemaps.OpenStreetMap.Mapnik))

# Função para exibir o popup no mapa para o GenomeID selecionado
def display_metadata_popup(marker):
    with output:
        output.clear_output()
        display(HTML(marker.popup_content))

# Função auxiliar para configurar o clique do marcador
def add_marker_click(marker):
    def marker_click_handler(*args, **kwargs):  # Modificado para aceitar argumentos extras
        display_metadata_popup(marker)
    marker.on_click(marker_click_handler)

# Criar e adicionar marcadores no mapa com base nos dados da tabela
output = Output()
genome_markers = {}

for idx, row in metadata_df.iterrows():
    latitude = row['Latitude']
    longitude = row['Longitude']
    genome_id = row['Assembly']
    
    if pd.notnull(latitude) and pd.notnull(longitude):
        popup_content = f"""
        <b>GenomeID:</b> {genome_id}<br>
        <b>Location:</b> {row['Location']}<br>
        <b>Biome Distribution:</b> {row['BiomeDistribution']}<br>
        <b>Latitude:</b> {latitude}<br>
        <b>Longitude:</b> {longitude}
        """
        
        # Criar o marcador
        marker = CircleMarker(location=(latitude, longitude), radius=7, color="blue", fill_opacity=0.7)
        marker.popup_content = popup_content
        marker.genome_id = genome_id
        
        # Configurar o evento de clique no marcador
        add_marker_click(marker)
        
        m.add_layer(marker)
        genome_markers[genome_id] = marker

# Função para destacar a linha na tabela e exibir o popup ao clicar
def highlight_marker(genome_id):
    if genome_id in genome_markers:
        marker = genome_markers[genome_id]
        m.center = marker.location
        display_metadata_popup(marker)

# Criar o gráfico de distribuição de biomas
biome_counts = metadata_df['BiomeDistribution'].value_counts()
fig = px.pie(biome_counts, names=biome_counts.index, values=biome_counts.values, title="Biome Distribution")
fig_html = to_html(fig, full_html=False)
graph = HTML(fig_html)

# Criar a tabela e organizar o layout
interactive_table = create_interactive_table()
display(VBox([HBox([interactive_table, graph]), m, output]))

