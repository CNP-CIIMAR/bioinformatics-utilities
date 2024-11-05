# Autor: Leandro de Mattos Pereira
# 06 novembro 2024
# Actual: CNP Team - Leam Lab 

import subprocess
import sys
import time
from ete3 import NCBITaxa
from Bio import Entrez
import xml.etree.ElementTree as ElementTree

# Inicializa o NCBITaxa e define o email para o Entrez
ncbi = NCBITaxa()
Entrez.email = 'seu.email@dominio.com'  # Substitua pelo seu email real

if len(sys.argv) < 3:
    print("Uso: script.py <input_file_with_assembly_ids> <output_file>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]
filtered_output_file = "filtered_" + output_file  # Arquivo filtrado adicional

# Contadores para rastrear informações recuperadas
biome_count = 0
lat_lon_count = 0

# Função para obter a linhagem taxonômica com base no Organism Tax ID
def get_lineage(tax_id):
    try:
        lineage = ncbi.get_lineage(tax_id)
        names = ncbi.get_taxid_translator(lineage)
        lineage_str = "; ".join([names[t] for t in lineage])
        return lineage_str
    except Exception as e:
        print(f"Erro ao obter linhagem para Tax ID {tax_id}: {e}")
        return "Desconhecido"

# Função para buscar metadados adicionais do NCBI BioSample
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

                # Extrai atributos específicos do XML
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

                # Atualiza contadores para as informações recuperadas
                global biome_count, lat_lon_count
                if biome_description:
                    biome_count += 1
                if latitude is not None and longitude is not None:
                    lat_lon_count += 1

                return {
                    'Location': geographic_location_name or "Unknown",
                    'BiomeDistribution': categorize_biome(biome_description),
                    'Latitude': latitude if latitude is not None else "",
                    'Longitude': longitude if longitude is not None else ""
                }
    except Exception as e:
        print(f"Erro ao buscar dados para o assembly {assembly_accession}: {e}")
    return {
        'Location': "Unknown",
        'BiomeDistribution': "Unknown",
        'Latitude': "",
        'Longitude': ""
    }

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

# Função para categorizar o tipo de bioma de acordo com as descrições GOLD
def categorize_biome(biome_description):
    if biome_description is None:
        return 'Unknown'
    biome_description = biome_description.lower()
    # Categorias de bioma com base em descrições GOLD
    if any(term in biome_description for term in ['soil', 'forest', 'desert', 'savanna']):
        return 'Terrestrial'
    elif any(term in biome_description for term in ['marine', 'sea', 'ocean', 'coastal']):
        return 'Marine'
    elif any(term in biome_description for term in ['lake', 'freshwater', 'river', 'pond']):
        return 'Freshwater'
    elif any(term in biome_description for term in ['waste', 'wastewater', 'sewage']):
        return 'Wastewater'
    elif any(term in biome_description for term in ['host', 'symbiont', 'root', 'nodule']):
        return 'Host-Associated'
    elif 'hypersaline' in biome_description:
        return 'Extreme - Hypersaline'
    elif 'hot spring' in biome_description or 'thermal' in biome_description:
        return 'Extreme - Thermal'
    elif 'acidic' in biome_description or 'alkaline' in biome_description:
        return 'Extreme - Acidic/Alkaline'
    elif any(term in biome_description for term in ['reef', 'coral']):
        return 'Reef'
    elif 'environmental sample' in biome_description:
        return 'Environmental Sample'
    else:
        return 'Other'

# Escreve o cabeçalho nos arquivos de saída
with open(output_file, 'w') as out_file, open(filtered_output_file, 'w') as filtered_out_file:
    header = (
        "Assembly Accession\tOrganism Name\tOrganism Common Name\tOrganism Tax ID\tLineage\t"
        "Assembly Level\tBioProject Accession\tBioSample Accession\tGC Percent\t"
        "Total Sequence Length\tSequencing Technology\tRelease Date\tCollection Date\t"
        "BioSample Description\tLocation\tBiomeDistribution\tLatitude\tLongitude\n"
    )
    out_file.write(header)
    filtered_out_file.write(header)

# Lê o arquivo de entrada e processa cada Assembly ID
with open(input_file, 'r') as in_file:
    next(in_file)  # Pula o cabeçalho do arquivo de entrada

    for line in in_file:
        assembly_id = line.strip().split('\t')[0]
        command = (
            f"./datasets summary genome accession {assembly_id} --as-json-lines | "
            f"./dataformat tsv genome --fields accession,organism-name,organism-common-name,"
            f"organism-tax-id,assminfo-level,assminfo-bioproject,assminfo-biosample-accession,"
            f"assmstats-gc-percent,assmstats-total-sequence-len,assminfo-sequencing-tech,"
            f"assminfo-release-date,assminfo-biosample-collection-date,"
            f"assminfo-biosample-description-title"
        )

        # Executa o comando e captura a saída
        result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        if result.returncode == 0:
            output_lines = result.stdout.decode('utf-8').splitlines()

            if output_lines and output_lines[0].startswith("Assembly Accession"):
                output_lines = output_lines[1:]
            
            for output_line in output_lines:
                fields = output_line.split('\t')
                tax_id = fields[3]
                lineage = get_lineage(int(tax_id)) if tax_id.isdigit() else "Desconhecido"
                
                # Busca metadados adicionais do BioSample
                metadata = fetch_biosample_metadata(assembly_id)
                
                # Constrói a linha enriquecida com metadados adicionais ao final
                enriched_line = (
                    "\t".join(fields[:4]) + f"\t{lineage}\t" +
                    "\t".join(fields[4:]) + f"\t{metadata['Location']}\t{metadata['BiomeDistribution']}\t"
                    f"{metadata['Latitude']}\t{metadata['Longitude']}"
                )
                
                # Escreve no arquivo principal
                with open(output_file, 'a') as out_file:
                    out_file.write(enriched_line + "\n")
                
                # Se ambos BiomeDistribution e Latitude/Longitude estiverem preenchidos, escreve no arquivo filtrado
                if metadata['BiomeDistribution'] != "Unknown" and metadata['Latitude'] and metadata['Longitude']:
                    with open(filtered_output_file, 'a') as filtered_out_file:
                        filtered_out_file.write(enriched_line + "\n")
        else:
            print(f"Erro ao processar {assembly_id}: {result.stderr.decode('utf-8')}")

# Imprime o resumo de informações recuperadas
print(f"Total de Genome IDs com BiomeDistribution preenchido: {biome_count}")
print(f"Total de Genome IDs com Latitude e Longitude preenchidos: {lat_lon_count}")
print(f"Total de Genome IDs salvos na tabela filtrada: {biome_count if biome_count == lat_lon_count else min(biome_count, lat_lon_count)}")
