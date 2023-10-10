import requests
import time
import argparse
import pandas as pd
#Author: Leandro de Mattos Pereira
# CNP -Team
#Team Leader. Pedro Leão
def submit_mapping_request(ids_list):
    """Submit the mapping request to UniProt."""
    url = "https://rest.uniprot.org/idmapping/run"
    data = {
        'from': "RefSeq_Protein",
        'to': "UniProtKB",
        'ids': ",".join(ids_list)
    }
    headers = {
        'Content-Type': 'application/x-www-form-urlencoded'
    }
    response = requests.post(url, data=data, headers=headers)
    if response.status_code == 200:
        return response.json().get('jobId', None)
    else:
        response.raise_for_status()

def get_mapping_status(job_id):
    """Check the status of the mapping job."""
    url = f"https://rest.uniprot.org/idmapping/status/{job_id}"
    response = requests.get(url)
    if response.status_code == 200:
        status = response.json().get('jobStatus', None)
        if not status:
            print(f"Unexpected response: {response.json()}")
            return None
        return status
    else:
        response.raise_for_status()

def fetch_mapping_results(job_id):
    """Fetch the mapping results including additional values with handling of missing data."""
    url = f"https://rest.uniprot.org/idmapping/uniprotkb/results/{job_id}"
    response = requests.get(url)
    if response.status_code == 200:
        json_data = response.json()
        from_ids = []
        to_ids = []
        entry_types = []
        descriptions = []
        primary_accessions = []
        proteomes = []
        id_values = []
        lineage = []
        cofactor = []
        interpro = []
        supfam = []
        sequence = []
        go_term = []
        rhea_reaction_id = []
        go_term_rhea = []
        ec_classes = []
        kegg = []

        for result in json_data['results']:
            from_ids.append(result['from'])
            to_info = result['to']
            to_ids.append(to_info['uniProtkbId'])
            entry_types.append(to_info['entryType'])
            descriptions.append(to_info['proteinDescription']['submissionNames'][0]['fullName']['value'])

            # Extract additional values with handling of missing data
            primary_accessions.append(to_info.get('primaryAccession', ''))
            proteomes.append(to_info.get('proteomes', ''))
            id_values.append(to_info.get('id', ''))
            lineage.append(to_info.get('lineage', ''))
            cofactor.append(to_info.get('cofactor', ''))
            interpro.append(to_info.get('interpro', ''))
            supfam.append(to_info.get('supfam', ''))
            sequence.append(to_info.get('sequence', ''))
            go_term.append(to_info.get('goTerm', ''))
            rhea_reaction_id.append(to_info.get('rheaReactionId', ''))
            go_term_rhea.append(to_info.get('goTermRhea', ''))
            ec_classes.append(to_info.get('ecClasses', ''))
            kegg.append(to_info.get('kegg', ''))

        return {
            'From_ID': from_ids,
            'To_ID': to_ids,
            'Entry_Type': entry_types,
            'Description': descriptions,
            'Primary_Accession': primary_accessions,
            'Proteomes': proteomes,
            'ID': id_values,
            'Lineage': lineage,
            'COFACTOR': cofactor,
            'InterPro': interpro,
            'SUPFAM': supfam,
            'Sequence': sequence,
            'GoTerm': go_term,
            'RHEA_Reaction_ID': rhea_reaction_id,
            'GO_Term_RHEA': go_term_rhea,
            'EC_Classes': ec_classes,
            'Kegg': kegg,
        }
    else:
        response.raise_for_status()

def main(input_file, output_file):
    # Read RefSeq IDs from file
    with open(input_file, 'r') as file:
        ids_list = [line.strip() for line in file.readlines()]

    # Submit mapping request
    job_id = submit_mapping_request(ids_list)
    if not job_id:
        print("Failed to retrieve job ID from UniProt.")
        return

    # Poll for job status
    status = get_mapping_status(job_id)
    while status and status != "FINISHED":
        time.sleep(10)  # Wait for 10 seconds before polling again
        status = get_mapping_status(job_id)

    # Fetch results (moved this part inside main)
    results_dict = fetch_mapping_results(job_id)

    # Convert to DataFrame
    df = pd.DataFrame(results_dict)

    # Save to TSV (tab-separated values)
    df.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Map RefSeq Protein IDs to UniProtKB IDs")
    parser.add_argument("input_filepath", help="Path to the input file containing RefSeq Protein IDs")
    parser.add_argument("output_filepath", help="Path to the output file to save the mappings")
    args = parser.parse_args()

    main(args.input_filepath, args.output_filepath)
