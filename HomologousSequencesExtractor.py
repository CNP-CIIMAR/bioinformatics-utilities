import os
import argparse
import requests
import time
##Autor: Leandro de Mattos Pereira, Computational Biologist.
## CNP team 04/10/2023
## Pedro Leão team  Leader
def run_interproscan_api(input_fasta_path, output_result_path):
    """
    Execute InterProScan via API for a given fasta file.

    Args:
    - input_fasta_path: path to the input fasta file.
    - output_result_path: path to save the InterProScan result.
    """
    INTERPROSCAN_URL = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/run"
    headers = {"Content-Type": "application/x-www-form-urlencoded"}
    with open(input_fasta_path, 'r') as fasta_file:
        fasta_content = fasta_file.read()
    
    # Submit the job to InterProScan
    data = {
        "sequence": fasta_content,
        "format": "xml",  # You can change this to other supported formats if needed
    }
    response = requests.post(INTERPROSCAN_URL, headers=headers, data=data)
    job_id = response.text

    # Check job status
    status_url = f"https://www.ebi.ac.uk/Tools/services/rest/iprscan5/status/{job_id}"
    while True:
        response = requests.get(status_url)
        status = response.text
        if status == "FINISHED":
            break
        elif status == "ERROR" or status == "FAILED":
            raise Exception("InterProScan job failed.")
        else:
            time.sleep(30)  # Wait for 30 seconds before checking again

    # Fetch the result
    result_url = f"https://www.ebi.ac.uk/Tools/services/rest/iprscan5/result/{job_id}/xml"
    response = requests.get(result_url)
    with open(output_result_path, 'w') as result_file:
        result_file.write(response.text)


def run_mafft(input_file, output_file):
    os.system(f"mafft {input_file} > {output_file}")

def run_hmmbuild(input_file, output_file):
    os.system(f"hmmbuild {output_file} {input_file}")

def run_hmmpress(hmm_file):
    os.system(f"hmmpress {hmm_file}")

def run_hmmsearch(hmm_file, fasta_directory, ref_db, output_directory):
    for root, dirs, files in os.walk(fasta_directory):
        for file in files:
            if file.endswith(".fasta"):
                input_path = os.path.join(root, file)
                domtblout_path = os.path.join(output_directory, f"{file.split('.')[0]}_domtblout.tbl")
                cmd = f"hmmsearch --domtblout {domtblout_path} {hmm_file} {input_path}"
                os.system(cmd)

                # Check if the domtblout file exists before trying to open it
                if os.path.isfile(domtblout_path):
                    with open(domtblout_path, "r") as tbl_file:
                        lines = tbl_file.readlines()
                        seqs = [line.split()[0] for line in lines if not line.startswith("#")]

                    seqs_file_path = os.path.join(output_directory, "seqs_to_extract.txt")
                    with open(seqs_file_path, "w") as seqs_file:
                        seqs_file.write("\n".join(seqs))

                    # Extract sequences using esl-sfetch
                    extracted_seqs_path = os.path.join(output_directory, f"{file.split('.')[0]}_extracted.fasta")
                    os.system(f"esl-sfetch -f {ref_db} {seqs_file_path} > {extracted_seqs_path}")
                    interproscan_result_path = os.path.join(output_directory, f"{file.split('.')[0]}_interproscan.xml")
                    run_interproscan_api(extracted_seqs_path, interproscan_result_path)
                else:
                    print(f"Domtblout file not found for {input_path}.")

def main(args):
    # Verificar se o diretório de saída existe. Se não, crie-o.
    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)

    run_mafft(args.input_fasta, args.output_mafft)
    run_hmmbuild(args.output_mafft, args.output_hmm)
    run_hmmpress(args.output_hmm)
    run_hmmsearch(args.output_hmm, args.fasta_directory, args.ref_db, args.output_directory)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Script to run a series of bioinformatics commands.')
    parser.add_argument('--input_fasta', required=True, help='Path to the input fasta file.')
    parser.add_argument('--output_mafft', required=True, help='Path to the output mafft file.')
    parser.add_argument('--output_hmm', required=True, help='Path to the output hmm file.')
    parser.add_argument('--fasta_directory', required=True, help='Directory containing fasta files.')
    parser.add_argument('--ref_db', required=True, help='Path to the reference database fasta file.')
    parser.add_argument('--output_directory', required=True, help='Directory to save output files.')

    args = parser.parse_args()
    main(args)
