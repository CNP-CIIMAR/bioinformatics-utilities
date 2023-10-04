#Author: Leandro de Mattos Pereira
#CNP- team - Pedro Leão - Team Leader
import subprocess
import argparse
import os

def run_command(command):
    process = subprocess.Popen(command, shell=True)
    process.communicate()

def main(args):
    # 1 & 2) Build the alignment
    mafft_command = f"mafft --thread 10 --maxiterate 1000 --localpair {args.input_fasta} > {args.output_mafft}"
    run_command(mafft_command)

    # 3) Build hmm database
    hmmbuild_command = f"hmmbuild {args.output_hmm} {args.output_mafft}"
    run_command(hmmbuild_command)

    # Display statistics of the hmm database
    hmmstat_command = f"hmmstat {args.output_hmm}"
    run_command(hmmstat_command)

    # 5) hmmpress: binary compression of files and indexing
    hmmpress_command = f"hmmpress {args.output_hmm}"
    run_command(hmmpress_command)

    # 6) Perform hmmsearch on each FASTA file in the provided directory
    fasta_files = [f for f in os.listdir(args.fasta_directory) if f.endswith('.fasta')]
    for fasta in fasta_files:
        fasta_path = os.path.join(args.fasta_directory, fasta)
        out_file = os.path.join(args.fasta_directory, fasta.replace(".fasta", ".out"))
        hmmsearch_command = f"hmmsearch {args.output_hmm} {fasta_path} > {out_file}"
        run_command(hmmsearch_command)

    # 7) Generate table of hits
    domtblout_command = f"hmmsearch --domtblout {args.domtblout} {args.output_hmm} {args.ref_db}"
    run_command(domtblout_command)

    # 8 & 8.1) Index reference database
    index_db_command = f"esl-sfetch --index {args.ref_db}"
    run_command(index_db_command)

    # 9) Extract homologous sequences
    grep_awk_command = f"grep -v '^#' {args.domtblout} | awk '{{print $1}}' | esl-sfetch -f {args.ref_db} - > {args.homologous_output}"
    run_command(grep_awk_command)

    # 10 & 11) Retrieve subsequences (domains)
    retrieve_subseq_command = f"grep -v '^#' {args.domtblout} | awk '{{print $1\"/\"$20\"-\"$21, $20, $21, $1}}' | esl-sfetch -Cf {args.ref_db} - > {args.domain_output}"
    run_command(retrieve_subseq_command)

    # 12) Run interproscan
    interpro_command = f"./interproscan.sh -i {args.domain_output} -f tsv"
    run_command(interpro_command)

    # 13) Get domain for Phylogenetic tree
    grep_command = f"grep 'FAAL' *.tsv > {args.interpro_kegg}"
    run_command(grep_command)

    cut_bedtools_command = f"cut -f1,7,8 {args.interpro_kegg} > {args.bed_output} && bedtools getfasta -fi {args.ref_db} -bed {args.bed_output} -fo {args.final_output}"
    run_command(cut_bedtools_command)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process fasta files and produce outputs")
    parser.add_argument("--input_fasta", required=True, help="Input fasta file for mafft")
    parser.add_argument("--output_mafft", required=True, help="Output for mafft")
    parser.add_argument("--output_hmm", required=True, help="Output for hmmbuild")
    parser.add_argument("--fasta_directory", required=True, help="Directory containing multiple fasta files for hmmsearch")
    parser.add_argument("--domtblout", required=True, help="Output for domtblout")
    parser.add_argument("--ref_db", required=True, help="Reference database for various commands")
    parser.add_argument("--homologous_output", required=True, help="Output file for extracted homologous sequences")
    parser.add_argument("--domain_output", required=True, help="Output file for domain sequences")
    parser.add_argument("--interpro_kegg", required=True, help="Output for InterProScan's KEGG results")
    parser.add_argument("--bed_output", required=True, help="Output BED file")
    parser.add_argument("--final_output", required=True, help="Final output for extracted sequences based on BED file")

    args = parser.parse_args()
    main(args)
