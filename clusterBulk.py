# File: clusterBulk
# Authors: 
# Leandro de Mattos Pereira
# Pedro Leao - Team Leader Researchers - CNP Team
# Description: use ClusterPackage to Cluster all Sequences in dir Sequences
# ATTENTION Files must end with fasta, or modify the code!!!
import sys
import os


from clusterpackage import cluster, balance_hash, parse_blast, interpretate_hash, write_result, get_list, occur_error, get_sequence_of_fasta_file, get_amount_seq

Sequences = "/investgenomicaarea/Leandro/Metagenomica/KEGG_Leandro/KEGG_2020/KEGG_2020cp/Result"

if len(sys.argv) < 2:
    sys.exit("Usage: python clusterBulk.py <bits_threshold> [statisticfile]")

Threshold = sys.argv[1]
Statistics = sys.argv[2] if len(sys.argv) > 2 else None

try:
    files = os.listdir(Sequences)
except OSError:
    sys.exit("Error opening Sequences directory")

print("ATTENTION: Only files with the .fasta extension will be used")

for file in files:
    if file.endswith(".fasta"):
        cluster(os.path.join("Sequences", file), Threshold, Statistics)
