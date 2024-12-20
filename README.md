

# Python Bioinformatics Utilities:

_____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________


# Script 1: Subdirectory File Copy Script
This script copy_files_with_folder_name_gz.py is a Python utility for copying files from subdirectories to a new directory while adding the name of the subdirectory to the beginning of the file name.
It also supports decompressing and copying gzip-compressed files.**

## Usage:
To use this script, follow these steps:

1. Clone or download the repository to your local machine.

2. Make sure you have Python 3 installed.

3. Open a terminal or command prompt and navigate to the directory where you downloaded the script.

4. Run the script with the following command:

   ```python copy_files_with_folder_name_gz.py <source_directory> <destination_directory>```

If you need any help or have questions, feel free to open an issue in the GitHub repository.

# Script 2: Get genome ID and Taxonomic rank from Protein ID

get_genome_id_and_taxonomic_rank_from_protein.id.py


#  Script 3: File Comparison

This Python script compares the contents of two text files and identifies the different lines between them.

## Usage:
Make sure you have Python installed on your system. The script was developed using Python 3.

1. Open a terminal or command prompt.

2. Run the script supplying the names of the files you want to compare as command line arguments:

    ```
    python comparar_arquivos.py arquivo1.txt arquivo2.txt
    ```

    Be sure to replace `file1.txt` and `file2.txt` with the actual names of the files you want to compare. The files must be in the same directory as the script or provide the full path to them.

3. The script will compare the files and display the different lines found. If there are no differences, the message "The files have the same content."

## Example:
Suppose you have two text files: `file1.txt` and `file2.txt`. Their content is as follows:

**file1.txt**

PYQ20818.1 

**file2.txt**

PYQ20818.1

# Script 4: InterProScan to iTOL Converter

This script converts InterProScan results into an iTOL (Interactive Tree Of Life) compatible format.

## Requirements:
- Python 3.x
- Pandas (instalável via `pip install pandas`)

## Usage:
Run the script by providing the path to the input file containing the InterProScan results and the path to the output file that will be generated in iTOL compatible format.

Example of use:
```python interpro_to_itol.py arquivo_entrada.txt arquivo_saida.txt```

- Input file format
The input file must be a text file containing InterProScan results, where each line represents a domain annotation on a protein. Columns must be separated by tabs and in the following order:

Protein.accession\tSequence.length\tshape\tStart.location\tStop.location

*Protein.accession: Identificador da proteína.*
*Sequence.length: Comprimento da sequência da proteína.*
*shape: Forma do domínio.*
*Start.location: Posição de início do domínio na sequência.*
*Stop.location: Posição de término do domínio na sequência.*

The script will generate an output file in iTOL compatible format. The output file will have the following columns:

Protein.accession\tSequence.length\tshape\tStart.location\tStop.location\tcolor

Protein.accession: Identificador da proteína.
Sequence.length: Comprimento da sequência da proteína.


If you need any help or have questions, feel free to open an issue in the GitHub repository.

# Script 5: FASTA Sequence Filter
This is a Python script that filters a FASTA file containing DNA or protein sequences based on a list of sequence IDs.

## Requirements:
- Python 3.x
- Biopython library

## Installation:
1. Clone the repository:

 ```git clone https://github.com/mattoslmp/CNP-Ciimar.git ```
   
2.  Install the Biopython library:

  ```pip install biopython  ```

  ```python3 fasta_sequence_filter.py input.fasta output.fasta ids.txt --exclude ```

input.fasta: The path to the input FASTA file containing the sequences to be filtered.
output.fasta: The path to the output FASTA file where the filtered sequences will be saved.
ids.txt: The path to the file containing the sequence IDs to filter.
--exclude: Optional flag to indicate that the sequences with matching IDs should be excluded. Alternatively, you can use --keep to keep only the sequences with matching IDs.
Make sure you have the Biopython library installed (pip install biopython) before running the script.


**Example**:
To filter the sequences in input.fasta based on the IDs in ids.txt and generate the filtered sequences in output.fasta, excluding the matching IDs, run the following command:
  ```python script.py input.fasta output.fasta ids.txt --exclude``` 

**The script will provide information about the number of sequences in the input file and the number of sequences remaining after filtration.**

If you need any help or have questions, feel free to open an issue in the GitHub repository.

## Script 6: Taxonomic Rank Retrieval

**This is a Python script that utilizes the Biopython library to retrieve taxonomic classification information for protein sequences. It takes a list of protein accessions as input and returns the species name and corresponding taxonomic lineage.**

## Prerequisites
- Python 3.x
- Biopython library

## Installation
Make sure you have Python 3.x installed on your system.
Install the Biopython library by executing the following command:
pip install biopython

## Usage
Create a text file containing the desired protein accessions, with each accession on a separate line. For example, you can name the file protein_accessions.txt and include the accessions as follows:

NP_000001.1
NP_000002.1
NP_000003.1
Run the Python script script.py by providing the input file containing the protein accessions and the output file to save the results. For example:

python script.py protein_accessions.txt output.csv
Make sure to replace protein_accessions.txt with the path to your input file and output.csv with the desired path for the output file.

Wait for the script to execute. The results will be saved in the specified output file.

## Notes
If there are errors while retrieving the data for a specific protein accession, a message will be displayed indicating the error.
The script will save the results in a tab-separated values (.tsv) file for easy further analysis. You can adjust this setting by editing the line df.to_csv(output_file_path, sep='\t', index=False).


## Script 7:
- Run Barnap
This script is designed to automate the usage of the Barnap program, which predicts the location of ribosomal RNA (rRNA) genes in genomes. The script takes a directory containing multiple DNA sequence files in FASTA format and executes Barnap on each file, saving the results in a designated output directory.

## Installation:
Before running the script, make sure to install the Barnap program. Here are the installation instructions for different platforms:

## Conda
If you have Conda or Miniconda installed, you can use the following command to install Barnap:

 ```shell
conda install -c bioconda -c conda-forge barrnap
 ```
## Homebrew
For macOS users, install Homebrew, and then run the following command to install Barnap:

 ```shell
brew install brewsci/bio/barrnap
 ```
Source
To install the latest version directly from GitHub, follow these steps:

Change to your home directory:
 ```shell
cd $HOME
 ```
Clone the Barnap repository:
 
 ```shell
git clone https://github.com/tseemann/barrnap.git
 ```
- Change to the barrnap/bin directory:
 ```shell
cd barrnap/bin
 ```
- Add the bin directory to your PATH environment variable. You can add the following line to your shell configuration file (e.g., .bashrc, .bash_profile, or .zshrc):

 ```shell
export PATH="$HOME/barrnap/bin:$PATH"
 ```
## Usage:
**To use the script, follow these steps:**

- Ensure that the Barnap program is properly installed on your system (as described above).

- Save the provided Python script to a file, e.g., barnap_script.py.

- Open a terminal and navigate to the directory containing the script.

- Run the script with the following command, providing the directory path as an argument:

 ```shell
python run_barnap.py <directory>
 ```
- Replace <directory> with the path to the directory containing the DNA sequence files in FASTA format that you want to analyze.

**Output**
The script will create a directory named "16SrRNAmultigenomes" (if it doesn't already exist) in the same directory as the input files. The results of the Barnap analysis for each input file will be saved in this directory with the filename formatted as <input_filename>_processado.rrna.fa.

## Credits:
The Barnap program was developed by Torsten Seemann. For more information about Barnap, please refer to the official repository.
If you have any questions or need further assistance, please don't hesitate to reach out.


## Script 8: Get Genomes data from BV-BRC database

## Pré-requisitos
- Python 3.10

## Instalação

1. Clone o repositório:

2.  ```git clone https://github.com/mattoslmp/CNP-Ciimar.git ```
3.  Create one list of Genome ID of BV-BRC database
4.  Save the ID into the file named Genomes_id

Example of content of the file Genome_id:  
 
1033813.3  
1038927.31  
1038927.40  
1038927.41  
1038927.45  

6.  ```python3 get_genomes_by_id_bv-brcdb.py ```

## Script 9: ProteinHMM

A tool to perform hmmsearch searches on .faa files using HMM models from one provided directory and then convert the results to .tsv format.

proteinsearchhmm.py MMER Table Converter:

** This Python script performs hmmsearch searches on .faa files using provided HMM models and then converting HMMER output tables into a more user-friendly format .tsv for easy viewing and analysis.

# Requirements
- Python 3.6 or higher
- pandas library
- Installation
- Clone or download this repository to your local machine.
**Install the required dependencies using pip:
- pip install pandas

## Usage
Run the script using the following command:
 ```shell
python proteinhmm.py <models_dir> <fastas_dir> <output_dir>
 ```
<models_dir> is the directory containing the HMM models (.hmm).
<fastas_dir> is the directory containing the FASTA files (.faa).
<output_dir> is the directory where results will be saved.

# Example
Suppose you have your models in /home/user/models/, your .faa files in /home/user/fastas/, and you wish to save the results in /home/user/results/. You would run:
 ```shell
python proteinsearchhmm.py /home/user/models/ /home/user/fastas/ /home/user/results/
 ```
## Features
Searches FASTA sequences using HMM models with the hmmsearch command.
Converts HMMER's tblout output format into a .tsv format.
Extracts species from each row of the "query_description" column using regular expressions.

## Supported Formats
The script supports two HMMER output formats:

tblout: The output format generated by hmmsearch or hmmscan. The converted TSV file will contain the following columns:

- target_name: The name of the target sequence or profile.
- target_accession: The accession of the target sequence or profile.
- query_name: The name of the query sequence or profile.
- query_accession: The accession of the query sequence or profile.
- full_sequence_e-value: E-value of the full sequence.
- full_sequence_score: Score of the full sequence.
- full_sequence_bias: Bias of the full sequence.
- best_domain_e-value: E-value of the best domain.
- best_domain_score: Score of the best domain.
- best_domain_bias: Bias of the best domain.
- domain_number_estimation_exp: Expected number of domains.
- domain_number_estimation_reg: Regularized number of domains.
- domain_number_estimation_clu: Clustering number of domains.
- domain_number_estimation_ov: Overlapping number of domains.
- domain_number_estimation_env: Enveloping number of domains.
- domain_number_estimation_dom: Domain number.
- domain_number_estimation_rep: Repeated number of domains.
- domain_number_estimation_inc: Inclusion number of domains.
- query_description: Description of the query.
- domtblout: The output format generated by hmmscan with --domtblout option. The converted TSV file will contain the following columns:

- target_name: The name of the target sequence or profile.
- target_accession: The accession of the target sequence or profile.
- target_length: Length of the target sequence.
- query_name: The name of the query sequence or profile.
- query_accession: The accession of the query sequence or profile.
- query_length: Length of the query sequence.
- full_sequence_e-value: E-value of the full sequence.
- full_sequence_score: Score of the full sequence.
- full_sequence_bias: Bias of the full sequence.
- this_domain_number: Number of this domain.
- total_domains: Total number of domains.
- this_domain_e-value: E-value of this domain.
- this_domain_i-value: I-value of this domain.
- this_domain_score: Score of this domain.
- this_domain_bias: Bias of this domain.
- hmm_coord_from: Start coordinate of the HMM alignment.
- hmm_coord_to: End coordinate of the HMM alignment.
- ali_coord_from: Start coordinate of the alignment in the target sequence.
- ali_coord_to: End coordinate of the alignment in the target sequence.
- env_coord_from: Start coordinate of the envelope.
- env_coord_to: End coordinate of the envelope.
- acc: Accession of the query sequence or profile.
- query_description: Description of the query.
- Species

## Notes
The script was developed using HMMER v3.2.1.

The assert_acceptable_arguments function is used to validate the input arguments. If an invalid argument is provided for program or format, a ValueError will be raised.

**The script utilizes the pandas library for data manipulation and provides a convenient way to extract species information from the query_description column using regular expressions.**

## Script 9.1: recov_fasta_protein_hmmer.py Recovery fasta from proteinsearchhmm output FASTA Sequence Extractor

## Description
This Python script is designed to extract specific FASTA sequences from a FASTA file (.fasta or .faa) based on identifiers (IDs) provided in a TSV file. The script reads the IDs from the TSV file, searches for these IDs in the FASTA file, and writes the corresponding sequences to a new FASTA file.

## Requirements
- Python 3.x

## Usage
To use this script, you need a TSV file containing the desired IDs in the first column and a FASTA file from which the sequences will be extracted.

### Command
```bash
python recov_fasta_protein_hmmer.py <TSV_FILE> <FASTA_FILE> <OUTPUT_FILE>
 ```
- `TSV_FILE`: Path to the TSV file containing IDs in the first column.
- `FASTA_FILE`: Path to the FASTA file from which sequences will be extracted.
- `OUTPUT_FILE`: Path to the output FASTA file that will contain the extracted sequences.

### Example

```bash
python extract_fasta.py identifiers.tsv sequences.fasta output.fasta
```

## Key Functions
1. **Reading the TSV File**: Extracts IDs from the first column of the TSV file.
2. **Searching in the FASTA File**: Searches for corresponding sequences in the FASTA file using the extracted IDs.
3. **Writing Sequences**: Writes all found sequences into a single output FASTA file.

## Script 10: get_taxonomic_rank_specie.py:

**Este script em Python é usado para baixar os níveis taxonômicos de uma lista de espécies e salvar os resultados em um arquivo TSV.**

```bash
python get_taxonomic_rank_specie.py lista_species
```
This Python script is used to download the taxonomic ranks of a list of species and save the results in a TSV file.

## Requirements

Make sure you have the following software installed in your environment before running the script:

- Python 3 (version 3.6 or higher)
- `ete3` package (you can install it using `pip install ete3`)

## How to Use
1. Download the `get_taxonomic_rank_specie.py` file to your local directory.
2. Open a terminal or command prompt and navigate to the directory where you downloaded the `get_taxonomic_rank_specie.py` file.
3. Ensure you have a prepared species list file. The file should have one species per line.
4. Run the following command to start the script and download the taxonomic ranks:

```bash
   python get_taxonomic_rank_specie.py species_list_file.txt output_file.tsv
 ```

**
# Script 11: filter_dup_gca_gcf_keep_gcf.py

**"This is a Python script for filtering Tab-separated values (TSV) obtained from NCBI datasets/genomes (https://www.ncbi.nlm.nih.gov/datasets/genome/). The script creates a new table and retains only genomes with GCF (REFSEQ) when this specie genome have both GCF (format for RefSeq, NCBI-derived assembly accessions) and GCA (format for GenBank primary assembly accessions), or keeps only GCA in the output table when there is no GCF available.**

## Requirements

- Python 3.x

## How to Use
1. Make sure you have Python 3.x installed on your system.
2. Download this repository or copy the contents of the `script.py` file.
3. Open a terminal or command prompt.
4. Navigate to the directory where the `script.py` file is located.
5. Run the following command:

```bash
python filter_dup_gca_gcf_keep_gcf.py <table_path> <output_path>
 ```
## Comments
- Make sure the table is in the proper format, with columns acoustic correctly by tab (\t).
- Make sure you have the necessary permissions to write the output file to the specified path.
- Be sure to replace the `<path_to_table>` and `<path_to_output>` sections with the following information

# Script 12: check_genome_folder_subdirectories.py

## Subdirectory Verifier
## Description
This Python script verifies the presence of subdirectories based on a provided list of genomes. It compares the names of subdirectories in a specified directory with the names listed in a file. The script then generates a .txt file containing information on which subdirectories are present and which are not.

## Requirements
Python 3.x
Usage

# Installation
Clone the repository to your local computer:
```bash
git clone https://github.com/CNP-CIIMAR/bioinformatics-utilities
 ```
## Execution
- Navigate to the directory where the check_subdirectories.py script is located and run the following command:
```bash
python check_subdirectories.py /path/to/directory /path/to/coluna1_genomas_download_19_set2023
 ```
# Script 13: Homologous Sequences Extractor: HomologousSequencesExtractor.py

This script serves as a tool to run a series of bioinformatics commands, which include performing alignments with MAFFT, constructing HMM models with HMMER, and searching sequences using HMMER. Additionally, the script also performs analyses on the extracted sequences using InterProScan via its API.

# Prerequisites
Python 3.x
requests library (for the InterProScan API)
MAFFT
HMMER

# InterProScan (optional, if local execution is desired)

## How to Use
To use this script, you need to provide several arguments:
--input_fasta: Path to the input fasta file (Fasta sequences which you want to build a Markov model)
--output_mafft: Path to the MAFFT output file.
--output_hmm: Path to the HMM output file.
--fasta_directory: Directory containing fasta files.
--ref_db: Path to the reference database fasta file.
--output_directory: Directory to save output files.

## Usage example:
```bash
python3 HomologousSequencesExtractor.py --input_fasta your_input_file.fasta --output_mafft output_mafft_synthetic.fasta --output_hmm output_hmm_synthetic.hmm --fasta_directory ./test_fastas/ --ref_db reference_database.fasta --output_directory output_test
 ```
Functions in the code:
run_mafft: Performs sequence alignment using MAFFT.
run_hmmbuild: Constructs an HMM model using HMMER.
run_hmmpress: Presses the HMM file to prepare it for searches.
run_hmmsearch: Searches sequences in a directory using the HMM model.
run_interproscan_api: Executes analyses on extracted sequences using InterProScan via its API.
Notes
When using the run_interproscan_api function, be aware that the InterProScan API might have limitations regarding file size or the number of requests. Additional adaptations might be necessary based on the user's requirements.

# Script 14: UniProtKB ID Mapper: protein_mapping_refseq_uniprotkb.py

This script maps RefSeq Protein IDs to UniProtKB IDs using the UniProt REST API. It fetches various attributes for each protein, such as entry type, description, lineage, and other relevant data.

Dependencies
requests
time
argparse
pandas
# You can install these with:
```bash
pip install requests pandas 
```
## How to Use
**To run the script, you need an input file containing a list of RefSeq Protein IDs and an output file path where the results will be saved in TSV format.**
## Example
# Given an input file named protein_list.txt with the following content:
- Protein.accession
- PYQ20818.1
- PYQ09033.1
- MCP4380582.1
- MBX7220068.1
- MBT5902101.1
- WP_247365613.1
- WP_052809540.1
- KIH99630.1
- HBH72471.1
- MCE9613943.1
- OGL16068.1
- WP_015206076.1
- WP_002624887.1
- BDT34836.1
- NEZ58459.1
- WP_249268147.1
- WP_002799844.1
- WP_264323877.1


# You can map the IDs to UniProtKB IDs and fetch the associated data with:
```bash
python protein_mapping_refseq_uniprotkb.py protein_list.txt output.tsv
 ```
## Data Columns

# The script fetches and outputs the following columns for each RefSeq Protein ID from UniprotkB:

1. **From_ID**: The original RefSeq Protein ID.
2. **To_ID**: The corresponding UniProtKB ID.
3. **Entry_Type**: The type of UniProtKB entry, e.g., Swiss-Prot or TrEMBL.
4. **Description**: Full name of the protein.
5. **Primary_Accession**: The primary accession number of the protein in UniProtKB.
6. **Proteomes**: Information related to the proteomes the protein is a part of.
7. **ID**: Identifier values associated with the protein.
8. **Lineage**: Taxonomic lineage of the protein.
9. **COFACTOR**: Cofactor associated with the protein.
10. **InterPro**: InterPro identifiers associated with the protein.
11. **SUPFAM**: SUPFAM identifiers associated with the protein.
12. **Sequence**: Amino acid sequence of the protein.
13. **GoTerm**: Gene Ontology terms associated with the protein.
14. **RHEA_Reaction_ID**: RHEA reaction identifiers associated with the protein.
15. **GO_Term_RHEA**: GO terms that are linked with RHEA reactions for the protein.
16. **EC_Classes**: Enzyme Commission numbers for the protein.
17. **Kegg**: KEGG database identifiers for the protein.

# Script 15: parser_protein_mapping_output_api_refseq_to_uniprokb.py

**Processing text output  generate from script: protein_mapping_refseq_uniprotkb.py to Excel and TSV friendly format**

This script processes a given text file, structures the data into a pandas DataFrame, and then exports the data to both Excel and TSV formats.

# Features

- Parses text data with tab-separated values.
- Extracts JSON formatted data from the "Sequence" column, and appends it as separate columns.
- Handles errors gracefully, emitting warnings for problematic rows.
- Outputs the structured data to both Excel and TSV formats.

# Requirements:
. pandas
- You can install the required library using:

```bash
pip install pandas
 ```

# Script 16: dbscan_analysis_graph_heatmap.py Overview Results Run_DBCAN4 Heatmap Generator

This is a Python script that generates heatmaps based on data collected from overview files generated by run_dbcan program (https://github.com/linnabrown/run_dbcan) storage in several directory.

## Requirements
- Python 3.x
- Pandas
- Seaborn
- Matplotlib

## How to Use
1. Clone this repository or download the `heatmap_generator.py` script to your system.
2. Execute the script by providing the directory containing the overview files as an argument. For example:
```bash
python dbscan_analysis_graph_heatmapv2.py /path/to/your/directories/withoverviewfiles
 ```
This will generate heatmaps and count tables for different types of data found in the overview files.

- The results will be saved in the same directory where the script was executed.
# Features

- Generates heatmaps for the types of data found in the overview files (DIAMOND, dbCAN_sub, and HMMER).

- Saves count tables in tab-separated values (TSV) format for further analysis with collums: Subfamily       contagem        Specie
- The name of species should be append in the output file named overview of run_dbcan4, after the name **overview** , should be include the name: **_output_Anabaena_aphanizomenioides_LEGE_00250.txt**, then the complete name should be: **overview_output_Anabaena_aphanizomenioides_LEGE_00250.txt** (Example).
- **OBS: It is mandatory that you include the species name in each file overview generated by run_dbscan4**
- Columns in tables and labels on the x-axis of heatmaps are sorted in ascending order.

## License

- This project is licensed under the MIT License. See the LICENSE file for more information.

**If the Bioinformatics scripts utilities were useful, please give proper credit to the authors or include the following citation: "CNP-CIIMAR GitHub [https://github.com/CNP-CIIMAR]."**

# Script 16: dbscan_analysis_graph_barplot.py Overview Results Run_DBCAN4 Barplots Generator

## Overview
This repository contains a Python script that generates horizontal stacked bar plots from data extracted from overview files generated by run_dbcan program (https://github.com/linnabrown/run_dbcan) storage in several directory.

The script is specifically designed to visualize counts of subfamilies across different species in a clear and comprehensible manner using vibrant and distinct colors.

## Prerequisites
Python 3.x
Pandas
Matplotlib
Seaborn
Argparse
You can install the required packages using pip:

```bash
pip install pandas matplotlib seaborn argparse
 ```
## Functions:
The script defines the following primary functions:

- collect_data_from_overviews(directory)
- The name of species must be append in the output file named overview of run_dbcan4, after the name **overview** , must be include the name: **_output_Anabaena_aphanizomenioides_LEGE_00250.txt** (specie name), then the complete name should be: **overview_output_Anabaena_aphanizomenioides_LEGE_00250.txt** (Example).
- **OBS: It is mandatory that you include the species name in each file overview generated by run_dbscan4**
- Columns in tables and labels on the x-axis of heatmaps are sorted in ascending order.

- Takes a directory path as input. This directory have all subdirectories with overview files 
Iterates through all files in the directory and its subdirectories, collecting and organizing data related to 'DIAMOND', 'dbCAN_sub', and 'HMMER'.
Returns three dictionaries, each containing the organized data for 'DIAMOND', 'dbCAN_sub', and 'HMMER', respectively.
plot_stacked_bar(data, title, file_prefix, top_n=40)

- Takes a data dictionary, a title for the plot, a file prefix for the saved plot file, and an optional parameter top_n (default is 40) as input.
Transforms, filters, and plots the data as a horizontal stacked bar chart.
Saves the plot as a .png and .svg file.

- main()
Parses command-line arguments to get the directory path.
Calls collect_data_from_overviews and plot_stacked_bar functions to generate plots for 'DIAMOND', 'dbCAN_sub', and 'HMMER' data.
Command-Line Usage
To use the script from the command line, navigate to the script's directory and type:

```bash
python dbscan_analysis_graph_barplot.py /path/to/directory
 ```
Replace script_name.py with the actual name of the script and /path/to/directory with the path to the directory containing the overview files.

## Example Outputs
The script generates horizontal stacked bar plots with distinct colors for each category and saves them in the .png and .svg format. It primarily focuses on showcasing the distribution and counts of different subfamilies across various species, ensuring the visualizations are clear, concise, and interpretable.

## Script 17: concatenate_results_dbcan.py concatenate multiple dbsub.out files generated by run_dbscan 4 (https://github.com/linnabrown/run_dbcan)
## Overview
This Python script is built to concatenate multiple dbsub.out files, which are generated by run_dbcan4, a tool utilized to annotate sequences with CAZy enzymes. Each dbsub.out file contains annotations specific to a species. The script aggregates the data from all files into a single TSV (Tab Separated Values) file and additionally appends a column to denote the species name.

## Prerequisites
-  Python 3
Pandas library

## Input Data Format
The dbsub.out file should encompass the following columns:
dbCAN subfam    Subfam Composition    Subfam EC    Substrate    Profile Length    Gene ID    Gene Length    E Value    Profile Start    Profile End    Gene Start    Gene End    Coverage

## Running the Script

To use the script, simply execute it via the command line, providing the root directory (which contains the subdirectories with dbsub.out files) as an argument:
```bash
python3 concatenate_results_dbcan.py [directory_path_with_run_dbscan_output]
 ```
- Replace concatenate_results_dbcan.py with the actual name of the script file and [directory_path_with_run_dbscan_output] with the path to the root directory.

## Output
- The script outputs a single TSV file named final_output_dbsub.tsv, which contains concatenated data from all dbsub.out files and includes an additional column specifying the species name, extracted from the directory name.

## Script 18: concatenate_results_dbcan_diamond.py concatenate multiple diamond.out files generated by run_dbscan 4 (https://github.com/linnabrown/run_dbcan)

## Overview
This Python script is built to concatenate multiple diamond.out files, which are generated by run_dbcan4, a tool utilized to annotate sequences with CAZy enzymes. Each diamond.out file contains annotations specific to a species. The script aggregates the data from all files into a single TSV (Tab Separated Values) file and additionally appends a column to denote the species name.

## Prerequisites
-  Python 3
Pandas library

## Input Data Format
The diamond.out file should encompass the following columns:
Gene ID   CAZy ID   % Identical   Length   Mismatches   Gap Open   Gene Start   Gene End   CAZy Start   CAZy End   E Value   Bit Score   Specie

## Running the Script
To use the script, simply execute it via the command line, providing the root directory (which contains the subdirectories with dbsub.out files) as an argument:
```bash
python3 concatenate_results_dbcan_diamond.py [directory_path_with_run_dbscan_output]
 ```
## Output
- The script outputs a single TSV file named final_output_diamond.tsv, which contains concatenated data from all diamond.out files and includes an additional column specifying the species name, extracted from the directory name.

## Script 19: concatenate_results_dbcan_diamond_hmmer.py concatenate multiple hmmer.out files generated by run_dbscan 4 (https://github.com/linnabrown/run_dbcan)

## Overview
This Python script is built to concatenate multiple hmmer.out files, which are generated by run_dbcan4, a tool utilized to annotate sequences with CAZy enzymes. Each hmmer.out file contains annotations specific to a species. The script aggregates the data from all files into a single TSV (Tab Separated Values) file and additionally appends a column to denote the species name.

## Prerequisites
-  Python 3
Pandas library

## Input Data Format
The hmmer.out file should encompass the following columns:
HMM Profile   Profile Length   Gene ID   Gene Length   E Value   Profile Start   Profile End   Gene Start   Gene End   Coverage Specie   E Value   Bit Score   Specie

## Running the Script
To use the script, simply execute it via the command line, providing the root directory (which contains the subdirectories with dbsub.out files) as an argument:
```bash
python3 concatenate_results_dbcan_diamond_hmmer.py [directory_path_with_run_dbscan_output]
 ```

## Output
- The script outputs a single TSV file named final_output_hmmer.tsv, which contains concatenated data from all diamond.out files and includes an additional column specifying the species name, extracted from the directory name.

## Script 20: move_genome_files_from_list.py

## Move Matching Partial Genome Files names

This repository contains a Python script for moving Genome files with end sufix .fna based on partial names found in a text file. The script is useful for organizing files in directories, especially when dealing with large datasets or project files.

## Description

The `move_genome_files_from_list.py` script reads a text file containing partial file names. It then searches for files that match these partial names in a provided directory and moves these files to a new destination directory.

## Requirements

- Python 3.x

## Installation

No specific installation is required, as the script uses only standard Python libraries.

## Usage

To use the script, you need to provide three arguments:

1. The path to the source directory where the files are located.
2. The path to the `.txt` file containing the partial names of the files.
3. The name of the destination directory where the matching files will be moved.

Example usage:

```bash
python move_genome_files_from_list.py /path/to/source_directory /path/to/exclude_genome.txt genomes_missing_run
 ```

## Script 20: Get Species from Genome ID
# Genome Information Fetcher
This Python script automates the task of collecting genome information from a list of access IDs and generating a TSV (Tab-Separated Values) file with that data. It utilizes the NCBI's command-line tools `datasets` and `dataformat` to fetch the information.

## Requirements
- Python 3.x
- Access to shell/terminal
- NCBI's command-line tools `datasets` and `dataformat` installed and configured in your environment

## How to Use
To use this script, you will need a file containing a list of genome access IDs, one per line. Then, run the script in the following manner:

    ```
    python script.py <path_to_genome_list> <output_file_name>
    ```

## Example of Input File (genome_list) :
- GCF_000001405.39
- GCF_000002985.6
- GCF_000001635
    
## Functions

### `read_genome_ids(file_path)`
**Description**: Reads genome IDs from a file and returns them as a list.

**Parameters**:
- `file_path`: Path to the file containing the access IDs.

### `generate_tsv(genome_list, temp_file_name)`
**Description**: Generates a temporary TSV file with the specified genomes' information.

**Parameters**:
- `genome_list`: List of access IDs of the genomes.
- `temp_file_name`: Name of the temporary TSV file to be created.

### `finalize_output(temp_file_name, output_file_name)`
**Description**: Filters and finalizes the output file, removing the temporary file in the process.

**Parameters**:
- `temp_file_name`: Name of the temporary TSV file.
- `output_file_name`: Name of the final TSV file.

## Notes
This script is an example of how to automate the collection of genome data using NCBI tools. Ensure that the NCBI command-line tools `datasets` and `dataformat` are correctly installed and accessible in your execution environment.

## License
Include here the license under which the code is made available, for example, MIT, GPL, etc.

## Script 21: get_genome_id_and_taxonomic_rank_from_protein.id.py

Protein Accession to Genome ID and Taxonomic Ranking.
This script is designed to retrieve genome IDs and taxonomic rankings for a list of protein accessions. It utilizes NCBI's E-utilities for fetching the relevant data and processes it to generate two output files: one containing the genome IDs and the other containing taxonomic information.

## Prerequisites
- Before running this script, ensure you have the following installed and set up:
- Python 3.x
- Biopython
- Pandas library
- NCBI's EDirect utilities (for the `efetch` command)

Additionally, you will need an internet connection to access NCBI's databases.

## Installation
No installation is needed. Just ensure all prerequisites are met.

## Usage
python script.py <input_filename> <genome_output_filename> <taxonomic_output_filename>

To use this script, you need to provide it with an input file containing a list of protein accession numbers, one per line. Then, specify the names for the two output files: one for the genome IDs and the other for the taxonomic information.

- `<input_filename>`: The file containing a list of protein accession numbers.
- `<genome_output_filename>`: The file where the retrieved genome IDs will be saved.
- `<taxonomic_output_filename>`: The file where the taxonomic information will be saved.


## Output Format

- The genome IDs output file will contain two columns: Protein Accession and Genome Accession.
- The taxonomic information output file will contain four columns: Protein Accession, Genome Accession, Species, and Lineage.

## Important Notes

- Ensure you replace `'your_email@example.com'` in the script with your actual email address as required by NCBI's E-utilities policy.
- The script assumes you have the `efetch` utility accessible from your command line and that it's properly configured to work with NCBI's services.
- The script does not handle all possible errors and exceptions. Please use it as a starting point and customize it as needed for robust error handling and specific use cases.

## Contributing

Contributions to enhance the functionality, improve error handling, or extend the script's capabilities are welcome. Please feel free to fork the repository, make your changes, and submit a pull request.

## License

This script is provided "as is", without warranty of any kind. You are free to use, modify, and distribute it as you see fit.

## Contact

For any questions or suggestions, please open an issue on the GitHub repository page.


## Script 21.1: taxo_table_processor.py

This Python script processes a given CSV file containing taxonomic data, dividing the Lineage column into multiple taxonomic levels, and then aggregates the data at each level. It supports generating output files both with and without protein accession numbers.**

## Prerequisites
Before running this script, you need to have Python installed on your system along with the following Python libraries:
- pandas
- argparse
- These dependencies can be easily installed using Conda, a popular package and environment management system.

## Installation
First, ensure you have Conda installed on your system. If you do not have Conda installed, follow the instructions on the official Conda installation guide.

Once Conda is installed, you can create a new environment and install the required libraries using the following commands:

# Create a new Conda environment
```
conda create --name taxo_processor python=3.8
```
# Activate the environment
```
conda activate taxo_processor
```
# Install pandas
```
conda install -c anaconda pandas
```
# No need to install argparse as it is part of the Python standard library
# Usage

**To use the script, you must provide a path to a CSV file as an argument which is generated by Script 21: get_genome_id_and_taxonomic_rank_from_protein.id.py. The CSV file should have a column named Lineage containing taxonomic data separated by semicolons (;).**

```python taxo_table_processor.py <path_to_your_csv_file>
```
## Output
The script will generate multiple TSV files, each corresponding to a different taxonomic level. For the species level, it will generate two files: one with and one without protein accession numbers. The files are named according to the taxonomic level and whether they include protein accession numbers.

## Script 22: run_interproscan.py

This Python script automates the process of running InterProScan for sequence analysis, specific results filtering with `grep`, data extraction with `cut`, and sequence retrieval with `bedtools`.

## Features
- Creates a virtual environment for dependency isolation.
- Installs necessary dependencies, including Java JDK, Perl, and bedtools.
- Checks for the presence of the `interproscan.sh` script and downloads/unpacks it if not present.
- Runs InterProScan with specific parameters for analyzing `.fasta` files.
- Filters and processes results to extract specific sequences.

## Prerequisites
Before running this script, ensure you have `python3` installed on your system, as well as access to `sudo` for installing necessary packages.

## Setup and Execution
1. Clone this repository or download the script to your local system.
2. Open a terminal and navigate to the directory where the script is located.
3. Make the script executable (if necessary) with:
   ```bash
   chmod +x script_name.py
# Run the script with the following command:
```
python3 run_interproscan.py path_to_your_fasta_file.fasta
```
Replace run_interproscan.py with the actual file name of the script and path_to_your_fasta_file.fasta with the full path of the .fasta file you want to analyze.

# Important Notes
This script assumes you are using a Debian/Ubuntu-based system for package installation via apt.
Ensure your user has appropriate permissions to execute commands with sudo without manual interaction for password entry.
Troubleshooting
If you encounter any issues related to the execution of InterProScan or related tools, check the following:
The specified .fasta file exists and is accessible by the script.
All required components (Java, Perl, Python, bedtools) have been correctly installed.
The interproscan.sh script is in the expected directory and is executable.
For more information on InterProScan setup and options, visit the InterProScan GitHub page.

## Considerations
- **Location and Permissions**: This README assumes that the user has basic system permissions and directory navigation knowledge.
- **Customization**: You may need to adjust the path or specific parameters of InterProScan based on your setup and needs.
- **External Link**: I included a link to the InterProScan GitHub page for users who want more detailed information about the tool.

This README should be placed in the root of the repository where the script is located for easy access by other users who wish to use or contribute to the project.

## Contributing
Contributions to this project are welcome. Please fork the repository and submit a pull request with your suggested changes.

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Script 23: check_genomes.py  Genome File Matcher
## Overview
This Python script checks for the presence of genome file IDs from a list in a specified directory and moves matched files to a new directory. It generates a detailed report including which files were matched, unmatched, and exclusive to the directory. This tool is useful for managing large genomic datasets, allowing researchers to organize their files more effectively.

## Features
- Reads genome IDs from a specified text file.
- Checks for the presence of these IDs in a given directory.
- Moves matched genome files to a new user-defined directory.
- Outputs a tabular report listing matches and unmatched files.
- Provides summary statistics of the matching operation.

## Installation
- No installation is necessary, but you need Python installed on your machine. The script runs with Python version 3.6 or higher. It requires `pandas` and `tabulate` libraries to be installed. You can install these dependencies via pip if you do not have them:
```bash
pip install pandas tabulate
```
## Usage
To use this script, you need to provide four command-line arguments: the path to the genome list file, the directory to check for these genomes, the output file for the report, and the directory where matched files will be moved.
- Command-Line Syntax:
```
python check_genomes.py <path_to_genome_list_file> <directory_to_check> <output_report_file> <new_directory_for_matched_files>
```
## Example: 

```bash
python check_genomes.py /path/to/genomes.txt /data/genomes /results/match_report.txt /data/matched_genomes
```
# Genome File Matcher

## Overview
This Python script checks for the presence of genome file IDs from a list in a specified directory and moves matched files to a new directory. It generates a detailed report including which files were matched, unmatched, and exclusive to the directory. This tool is useful for managing large genomic datasets, allowing researchers to organize their files more effectively.

## Features
- Reads genome IDs from a specified text file.
- Checks for the presence of these IDs in a given directory.
- Moves matched genome files to a new user-defined directory.
- Outputs a tabular report listing matches and unmatched files.
- Provides summary statistics of the matching operation.

## Installation
No installation is necessary, but you need Python installed on your machine. The script runs with Python version 3.6 or higher. It requires `pandas` and `tabulate` libraries to be installed. You can install these dependencies via pip if you do not have them:

```bash
pip install pandas tabulate
```
## Usage
To use this script, you need to provide four command-line arguments: the path to the genome list file, the directory to check for these genomes, the output file for the report, and the directory where matched files will be moved.
## Output

*The script will generate a tabular report in the specified output file and print the table to the console. The report includes:*

- Genomes from the input file.
- Matched genomes in the directory.
- Genomes exclusive to the directory.

## It will also print summary statistics including the total number of matches and the number of files exclusively found in the specified directory.

## Script 24:  proteinHMM.py 

This Python script automates the process of running `hmmsearch` against a collection of FASTA files using multiple HMM models. It is designed to streamline the identification of protein domains within large genomic datasets.

## Prerequisites

Before running this script, ensure you have the following installed:
- Python (version 3.6 or later)
- `pandas` library
- HMMER suite (specifically `hmmsearch`)

## Installation

Clone this repository to your local machine using:
```bash
git clone https://github.com/CNP-CIIMAR/bioinformatics-utilities.git
```

cd bioinformatics-utilities
```
python hmm_search_pipeline.py <models_dir> <fastas_dir> <output_dir>
```
## Arguments: 

- models_dir: Directory containing the HMM models with .hmm extension.
  **Examples of files inside of models_dir: AMP_L.fasta  AMP_L_mafft.fasta  AMP_L_mafft.hmm  AMP_L_mafft.hmm.h3f  AMP_L_mafft.hmm.h3i  AMP_L_mafft.hmm.h3m  AMP_L_mafft.hmm.h3p** 
- fastas_dir: Directory containing the FASTA files with .faa extension.
  **Example of files inside of models_dir**: GCF_002608075.1_protein.faa GCF_002608225.1_protein.faa GCF_002760395.1_protein.faa
  The files of predicted proteomes were obtained using the API of NCBI-named datasets (https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/)  
- output_dir: Directory where the output tables (*.tsv) will be saved.

## Features
- Automatic Directory Creation: If the output directory does not exist, it will be created automatically.
- Batch Processing: The script processes all models and FASTA files in the specified directories.
- Output Management: For each model and FASTA file combination, a table with search results is generated, including a concatenated final table for all FASTA files per model.
- Robust Error Handling: Includes error checks and exception handling during the hmmsearch execution.
- Output Description
- The output files are tab-separated values (TSV) files, each named according to the model and FASTA file processed. Additionally, a concatenated TSV file for each model containing all results across processed FASTA files is created.

- Each output file includes the following columns organized in a multi-level header format:
- identifier: Contains target and query identifiers and descriptions.
- full_sequence: Metrics for the full sequence match.
- best_domain: Metrics for the best domain match.
- domain_number_estimation: Estimates of domain numbers in various categories.
- Specie: Extracted species information from the query description.

## Contributing
- Contributions to this project are welcome. Please fork the repository and submit a pull request with your enhancements.

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Script 25: split_names_fasta.py

This script processes a protein table and a multi-FASTA file to group and save protein sequences with similar descriptions into separate FASTA files.

## Features
Reads a table of protein entries and a multi-FASTA file.
Filters protein sequences based on entries in the table.
Groups protein sequences by 90% similarity in protein names.
Saves grouped sequences into separate FASTA files in a specified output directory.

## Requirements
Python 3.x
pandas
Biopython
Installation
To install the required packages, you can use pip:

```bash
pip install pandas biopython
```
##Usage
Run the script with the following command:

python split_names_fasta.py <table_file> <fasta_file> <output_dir>
table_file: Path to the input table file in CSV format.
fasta_file: Path to the input multi-FASTA file.
output_dir: Path to the output directory where the grouped FASTA files will be saved.
Example

```bash
python split_names_fasta.py proteins_table.csv proteins.fasta output_directory
```
This command processes the proteins_table.csv and proteins.fasta files, grouping sequences with similar protein names, and saves the results in the output_directory.

##Script Details

The script follows these steps:

Read the Table: Loads the protein table using pandas.
Read the Multi-FASTA File: Loads the sequences from the FASTA file using Biopython.
Filter Sequences: Filters the sequences based on the entries in the table.
Group by Similarity: Groups the sequences by 90% similarity in protein names.
Save to Output Directory: Saves the grouped sequences into separate FASTA files in the specified output directory.

## Example Table Format
The input table should be in CSV format with columns such as:

Entry,Reviewed,Entry Name,Protein names,Gene Names,Organism,Length,PubMed ID
A0A1W6GW32,reviewed,STPS1_SALMI,(-)-5-epieremophilene synthase STPS1 (EC 4.2.3.199) (Sesquiterpene synthase 1) (SmSTPS1),STPS1,Salvia miltiorrhiza (Chinese sage),546,28487717
...
## Notes
Ensure that the column Protein names in the table does not contain any missing values. The script currently ignores rows with missing values in this column.
The output files will have sanitized names to remove any non-alphanumeric characters from the protein names.

## Contributing
- Contributions to this project are welcome. Please fork the repository and submit a pull request with your enhancements.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Script 26: # genome_filter_download_scaffold_complete_gbff_by_metadata_table.py 

This repository contains a Python script to filter a genomic dataset and download specific genome assemblies based on the filtered results. The script reads a tab-delimited input file, filters the rows according to specified conditions, saves the filtered data to an output file, and then downloads genome assemblies for the filtered results.

## Requirements

- Python 3
- pandas library
- Anaconda (for the datasets command-line tool)

## Installation
1. Clone this repository:
   ```bash 
    git clone <repository_url>
    cd <repository_directory>
    ```
2. Install the required Python packages:
    ```sh
    pip install pandas
    ```
3. Ensure you have the `datasets` command-line tool installed via Anaconda:
    ```sh
    conda install -c bioconda ncbi-datasets-cli
    ```
## Usage
To use the script, run the following command:
```bash 
python script.py <input_file> <output_file> <output_directory>
```
```bash
python script.py input.tsv output.tsv genomes/
```
## Script Description
filter_table(input_file, output_file)
Reads the input file, filters the rows where:

## Assembly Level is either 'Scaffold' or 'Complete'
Lineage contains 'Cyanobacteriota', 'Pseudomonadota', 'Myxococcota', 'Actinomycetota', or 'Eukaryota'
## Saves the filtered results to the specified output file.

download_genomes(df, output_dir)
Downloads genome assemblies listed in the filtered dataframe:

## Creates the output directory if it doesn't exist
Downloads genome assemblies using the datasets command-line tool
Saves the downloaded files in the specified output directory

## Script 26: assign_colors_taxonomy_taxonomic_rank.py

Taxonomic Order Color Assignment.  This repository contains a Python script to read a taxonomic table - results of taxonomic rank python code, assign unique colors to each taxonomic order, and handle missing order values by using the closest available taxonomic rank or labeling them as "Unknown order".

## Features
Assigns a unique color to each taxonomic order.
Handles missing values in the 'order' column by using the closest available taxonomic rank or assigning "Unknown order".
Outputs a new table with columns 'Genome ID', 'order', and 'color'.
## Requirements:
- Python 3.x
- pandas

# Installation
# Clone the repository:
```sh
git clone https://github.com/yourusername/taxonomic-order-color-assignment.git
```
```sh
cd taxonomic-order-color-assignment
```

## Install the required Python packages:

```sh
pip install pandas
```
## Usage
Prepare your input CSV file with the appropriate taxonomic columns.

Run the script:
```sh
python3 assign_colors_taxonomy_taxonomic_rank.py path_to_your_input_file.csv
```
## The script will generate an output file named output_with_colors.csv in the same directory.

# Input Format

- The input CSV file should have the following columns (tab-separated):

| Genome ID                            | user_taxa                              | superkingdom | kingdom | superphylum | phylum          | subphylum | superclass | class         | subclass | superorder | order       | suborder | superfamily | family       | subfamily | genus |
|--------------------------------------|----------------------------------------|--------------|---------|-------------|-----------------|-----------|------------|---------------|----------|------------|-------------|----------|-------------|--------------|-----------|-------|

## Example

| Genome ID                            | user_taxa                              | superkingdom | kingdom | superphylum | phylum          | subphylum | superclass | class         | subclass | superorder | order       | suborder | superfamily | family       | subfamily | genus |
|--------------------------------------|----------------------------------------|--------------|---------|-------------|-----------------|-----------|------------|---------------|----------|------------|-------------|----------|-------------|--------------|-----------|-------|
| GCF_030382115.1_ASM3038211v1_genomic | Nostoc sp. GT001                       | Bacteria     | NA      | NA          | Cyanobacteriota | NA        | NA         | Cyanophyceae  | NA       | NA         | Nostocales  | NA       | NA          | Nostocaceae  | NA        | Nostoc |
| GCA_015206945.1_ASM1520694v1_genomic | Nostocales cyanobacterium LEGE 12452   | Bacteria     | NA      | NA          | Cyanobacteriota | NA        | NA         | Cyanophyceae  | NA       | NA         | Nostocales  | NA       | NA          | NA           | NA        | NA    |

## Output
The output CSV file output_with_colors.csv will contain the following columns:

Genome ID
order
color
Example:

## Output

The output CSV file `output_with_colors.csv` will contain the following columns:

| Genome ID                            | order       | color   |
|--------------------------------------|-------------|---------|
| GCF_030382115.1_ASM3038211v1_genomic | Nostocales  | #cd5c5c |
| GCA_015206945.1_ASM1520694v1_genomic | Nostocales  | #cd5c5c |

### Example

```csv
Genome ID,order,color
GCF_030382115.1_ASM3038211v1_genomic,Nostocales,#cd5c5c
GCA_015206945.1_ASM1520694v1_genomic,Nostocales,#cd5c5c


## Script 27: Enzyme Search Tool: enzymesearchtool.py

This tool enzymesearchtool.py is designed to search for specific enzyme names within .gbk files located in a directory and its subdirectories. It extracts relevant information and annotations, and saves the results to an output file.

```sh
enzymesearchtool.py
```
## Requirements
Python 3.6+
No additional libraries are required beyond the Python standard library.

## Usage
Command Line Interface
To run the script, use the following command:

```sh
python enzymesearchtool.py <diretorio> <arquivo_enzimas> <arquivo_saida>
```

- <diretorio>: The root directory containing all subdirectories with AntiSMASH results.
- <arquivo_enzimas>: A text file containing the list of enzyme names or annotations, one per line.
- <arquivo_saida>: The name of the output file where the results will be saved.

## Example

```sh
python script.py ./results lista.enzimas Resultado.tsv
```

## Functionality

-  Read Enzyme List
-  The script reads a text file containing enzyme names, specified by the user.
-  Extract Species Name
-  Extracts the species name from the file name using a regular expression pattern.

## Search for Enzymes

Iterates through all subdirectories and files in the specified directory.
Searches for enzyme names within .gbk files.
If an enzyme is found, it extracts the corresponding annotation.

## Save Results

Saves the search results to the specified output file.

## Script 28: get_genome_taxonomy_from_protein_id.py

Essential the download of program datasets from NCBI using this link: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/ and change the path of the program in the code get_genome_taxonomy_from_protein_id.py
Genome and Taxonomic Information Retrieval Script.

# Description

- This Python script get_genome_taxonomy_from_protein_id.py retrieves genome and taxonomic information for a list of protein accessions. It utilizes the BioPython library to interact with NCBI's databases and subprocess to execute command-line operations for efetch.
- The script reads a file containing protein accession numbers, fetches - the corresponding genome IDs, and then retrieves taxonomic information for each protein. The results are saved in two output files.

## Usage
# Command Line

```bash
python get_genome_taxonomy_from_protein_id.py <input_filename> <genome_output_filename> <taxonomic_output_filename>
```
## Arguments

- <input_filename>: Path to the input file containing protein accession numbers (one per line).
- <genome_output_filename>: Path to the output file where genome IDs will be saved.
- <taxonomic_output_filename>: Path to the output file where taxonomic information will be saved.

# Setup
## Install Dependencies:

## Ensure you have Python installed.

# Install the required Python libraries using pip:
```bash
pip install pandas biopython
```
Set Your Email for Entrez:

# Update the Entrez.email variable in the script with your email address. This is required by NCBI to identify the user.

## Functionality
# Fetch Genome Information:

The script reads protein accession numbers from the input file.
For each accession number, it uses efetch to fetch genome IDs (RefSeq or INSDC).

The results are saved in the specified genome output file.

# Fetch Taxonomic Information:
- The script reads the genome output file.
- For each protein accession, it fetches the species and taxonomic lineage using Entrez.
- The results are stored in a pandas DataFrame and saved as a tab-separated values (TSV) file in the specified taxonomic output file.

## Example

- Input File (input.txt)

- NP_000507.1
- NP_001123456.1

# Each line one accession number of protein

```bash
python get_genome_taxonomy_from_protein_id.py input.txt genome_output.txt taxonomic_output.txt
```
## Output Files

- genome_output.txt: Contains the protein accession and corresponding genome ID.


# taxonomic_output.txt: Contains the protein accession, genome ID, species, and taxonomic lineage.

Protein Accession  Genome Accession  Species  Lineage
NP_000507.1  NC_000001.11  Homo sapiens  Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini; Catarrhini; Hominidae; Homo
NP_001123456.1  NC_000002.12  Mus musculus  Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Rodentia; Sciurognathi; Muridae; Murinae; Mus


## Script 29: organize_protein_genome_metadata.py 

The Protein Genome Finder script reads a list of protein accession numbers from an input file, searches for matching records in a provided genome table, and creates an output table with the corresponding genome information. If a protein accession is not found in the genome table, the script notes "not genome found" in the output.

## Usage

```bash
python organize_protein_genome_metadata.py <input_file> <table_file> <output_file>
```

# Arguments
- <input_file>: Path to the input file containing protein accession numbers (one per line).
- <table_file>: Path to the tab-separated genome table file that contains genome information with "Protein Accession" as one of the columns.
- <output_file>: Path to the output file where the resulting table will be saved.

# Setup
# Install Dependencies:

# Ensure you have Python installed.

#Install the required Python library using pip:

```bash
pip install pandas
```
## Functionality
# Read Protein List:
# The script reads protein accession numbers from the input file.
# Create Output Table:

# For each protein in the list, it searches for a matching record in the genome table.

- If a match is found, the script appends the corresponding genome information to the output data.
- If no match is found, it appends "not genome found" to the output data.
- The results are stored in a pandas DataFrame and saved as a tab-separated values (TSV) file in the specified output file.

# Example

# Input File (input.txt)

- NP_000507.1
- NP_001123456.1

# Genome Table (`genome_table.tsv`)

| Protein Accession | Genome Accession | Species       | Lineage                                                                                                                                                             |
|-------------------|------------------|---------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| NP_000507.1       | NC_000001.11     | Homo sapiens  | Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini; Catarrhini; Hominidae; Homo          |

# Command
``` bash
python script.py input.txt genome_table.tsv output.txt
```
# Output File (output.txt)

| Protein Accession | Genome Accession   | Species         | Lineage                                                                                                                                                             |
|-------------------|--------------------|-----------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| NP_000507.1       | NC_000001.11       | Homo sapiens    | Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini; Catarrhini; Hominidae; Homo          |
| NP_001123456.1    | not genome found   | not genome found | not genome found                                                                                                                                                   |


## Script 30: filter_table_genome_quality_generate_plot.py

# Genome Quality Filter and Classifier

# Description

This Python script filters and classifies genomes based on their completeness and contamination metrics. It reads a genome table from a CSV file, filters the genomes according to specified thresholds, classifies them into quality categories, and plots the results. The plot is saved in multiple formats.

## Usage

# Command Line

```bash
python filter_table_genome_quality_generate_plot.py <input_file> <output_file> <base_filename> <completeness_threshold> <contamination_threshold>
```
## Arguments

- <input_file>: Path to the input CSV file containing the genome table.
- <output_file>: Path to the output CSV file where the filtered genome table will be saved.
- <base_filename>: Base name of the file for saving the plots.
- <completeness_threshold>: Minimum value of completeness.
- <contamination_threshold>: Maximum value of contamination.

## Setup
# Install Dependencies:
# Ensure you have Python installed.

Install the required Python libraries using pip:

```bash
pip install pandas matplotlib argparse
```
## Functionality

## Filter Genomes:

The script reads the genome table from the input CSV file.
It filters genomes based on the provided completeness and contamination thresholds.
The filtered table is saved to the specified output file.

## Classify Genomes:

# Genomes are classified into four categories:

- High-quality draft: >90% complete, <5% contamination
- Medium-quality draft: ≥50% complete, <10% contamination
- Low-quality drafts: <50% complete, <10% contamination

## Genomes meeting custom thresholds: ≥completeness_threshold% complete, ≤contamination_threshold% contamination

# Plot Genome Quality:

- The script creates a bar plot showing the number of genomes in each quality category.
- The plot is saved in PNG, SVG, and JPEG formats.

# Example

```bash
python filter_table_genome_quality_generate_plot.py genomes.csv filtered_genomes.csv genome_quality 90 5
```
- Input File (`genomes.csv`)

| Protein Accession | Genome Accession | Species       | Lineage      | Completeness | Contamination |
|-------------------|------------------|---------------|--------------|--------------|---------------|
| NP_000507.1       | NC_000001.11     | Homo sapiens  | ...          | 95           | 2             |
| NP_001123456.1    | NC_000002.12     | Mus musculus  | ...          | 85           | 1             |
| ...               | ...              | ...           | ...          | ...          | ...           |

- Output File (`filtered_genomes.csv`)

| Protein Accession | Genome Accession | Species       | Lineage      | Completeness | Contamination |
|-------------------|------------------|---------------|--------------|--------------|---------------|
| NP_000507.1       | NC_000001.11     | Homo sapiens  | ...          | 95           | 2             |
| ...               | ...              | ...           | ...          | ...          | ...           |

# Plot Output

The plot will be saved as genome_quality.png, genome_quality.svg, and genome_quality.jpeg.

## Script 31: get_specie_name_lineage_from_genome_id.py

This Python script reads a list of genome IDs from an input file, processes each ID to extract a specific prefix, retrieves detailed information about each genome using external commands, and saves the results in a TSV (Tab-Separated Values) file. The script handles errors gracefully and ensures only relevant lines are included in the final output.

## Usage

# Command Line

```bash
python get_specie_name_lineage_from_genome_id.py <genome_list_file> <output_file_name>
```
## Arguments

<genome_list_file>: Path to the input file containing genome IDs (one per line).
<output_file_name>: Path to the output TSV file where the detailed genome information will be saved.

## Setup

# Ensure Dependencies Are Met:

This script relies on external commands (datasets and dataformat) being available in your system's PATH. Make sure these commands are installed and accessible. The dataset and dataformat program can be obtained from this NCBI link: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/ 


## Ensure you have Python installed.

# Functionality

Read Genome IDs:

## The script reads genome IDs from the specified input file.

# Extract Prefix:

- For each genome ID, the script extracts the prefix up to the second _.

# Generate Intermediate List:

Creates an intermediate list of genome prefixes.

# Generate TSV:

- Uses the datasets command to fetch detailed information about each genome.
- Writes the output to a temporary TSV file.

# Finalize Output:

Copies the header from the temporary file to the final output file.
Filters relevant lines and appends them to the final output file using grep.
Removes the temporary file.

## Example

# Command

```bash
python get_specie_name_lineage_from_genome_id.py genome_ids.txt genome_info_output.tsv
```
# Input File (genome_ids.txt)

- GCA_000001405.15_GRCh38

- GCA_000001635.9_GRCm39

# Intermediate List (intermediate_list_YYYYMMDD_HHMMSS.txt)

- GCA_000001405.15
- GCA_000001635.9

# Temporary File (`temp_output_YYYYMMDD_HHMMSS.tsv`)

| Assembly Accession | Organism Common Name | Organism Name |
|--------------------|----------------------|---------------|
| GCA_000001405.15   | Human                | Homo sapiens  |
| GCA_000001635.9    | Mouse                | Mus musculus  |

# Output File (`genome_info_output.tsv`)

| Assembly Accession | Organism Common Name | Organism Name |
|--------------------|----------------------|---------------|
| GCA_000001405.15   | Human                | Homo sapiens  |
| GCA_000001635.9    | Mouse                | Mus musculus  |


## Script 32: extract_genome_id_main.py

# Genome ID Prefix Extractor

# Description

This Python script reads a list of genome IDs from an input file, extracts the part of each genome ID that precedes the second underscore (`_`), and writes the results to a new output file.

## Usage

# Command Line

```bash
python extract_genome_id_main.py <input_file> <output_file>
```
# Arguments

- <input_file>: Path to the input file containing genome IDs (one per line).
- <output_file>: Path to the output file where the extracted prefixes will be saved.

## Setup

Install Python:

# Ensure you have Python installed on your system.

Save the Script:

- Save the provided script to a file, for example extract_prefixes.py.

# Functionality

# Read Genome IDs:

- The script reads genome IDs from the specified input file.

# Extract Prefix:

For each genome ID, the script extracts the prefix up to the second underscore (_).
Write Prefixes to File:

The script writes the extracted prefixes to the specified output file.

# Example

```bash
python extract_genome_id_main.py genome_ids.txt prefixes_output.txt
```

## Input File (genome_ids.txt)

- GCF_030382115.1_ASM3038211v1_genomic
- GCA_015206945.1_ASM1520694v1_genomic
- GCF_002631755.1_ASM263175v1_genomic
- GCF_002246015.1_ASM224601v1_genomic
- GCF_017114965.1_ASM1711496v1_genomic
- GCF_017115065.1_ASM1711506v1_genomic
- GCA_003326195.1_ASM332619v1_genomic
- GCF_017115095.1_ASM1711509v1_genomic
- GCF_017114945.1_ASM1711494v1_genomic
- GCF_017114995.1_ASM1711499v1_genomic

## Output File (prefixes_output.txt)

- GCF_030382115.1
- GCA_015206945.1
- GCF_002631755.1
- GCF_002246015.1
- GCF_017114965.1
- GCF_017115065.1
- GCA_003326195.1
- GCF_017115095.1
- GCF_017114945.1
- GCF_017114995.1

## Contributing

- Contributions to this project are welcome. Please fork the repository and submit a pull request with your enhancements.


## Script 33: genome-downloader.py

This Python script allows you to download GBFF (GenBank Flat File) files of genomes from a list of accession IDs. It uses the `datasets` tool provided by NCBI to perform the download.

## Requirements
- Python 3.x
- The NCBI `datasets` tool must be installed in the Python environment you are using. Make sure the path to the `datasets` executable is correct in the script.

## Installation

1. Clone este repositório:

```bash
git clone https://github.com/seuusuario/genome-downloader.git
```
```bash
 cd genome-downloader
```

Install the necessary dependencies.
This script requires that the NCBI datasets tool is available in your Python environment. If you haven't installed it yet, you can do so via conda:

```bash
conda install -c bioconda ncbi-datasets-cli
```

- Make sure the path to the datasets executable is correct in the genome_downloader.py script. By default, the path is /home/mattoslmp/anaconda3/envs/biopython/bin/datasets, but you can adjust it as needed.

## Usage

- Prepare a text file containing the genome IDs you want to download. Each line of the file should contain a genome accession ID. Example of a genome_ids.txt file:

- GCF_000001405.39
- GCF_000002265.5
- GCF_000002865.1

Run the script by passing the file with the IDs and the output directory as parameters:

```bash
python genome_downloader.py genome_ids.txt output_directory/
```
Where:

- genome_ids.txt is the file containing the genome IDs.

output_directory/ is the directory where the GBFF files will be saved.

The script will download each genome listed in the file and save the .zip files in the specified directory. Each .zip will contain the corresponding GBFF file.

## Example of Usage
Here is a complete example of how to use the script:
```bash
python genome_downloader.py genome_ids.txt ./genomes/
```
This command will download the genomes listed in genome_ids.txt and save the .zip files containing the GBFFs in the ./genomes/ directory.

# Common Errors

Error: "Failed to download ...": This can happen if the genome ID is invalid or if there are internet connectivity issues. Check the ID and try again.

"datasets command not found": Make sure the datasets tool is installed and the path in the script is correct.

## Script 34:download_genome_ncbi_datasets.py  - "Dataset Download Script"

title: "Dataset Download Script"
description: >
  This script downloads genome datasets using the `datasets` command-line tool.
  It runs the download process in the background and logs the progress to a log file.

requirements:
  - "Python 3.x"
  - "`datasets` command-line tool installed and accessible in your system PATH."

file_structure:
  - "processo_datasets.log: A log file where the script records the process information."
  - "nohup_output.txt: A file that captures the output and error messages from the background process."

## usage:
# command: 

```bash
python download_genome_ncbi_datasets.py <input_file.txt> <output_file.zip>
```
  parameters:
    - input_file: "A text file containing the genome accession numbers you wish to download."
    - output_file: "The name of the zip file where the downloaded genome data will be stored."

how_it_works:
  - step: "**Logging Configuration**"
    description: >
      The script sets up logging to record the process's start time, level of severity,
      and messages in the `processo_datasets.log` file.
  - step: "**Executing Command in Background**"
    description: >
      The script constructs a command using the `datasets` tool to download genome data.
      It then runs this command in the background using `subprocess.Popen`,
      redirecting both standard output and error to `nohup_output.txt`.
  - step: "**Process Management**"
    description: >
      The script logs the process ID (PID) of the background command,
      allowing users to track the download process.
# command: 
```bash
    python download_genome_ncbi_datasets.py genome_accessions.txt downloaded_genomes.zip
```
  # description: 
- This command will download genome data based on the accession numbers listed in
- `genome_accessions.txt` and save it to `downloaded_genomes.zip`.

notes:
  - "Make sure that the `datasets` command-line tool is correctly installed and configured before running the script."
  - "Check the `processo_datasets.log` and `nohup_output.txt` files for details about the execution and any potential errors."

additional_tips:
  - "**Personalize o Nome do Script:**"
    description: >
      If your Python script has a specific name (e.g., `download_genomes.py`),
      replace `script.py` with the correct name in the usage examples.
  - "**Adicione Seções se Necessário:**"
    description: >
      Depending on your project's needs, you can add additional sections such as
      "Installation", "Contribution", "License", etc.
  - "**Verifique os Caminhos dos Arquivos:**"
    description: >
      Ensure that the paths to the log files (`processo_datasets.log` and
      `nohup_output.txt`) are correct and accessible in the environment where
      the script will be executed.

# Script 35: Create get_metadata_from_genome_id.py

This script processes genome assembly IDs to retrieve and enrich metadata using the NCBI Entrez system. 
The metadata includes lineage information, biome distribution, geographic location, and coordinates (latitude and longitude) for each genome. The output consists of two files:

1. A primary file with metadata for all input assembly IDs.
2. A filtered file containing only genomes with both `BiomeDistribution` and geographic coordinates.

## Prerequisites

- **Python Libraries**:
  - `ete3` for taxonomy data.
  - `Biopython` (`Bio`) for NCBI Entrez interactions.
  - `xml.etree.ElementTree` for XML parsing.

# To install the necessary libraries, run:
```bash
pip install ete3 biopython
```
# Usage
To run the script, use the following command in the terminal:

```bash
python get_metadata.py input_file.tsv output_file.tsv
```

# Usage
- Input File: The input file should be a tab-separated file containing genome assembly IDs with a header. Example format:

Assembly
GCA_000003745.2
GCA_000004075.3

# Running the Script: Run the script as follows:

```bash
python get_metadata.py <input_file_with_assembly_ids> <output_file>
```
# This command will produce:

- genomes_metadata.tsv: Main output file with metadata for all assembly IDs.
- filtered_genomes_metadata.tsv: Filtered output file with only assembly IDs where BiomeDistribution, Latitude, and Longitude are available.

## Features

- **Lineage Retrieval**: Uses the NCBI taxonomy database to fetch lineage based on organism taxonomic ID.
- **Biome Distribution and Location Extraction**: Identifies and categorizes biomes using keywords aligned with GOLD standards (e.g., terrestrial, marine, freshwater).
- **Latitude and Longitude Parsing**: Extracts geographic coordinates for samples with specified locations.
- **Output Files**:
  - Main output with all metadata.
  - Filtered output containing only rows with populated `BiomeDistribution`, `Latitude`, and `Longitude`.

# The script will print a summary at the end, showing:

- Total Genome IDs with BiomeDistribution filled.
- Total Genome IDs with Latitude and Longitude.
- Total Genome IDs in the filtered file with all required metadata.

## Main output with all metadata.

# Both output files will contain columns as shown below:

| Assembly Accession | Organism Name | Organism Common Name | Organism Tax ID | Lineage | Assembly Level | BioProject Accession | BioSample Accession | GC Percent | Total Sequence Length | Sequencing Technology | Release Date | Collection Date | BioSample Description | Location | BiomeDistribution | Latitude | Longitude |
|--------------------|---------------|----------------------|-----------------|---------|----------------|----------------------|---------------------|------------|-----------------------|------------------------|--------------|----------------|------------------------|----------|-------------------|----------|-----------|
| GCA_000003745.2    | Vitis vinifera | wine grape           | 29760           | root; cellular organisms; Eukaryota;... | Chromosome | PRJEA18785 | SAMEA2272750 | 34.5 | 485326422 |  | 2009-12-07 | | BioSample entry for genome collection GCA_000003745 | USA: Angelo Coast Range Reserve, CA | Soil | 39.74 | -123.63 |
| GCA_000004075.3    | Cucumis sativus | cucumber            | 3659            | root; cellular organisms; Eukaryota;... | Chromosome | PRJNA33619 | SAMN02953750 | 32.5 | 224801081 | PacBio RSII; PacBio Sequel; 10X Genomics; Hi-C; Illumina | 2019-11-15 | | Sample from Cucumis sativus | Peru: Oxygen minimum zone | Unknown |  |  |
| GCA_000004515.5    | Glycine max    | soybean             | 3847            | root; cellular organisms; Eukaryota;... | Chromosome | PRJNA19861 | SAMN00002965 | 34.5 | 978386919 | ABI 3739 | 2021-03-10 | | Glycine max cv. Williams 82 callus from plants grown in dark condition | Canada: Vancouver, Saanich Inlet | Marine | 48.36 | -123.3 |
| ...                | ...           | ...                  | ...             | ...     | ...            | ...                  | ...                 | ...        | ...                   | ...                    | ...          | ...            | ...                    | ...      | ...               | ...      | ...       |

# Notes
- NCBI Request Limits: The script adds pauses to comply with NCBI request limits. Consider running with an NCBI API key for higher request thresholds.
- Error Handling: If a genome assembly ID does not return any metadata, it is skipped, and an error message is logged.

## Script 36: genome_itol_table_update.py

## Genome IDs Table Updater

This script is designed to update a table based on a list of Genome IDs. It adds 'Label' and 'Color' columns to the records that match the provided Genome IDs.

## Features
Loads a list of Genome IDs from a text file.
Loads an input table containing a column with phylogenetic tree node IDs.
Adds 'Label' and 'Color' columns for matching genomes based on the provided IDs.
Saves the updated table to an output file.

- Requirements
- Python 3.6 or higher
- Python libraries:
- pandas
- argparse
- os
- sys

## Installation
Clone or download this repository.
Install the required Python packages:
```bash
pip install pandas
```

## Usage
To run the script, use the following command:

```bash
python update_itol_table.py -i <genome_ids.txt> -t <input_table.tsv> -o <output_table.tsv>
```


## Arguments:
-i, --input_ids: Path to the genome IDs text file (e.g., genome_ids.txt).
-t, --input_table: Path to the input table file (e.g., table.tsv).
-o, --output_table: Path to the output table file where the updated table will be saved (e.g., updated_table.tsv).

## Example:
```bash
python genome_itol_table_update.py -i genome_ids.txt -t input_table.tsv -o updated_table.tsv
```

## Description

The script performs the following tasks:

- Loads the genome IDs from the text file.
- Loads the input table and checks for the column 'Tree node ID'.
- Extracts and cleans the genome prefixes from the IDs.
- Matches the cleaned IDs with the Genome IDs provided in the input file.
- Adds 'Label' and 'Color' columns to the matched rows and updates the table.
Saves the updated table to the specified output file.
- Error Handling

The script checks if the input files exist and handles any errors during the file reading and writing processes.


## Script 37 name:mibig_downloader.py MIBIG Downloader  

description: >
  Um script Python para baixar, descompactar e organizar arquivos do repositório MIBIG (Biosynthetic Gene Cluster).
  O script extrai IDs de um arquivo de texto, baixa os arquivos correspondentes do MIBIG, 
  descompacta e organiza os arquivos .gbk em um diretório de saída.

version: "1.0.0"

dependencies:
  - Python >= 3.6
  - requests
  - zipfile
  - argparse

usage: |
  python mibig_downloader.py <arquivo_entrada> <diretorio_saida> [--url_base <url>] [--log_falhas <arquivo_log>]

parameters:
  - arquivo_entrada:
      description: Caminho para o arquivo de entrada (txt) contendo os IDs no formato BGCXXXXX.
  - diretorio_saida:
      description: Diretório onde os arquivos .zip serão salvos e descompactados.
  - url_base:
      description: URL base do repositório MIBIG. (opcional, padrão: https://mibig.secondarymetabolites.org/repository/)
  - log_falhas:
      description: Arquivo para registrar IDs que falharam no download. (opcional, padrão: falhas_download.txt)

example:
 ```bash
 python  python mibig_downloader.py lista_ids.txt ./saida/
 ```
  - command: python mibig_downloader.py lista_ids.txt ./saida/
    - description: 
    - Baixa e organiza os arquivos .gbk para os IDs especificados em lista_ids.txt no diretório ./saida/.

  - command: python mibig_downloader.py lista_ids.txt ./saida/ --url_base https://nova_url_do_mibig/repository/
  - description: 
      - Baixa e organiza os arquivos .gbk para os IDs especificados em lista_ids.txt no diretório ./saida/
      - usando uma nova URL base para o repositório MIBIG.

## Script 36:organize_big_slice.py

# organize_big_slice.py

## Overview
The `organize_big_slice.py` script is designed to organize antiSMASH result directories into datasets for BiG-SLiCE based on a specified taxonomic level, "Order." This allows researchers to streamline data organization for subsequent analysis in BiG-SLiCE, using either symbolic links or file copies. The script is customizable, enabling users to define specific taxonomy columns and indexes to suit various dataset formats.

## Requirements
To use this script, the following are required:
- Python 3
- Modules:
  - `os`
  - `pandas`
  - `shutil`
  - `argparse`
  - `collections.defaultdict`
  - `re`
  - `logging`

## Usage
To run this script, use the following syntax:

## Usage
To run this script, use the following syntax:
```bash
python3 organize_big_slice.py --bigslice_dir <BiG-SLiCE directory> --antismash_dir <antiSMASH directory> --taxonomy_table <taxonomy CSV path> --assembly_column <column name for Assembly Accession> --lineage_column <column name for Lineage> --order_index <Order index in lineage> [--use_symlinks] [--log_file <log file path>]

# Example:

 ```bash
python organize_big_slice.py --bigslice_dir /mnt/disk-sdc/Fasta/Table_organization/bigslice  --antismash_dir /mnt/disk-sdc/Fasta/Table_organization/Results_eukaryotes_not_plants/fungi/output_antismash --taxonomy_table /mnt/disk-sdc/Fasta/Table_organization/fungi_data.csv --assembly_column "Assembly Accession" --lineage_column "Lineage" --order_index 17
 ```
# After was created in the dir bigslice the files:
 
 ```bash
/mnt/disk-sdc/Fasta/Table_organization/bigslice$ ls
dataset_1  datasets.tsv  taxonomy
 ```
 ```plaintext
(base) mattoslmp@cnp:/mnt/disk-sdc/Fasta/Table_organization/bigslice$ more datasets.tsv
#Dataset name   Path to dataset folder  Path to taxonomy file               Description of the dataset
dataset_1       dataset_1               taxonomy/taxonomy_dataset_1.tsv     Dataset agrupado por Order: 1
 ```

 ```bash
dataset_1: dir with all subdir with antismash results:
genome_GCA_001593155.1  genome_GCA_023748825.1  genome_GCA_027602005.1  genome_GCA_910591825.1  genome_GCA_910592015.1
 ```
 ```bash
taxonomy$ ls
taxonomy_dataset_1.tsv
 ```

```plaintext
(base) mattoslmp@cnp:/mnt/disk-sdc/Fasta/Table_organization/bigslice/taxonomy$ more taxonomy_dataset_1.tsv
Genome folder name      Kingdom/Domain         Class               Order             Family              Genus            Species   Organism/Strain
genome_GCA_910592015.1/ Fungi incertae sedis   Glomeromycotina     Unknown           Diversisporales     Diversisporaceae Diversispora  Unknown
genome_GCA_027600645.1/ Fungi incertae sedis   Chytridiomycota     incertae sedis    Cladochytriales     Nowakowskiellaceae Nowakowskiella Unknown
```

### Arguments

- `--bigslice_dir`: Specifies the path to the BiG-SLiCE directory where datasets will be created.
- `--antismash_dir`: Specifies the path to the antiSMASH results directory.
- `--taxonomy_table`: Path to the taxonomy table in CSV format.
- `--assembly_column`: Column name in the taxonomy table containing Assembly Accession (default: "Assembly").
- `--lineage_column`: Column name in the taxonomy table containing Lineage (default: "Lineage").
- `--order_index`: Index (0-based) in the Lineage list where the taxonomic level "Order" is located.
- `--use_symlinks`: Optional flag that, if specified, uses symbolic links instead of copying files.
- `--log_file`: Path to the log file (default: "organize_big_slice.log").


# Functions
- parse_arguments()
- Parses and validates command-line arguments, ensuring all required inputs are provided.

setup_logging(log_file)
Sets up a logging system that outputs both to a log file and to the console, enabling monitoring of the process and tracking errors or warnings.

extract_assembly_accession(dir_name)
Extracts the Assembly Accession from a directory name using a regular expression. Examples include:

'Result_GCA_000007865.1.fasta_genomic.fna' → 'GCA_000007865.1'
'Result_GCA_910592015.1_Diversispora_eburnea_AZ414A_genomic.fna' → 'GCA_910592015.1'
get_order(lineage, order_index)
Extracts the taxonomic "Order" level from a lineage string based on a specified index. If the index is out of range or missing, it defaults to 'Unknown'.

# main()
Orchestrates the entire process, which includes:

# Parsing command-line arguments.
Setting up logging.
Loading and validating the taxonomy table.
Creating mappings between Assembly Accession and Order.
Organizing datasets by Order, either by copying or creating symbolic links for .gbk files.
Logging progress and errors to both the console and log file.
Writing organized dataset information to a datasets.tsv file.
Output
The script generates the following files and directories in the specified BiG-SLiCE directory:

# Organized Datasets: Datasets are organized by "Order" in subdirectories named dataset_X for each taxonomic group.

datasets.tsv: A summary file listing each dataset and its associated paths for easy reference in BiG-SLiCE.
Taxonomy Files: Taxonomy files for each dataset (e.g., taxonomy_dataset_X.tsv) that detail the taxonomic classification of each genome in the dataset.
Log File
The script logs all actions, including errors and warnings, in a specified log file and also displays them in the console. By default, it writes to organize_big_slice.log, but this can be customized with the --log_file argument.

# Example output 

| # Dataset name | Path to dataset folder           | Path to taxonomy file               | Description of the dataset        |
|----------------|----------------------------------|-------------------------------------|-----------------------------------|
| dataset_1      | /path/to/bigslice/dataset_1      | taxonomy/taxonomy_dataset_1.tsv     | Dataset grouped by Order: X       |
| dataset_2      | /path/to/bigslice/dataset_2      | taxonomy/taxonomy_dataset_2.tsv     | Dataset grouped by Order: Y       |
| ...            | ...                              | ...                                 | ...                               |

## After all the data is prepared, just run:

 ```bash
bigslice -i /mnt/disk-sdc/Fasta/Table_organization/bigslice --threshold 300 -t 16 bigslice_output
 ```

# OBS: Considering that the bigslice installation step was performed previously.

## Notes

Ensure the taxonomy CSV file includes the columns specified by --assembly_column and --lineage_column.
The order_index parameter should be adjusted according to the lineage structure in your data to correctly reference the "Order" level.
If symbolic links are used, the --use_symlinks flag should be included in the command line.
Troubleshooting
Missing Columns: If required columns are not found in the taxonomy table, the script will log an error and halt.
Invalid Path: The script verifies that the antiSMASH and taxonomy paths exist before proceeding. Any invalid paths will result in an error message.

## Script 37: download_ncbi_datasets.py 

# NCBI Datasets Downloader

![License](https://img.shields.io/badge/license-MIT-blue.svg)
![Python](https://img.shields.io/badge/python-3.6%2B-blue.svg)

## 📚 **Descrição**

O **NCBI Datasets Downloader** é um script Python projetado para facilitar o download automatizado de genomas, proteínas e CDS (Sequências Codificantes) do [NCBI](https://www.ncbi.nlm.nih.gov/) utilizando a ferramenta de linha de comando [`datasets`](https://www.ncbi.nlm.nih.gov/datasets/). O script permite que os usuários especifiquem quais tipos de dados desejam baixar, gerencie a descompactação dos arquivos baixados, renomeie os arquivos conforme necessário e organize os resultados em diretórios estruturados.

## 🌟 **Características**

- **Download Automatizado**: Baixe genomas, proteínas e CDS de múltiplos assemblies do NCBI de forma eficiente.
- **Inclusão Condicional**: Inclua tipos de dados (`genome`, `protein`, `cds`) apenas quando especificado pelo usuário.
- **Renomeação Inteligente**: Renomeia arquivos baixados para facilitar a identificação e organização.
- **Registro de Logs**: Mantém registros detalhados das operações realizadas, facilitando a resolução de problemas.
- **Relatório Final**: Gera um relatório resumindo os downloads bem-sucedidos e as falhas.
- **Fácil Integração**: Pode ser facilmente integrado em pipelines de análise genômica automatizados.

## 🛠️ **Requisitos**

- **Sistema Operacional**: Unix/Linux ou Windows com suporte ao comando `datasets`.
- **Python**: Versão 3.6 ou superior.
- **Ferramenta `datasets` do NCBI**: Deve estar instalada e acessível no sistema.

## ⚙️ **Instalação**

1. **Clone este repositório:**

   ```bash
   git clone https://github.com/seu-usuario/ncbi-datasets-downloader.git
   cd ncbi-datasets-downloader
   
## NCBI_Datasets_Downloader:
  Installation:
    Create_Virtual_Environment:
      Optional: true
      Steps:
        - command: python3 -m venv venv
        - command: source venv/bin/activate  # For Unix/Linux
        - command: venv\Scripts\activate     # For Windows
   ##  Install_Required_Dependencies:
      Details: "The script uses only standard Python libraries, so no additional installations are necessary."
    ## Install_NCBI_Datasets_Tool:
      # Instructions_Link: "https://www.ncbi.nlm.nih.gov/datasets/docs/v2"
      # Example_Installation:
        - Tool: Homebrew
        - Commands:
          - brew install ncbi-datasets
      # Verify_Installation:
        - Command: datasets --version

  # Usage:
```bash
    "download_ncbi_datasetsv2.py"
 ```
    # Syntax:
    ```bash
      - command: python3 download_ncbi_datasetsv2.py <genome_ids_file> <output_dir> [--include <types>] [--assembly-level <levels>]
      ```
    # Parameters:
      - genome_ids_file: "Path to the text file containing genome IDs (one per line)."
      - output_dir: "Directory where the files will be downloaded, unzipped, and organized."
    # Options:
      - include:
        - Types_Of_Data:
          - genome: "Downloads genomic sequences (.fna)."
          - protein: "Downloads protein sequences (.faa)."
          - cds: "Downloads coding sequences (.fna)."
        # Note: "By default, only genome is included. Specify other types only if necessary."
        # Examples:
          - Include_Only_Genome: "--include genome"
          - Include_Genome_And_Protein: "--include genome protein"
          - Include_All_Types: "--include genome protein cds"
      # assembly-level:
        # Description: "Limits downloads to specific assembly levels."
        # Options:
          - chromosome
          - complete
          - contig
          - scaffold
        # Format: "Separate multiple options with commas."
        # Example: "--assembly-level chromosome,complete"
        # Note: "Applicable only if genome is included."
    # Usage_Examples:
      - Download_Only_Genomes: "python3 download_ncbi_datasetsv2.py genome_ids.txt res"
      - Download_Genomes_And_Proteins: "python3 download_ncbi_datasetsv2.py genome_ids.txt res --include genome protein"
      - Download_Only_Proteins: "python3 download_ncbi_datasetsv2.py genome_ids.txt res --include protein"
      - Download_Genomes_With_Specific_Assembly_Levels: "python3 download_ncbi_datasetsv2.py genome_ids.txt res --include genome --assembly-level chromosome,complete"
      - Download_Genomes_Proteins_And_CDS: "python3 download_ncbi_datasetsv2.py genome_ids.txt res --include genome protein cds"

  # Directory_Structure_And_Files:
    # After_Execution: "The directory structure will look like this:"
    # Example_Structure:
      res/:
        genomes_obtained/:
          - "GCA_000236985.2_genomic.fna"
          - "GCA_000236985.2.faa"
          - "GCA_000236985.2_cds.fna"
        - "download_log.txt"
        - "download_report.txt"
        - "GCA_000236985.2.zip"
    # Details:
      genomes_obtained/: "Contains the renamed and organized .fna, .faa, and .cds.fna files."
      download_log.txt: "Detailed log of operations performed by the script."
      download_report.txt: "Summary of successful downloads and failures."
      genome_id.zip: "ZIP files downloaded from NCBI (can be removed after extraction, as needed)."

 # Genome_IDs_File_Format:
    Details: "Each line should contain a complete Assembly ID, possibly with suffixes like .fna, .faa, or .cds.fna."
    Example:
      - "GCA_000523875.1_ASM52387v1_genomic.fna"
      - "GCA_000001405.39_ASM140v39_genomic.faa"
      - "GCA_000022825.1_ASM22825v1_genomic.cds.fna"

  # Logs_And_Reports:
    # Operation_Log:
      File: "download_log.txt"
      Details: "Contains detailed records of operations performed, including downloads, unzips, file movements, and any errors encountered."
      Example_Log_Entry:
        - "2024-11-15 16:16:40,440 - INFO - Directory 'genomes_obtained' created or already exists: res/genomes_obtained"
        - "2024-11-15 16:16:40,440 - INFO - Data types included for download: genome"
        - "2024-11-15 16:16:40,440 - INFO - Starting download for: GCA_000236985.2 with option '--include genome'"
        - "2024-11-15 16:16:40,440 - INFO - Successfully downloaded: res/GCA_000236985.2.zip"
        - "2024-11-15 16:16:40,461 - INFO - Unzipped: res/GCA_000236985.2.zip"
        - "2024-11-15 16:16:40,461 - INFO - Files present in 'res/ncbi_dataset/data/GCA_000236985.2': ['GCA_000236985.2_ASM23698v2_genomic.fna']"
        - "2024-11-15 16:16:40,461 - INFO - Moving res/ncbi_dataset/data/GCA_000236985.2/GCA_000236985.2_ASM23698v2_genomic.fna to res/genomes_obtained/GCA_000236985.2_genomic.fna"
    # Final_Report:
      # File: "download_report.txt"
      # Details: "Provides a summary of the downloads performed, highlighting which were successful and which failed."
      # Example_Report:
        Total_Items_Requested: 3
        Successfully_Downloaded: 2
        Failed_Downloads: 1
        Successful_Items:
          - "GCA_000523875.1"
          - "GCA_000001405.39"
        Failed_Items:
          - "GCA_000022825.1_ASM22825v1_genomic.cds.fna"

  # Troubleshooting:
    Problem: "File Not Found After Extraction"
    Error_Message: "File '<genome_id>_genomic.fna' not found after extraction in: <path>"
    Possible_Solutions:
      - Check_Contents_Of_ZIP:
          Command: "unzip -l <zip_file>"
      - Verify_File_Naming:
          Details: "Ensure that the search patterns in the script match the actual names of the files inside the ZIP."
      - Check_Directory_Permissions:
          Commands:
            - "ls -ld <directory>"
            - "ls -l /usr/local/bin/datasets"
      - Re-run_Script:
          Details: "If the download or extraction failed, try re-running the script for that specific assembly_id."
    # Other_Common_Errors:
      Error: "Download Failure"
      Solutions:
        - "Check the internet connection."
        - "Confirm that the assembly_id is correct and available on NCBI."
        - "Ensure that the datasets tool is up to date."

  # Security_Considerations:
    File_Permissions: "Ensure that the directories and files used by the script have appropriate permissions to prevent unauthorized access."
    Input_Validation: "The script performs basic checks on assembly_ids, but ensure the IDs file is correctly formatted to avoid unexpected behaviors."

  # Contributions:
    Details: "Contributions are welcome! Feel free to open issues or submit pull requests to improve this project."
    Steps:
      - Fork_Repository: "Fork this repository."
      - Create_Branch: 
          Command: "git checkout -b feature/new-feature"
      - Commit_Changes:
          Command: "git commit -m 'Add new feature'"
      - Push_To_Branch:
          Command: "git push origin feature/new-feature"
      - Open_Pull_Request: "Submit a pull request on GitHub."

## Script 38: analise_CAL_domain_report_copy_gbk.py # Genome CAL Domain Analysis Script

## Description

This Python script automates the analysis of GenBank (`.gbk`) files within multiple subdirectories. It identifies specific protein domains (`CAL_domain` and optionally `AMP-binding`) within genomic features, calculates the total size of sequences, and organizes the results into a summary report. Additionally, it filters and copies relevant `.gbk` files to a designated directory for further analysis.

## Key Features

- **Domain Detection**: Identifies and counts occurrences of specified domains (`CAL_domain` and `AMP-binding`) within genomic features.
- **Genome ID Extraction**: Extracts genome IDs from subdirectory names following a specific naming convention.
- **Size Calculation**: Computes the total size of genomic sequences in Megabases (Mb) and Gigabases (Gb).
- **Filtered File Copying**: Copies `.gbk` files containing identified domains to a filtered directory with prefixed Genome IDs.
- **Comprehensive Logging**: Generates detailed logs of the processing steps and outcomes.
- **Summary Report**: Compiles a `summary.csv` file with aggregated data for easy reference.

## Requirements

- **Python**: Version 3.6 or higher.
- **Python Libraries**:
  - `Biopython`
  - `argparse`
  - `pathlib`
  - `csv`
  - `shutil`
  - `os`

## Installation

1. **Clone the Repository**:

    ```bash
    git clone https://github.com/your-username/your-repository.git
    cd your-repository
    ```

2. **Create a Virtual Environment (Optional but Recommended)**:

    ```bash
    python3 -m venv venv
    source venv/bin/activate  # For Unix/Linux
    venv\Scripts\activate     # For Windows
    ```

3. **Install Dependencies**:

    ```bash
    pip install biopython
    ```

## Usage

Run the script via the command line, providing the input directory containing `.gbk` files and specifying the log file location. Optional parameters allow for additional domain searches.

### Syntax

```bash
python analise_CAL_domain_report_copy_gbk.py <input_dir> <log_file> [--search-amp-binding]
```
# Arguments

<input_dir>: Required. Path to the input directory containing subdirectories with .gbk files.
<log_file>: Required. Path to the log file where processing details will be saved.
--search-amp-binding: Optional. If specified, the script will also search for the AMP-binding domain.

#Example Command
# Generated Outputs
```bash
python analise_CAL_domain_report_copy_gbk.py \
    /path/to/input_directory \
    /path/to/log_file.log \
    --search-amp-binding
```

# Log File (<log_file>.log):

- Contains detailed information about the processing of each subdirectory and .gbk file.
- Logs include domain counts, file sizes, and any errors encountered during processing.
- Summary Report (summary.csv):
- summary.tsv Copied to a filtrados_subdir_CAL directory located alongside the input directory.
- Files are prefixed with their respective Genome IDs for easy identification.

# Print log output (Example).

# Há espaço suficiente no disco para copiar os arquivos.

- Processamento concluído.
- Total de subdiretórios acessados: 25
- Total de subdiretórios com 'CAL_domain': 5
- Tamanho total de todos os arquivos .gbk encontrados: 0.27 Mb / 0.00 Gb
- Resumo salvo em: fungi/summary.csv
- Log detalhado salvo em: log_file

**Arquivos gbk filtrados copiados para o diretório: filtrados_subdir_CAL** 

# Located within the input directory.
summary.tsv:

| # Columns:
|------------------------------------------
| # Relatório de Correspondências para Bacteria
| - Assembly,CAL_domain,AMP-binding,Total_size_MB,Total_size_GB
| - GCA_910592015.1,0,0,0.00 Mb,0.00 Gb
| - GCA_027600645.1,0,0,0.00 Mb,0.00 Gb
| - GCA_021556665.1,0,0,0.00 Mb,0.00 Gb
| - GCA_910592205.1,0,0,0.00 Mb,0.00 Gb
| - GCA_910591965.1,0,0,0.00 Mb,0.00 Gb

## Script 39: merge_summary_metadata.py  Para executar esse codigo primeiro é necessario executar o codigo anterior: analise_CAL_domain_report_copy_gbk.py

# Merge Summary and Metadata Script

## Description

This Python script automates the process of merging two data files — `summary.csv` and `Supplementary_table2.tsv` — based on the common column `Assembly`. Additionally, it generates a detailed report on the matches found, specifically focusing on entries classified under the taxonomic category "Bacteria". The report includes counts of complete and incomplete matches, as well as a list of genome IDs (`Assembly`) that have a `CAL_domain` value of 0.

## Key Features

- **Data Merging**: Combines information from `summary.csv` with metadata from `Supplementary_table2.tsv`, retaining all rows from the metadata file.
- **Taxonomic Filtering**: Filters and analyzes data specific to entries classified as "Bacteria".
- **Detailed Reporting**:
  - Total number of "Bacteria" entries.
  - Number of complete matches (`Assembly`, `CAL_domain`, and optionally `AMP-binding`).
  - Number of incomplete matches.
  - List of `Assembly` IDs with `CAL_domain` equal to 0.

## Requirements

- **Python**: Version 3.6 or higher.
- **Python Libraries**:
  - `pandas`
  - `argparse`
  - `os`

## Installation

1. **Clone the Repository**:

   ```bash
   git clone https://github.com/your-username/your-repository.git
   cd your-repository
 ```
merge_summary_metadata.py
 ```
## Usage

Run the script via the command line, providing paths to the input and output files, along with optional parameters as needed.

 ```bash
python merge_summary_metadata.py --summary <path_to_summary.csv> --metadata <path_to_Supplementary_table2.tsv> --output <path_to_merged_summary.csv> [--search-amp-binding] [--report <path_to_report.txt>]
```

# Arguments
--summary: Required. Full path to the summary.csv file.
--metadata: Required. Full path to the Supplementary_table2.tsv file.
--output: Required. Path where the merged merged_summary.csv will be saved.
--search-amp-binding: Optional. If specified, includes the AMP-binding column in the final report.
--report: Optional. Path where the correspondence report will be saved.

OBS: metadata table can be obtained with the script download_ncbi_datasetsv2.py from this same repository.

## Example: 
 ```bash
python merge_summary_metadata.py \
    --summary /path/to/summary.csv \
    --metadata /path/to/Supplementary_table2.tsv \
    --output /path/to/merged_summary.csv \
    --search-amp-binding \
    --report /path/to/report_bacteria.txt
 ```
# Generated Outputs
 Merged File (merged_summary.csv):

- Contains all rows from the metadata file Supplementary_table2.tsv with data from summary.csv appended where available.
- Metadata columns are positioned after the Assembly column, followed by summary columns (CAL_domain, AMP-binding, Total_size_MB, Total_size_GB).
- If --search-amp-binding is not specified, the AMP-binding column will be excluded.

# Report (report_bacteria.txt) (if specified):
# Summary Counts:

- Total number of entries with Lineage containing "Bacteria".
- Number of complete correspondences (Assembly, CAL_domain, and AMP-binding if specified).
- Number of incomplete correspondences.
- List of Assembly IDs with CAL_domain equal to 0:
- Provides a detailed list of genomes with CAL_domain set to 0, facilitating the identification of specific cases that require further attention.


| # Relatório de Correspondências para Bacteria
|------------------------------------------
| - Total de linhas com 'Lineage' contendo 'Bacteria': 250
| - Quantidade de correspondências com informações completas ('Assembly', 'CAL_domain', 'AMP-binding'): 200
| - Quantidade de correspondências com informações incompletas: 50 

| # Relatório de Correspondências para Bacteria
|------------------------------------------
| - GCA_021323495.1
| - GCA_012345678.9
| - GCA_987654321.0


## License

All codes in this project is licensed under the MIT License. Leandro de Mattos Pereira built all the codes in this repository. 

contact: mattoslmp@gmail.com

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## More about:

[Junior Researcher, Leandro de Mattos Pereira](https://mattoslmp.github.io)

[CNP team, Pedro Leão, Researcher Leader](https://leaolab.wixsite.com/leaolab)
