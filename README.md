

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

Input file format
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

- Requirements
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

Installation
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
Prerequisites
Python 3.x
requests library (for the InterProScan API)
MAFFT
HMMER
InterProScan (optional, if local execution is desired)

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
## This script maps RefSeq Protein IDs to UniProtKB IDs using the UniProt REST API. It fetches various attributes for each protein, such as entry type, description, lineage, and other relevant data.

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
## The script fetches and outputs the following columns for each RefSeq Protein ID from UniprotkB:

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

# Protein Accession to Genome ID and Taxonomic Ranking

This script is designed to retrieve genome IDs and taxonomic rankings for a list of protein accessions. It utilizes NCBI's E-utilities for fetching the relevant data and processes it to generate two output files: one containing the genome IDs and the other containing taxonomic information.

## Prerequisites

Before running this script, ensure you have the following installed and set up:

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

### Taxonomic Table Processor

**This Python script processes a given CSV file containing taxonomic data, dividing the Lineage column into multiple taxonomic levels, and then aggregates the data at each level. It supports generating output files both with and without protein accession numbers.**

## Prerequisites
Before running this script, you need to have Python installed on your system along with the following Python libraries:

pandas
argparse
These dependencies can be easily installed using Conda, a popular package and environment management system.

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


## Usage

**To use the script, you must provide a path to a CSV file as an argument which is generated by Script 21: get_genome_id_and_taxonomic_rank_from_protein.id.py. The CSV file should have a column named Lineage containing taxonomic data separated by semicolons (;).**

```
python taxo_table_processor.py <path_to_your_csv_file>
```

## Output

The script will generate multiple TSV files, each corresponding to a different taxonomic level. For the species level, it will generate two files: one with and one without protein accession numbers. The files are named according to the taxonomic level and whether they include protein accession numbers.

## Script 22: Run_interproscan.py
# InterProScan Automation Script

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
python3 script_name.py path_to_your_fasta_file.fasta
```
Replace script_name.py with the actual file name of the script and path_to_your_fasta_file.fasta with the full path of the .fasta file you want to analyze.

Important Notes
This script assumes you are using a Debian/Ubuntu-based system for package installation via apt.
Ensure your user has appropriate permissions to execute commands with sudo without manual interaction for password entry.
Troubleshooting
If you encounter any issues related to the execution of InterProScan or related tools, check the following:

The specified .fasta file exists and is accessible by the script.
All required components (Java, Perl, Python, bedtools) have been correctly installed.
The interproscan.sh script is in the expected directory and is executable.
For more information on InterProScan setup and options, visit the InterProScan GitHub page.


### Considerations

- **Location and Permissions**: This README assumes that the user has basic system permissions and directory navigation knowledge.
- **Customization**: You may need to adjust the path or specific parameters of InterProScan based on your setup and needs.
- **External Link**: I included a link to the InterProScan GitHub page for users who want more detailed information about the tool.

This README should be placed in the root of the repository where the script is located for easy access by other users who wish to use or contribute to the project.

## Contributing
Contributions to this project are welcome. Please fork the repository and submit a pull request with your suggested changes.

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## More about:
[Junior Researcher, Leandro de Mattos Pereira](https://mattoslmp.github.io)

[CNP team, Dr. Pedro Leão, Researcher Leader](https://leaolab.wixsite.com/leaolab)


