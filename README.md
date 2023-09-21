

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


# Script 2: Get Genomes ID from List of Proteins

This script get_genomes_id_by_proteina_id.py retrieves genomes IDs from a list of protein accessions. It uses the NCBI `efetch` command to fetch information about each protein accession from the NCBI protein database.

## Usage:

1. Make sure you have Python 3 installed on your system.
2. Install the required dependencies:
- Python 3
- openpyxl

## How to Run:
- Make sure you have Python 3 installed.
- Open a terminal or command prompt.
- Navigate to the directory where the script is located.
- Run the following command:

```shell
python get_taxonomic_rank_from_protein.id.py <input_filename> <output_filename>
```

**Replace <input_filename> with the path to the file containing protein accessions and <output_filename> with the desired path to save the results.**

**Make sure you have the necessary permissions to read the input file and write to the output file.**

## Output Format

- The output file will be in a tabular format, with each line containing a protein accession and its corresponding genome ID, separated by a tab with informations about the Assembly and each protein
---
Note: The script assumes the existence of the `efetch` command and that it is available in your environment. It also relies on the openpyxl module, which may need to be installed if not already available in your Python environment.

If you need any help or have questions, feel free to open an issue in the GitHub repository.

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


## Script 8: HMMER Table Converter: convert_hmm_output_to_table.py

**This is a Python script for converting HMMER output tables into a more user-friendly format. It reads HMMER output files in either tblout or domtblout format and converts them into a tab-separated values (TSV) file.

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
python script.py <input_path> <output_path>
 ```
Replace <input_path> with the path to your HMMER output file, and <output_path> with the desired path for the converted TSV file.

**Example usage:**

 ```shell
python script.py input.tblout output.tsv
 ```
The script will convert the input file into a TSV file at the specified output path.

##Supported Formats
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
## Notes

The script was developed using HMMER v3.2.1.

The assert_acceptable_arguments function is used to validate the input arguments. If an invalid argument is provided for program or format, a ValueError will be raised.

**The script utilizes the pandas library for data manipulation and provides a convenient way to extract species information from the query_description column using regular expressions.**


# Script 9: get_taxonomic_rank_specie.py:

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
# Script 10: filter_dup_gca_gcf_keep_gcf.py

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

# Script 10:

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
## License

- This project is licensed under the MIT License. See the LICENSE file for more information.

**If the Bioinformatics scripts utilities were useful, please give proper credit to the authors or include the following citation: "CNP-CIIMAR GitHub [https://github.com/CNP-CIIMAR]."**

## More about:

[Junior Researcher, Leandro de Mattos Pereira](https://mattoslmp.github.io)

[CNP team, Dr. Pedro Leão, Researcher Leader](https://leaolab.wixsite.com/leaolab)





