import os
import sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from sklearn.manifold import TSNE
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO
from scipy.cluster.hierarchy import linkage, dendrogram
from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster import KMeans
from sklearn.metrics.pairwise import cosine_similarity
from gensim.models import Word2Vec
from sklearn.metrics import silhouette_score
import seaborn as sns
from transformers import BertTokenizer, BertModel
import torch

# Global variables
embedding = None
cluster_labels = None
fasta_ids = []
table_data = []

def save_sequences_and_print_info(alignment, cluster_labels, relevant_info):
    # Save sequences from each cluster into separate FASTA files
    for i, record in enumerate(alignment):
        cluster_label = cluster_labels[i]
        fasta_file = f'cluster_{cluster_label}.fasta'
        with open(fasta_file, 'a') as file:
            file.write(f'>{record.id}\n{record.seq}\n')

    # Print information about the sequence clusters
    for i, record in enumerate(alignment):
        cluster_label = cluster_labels[i]
        protein_accession = record.id
        protein_code = protein_accession.split()[0]  # Get only the code from Protein.accession
        matching_rows = relevant_info[relevant_info['Protein.accession'].str.startswith(protein_code)]
        if not matching_rows.empty:
            fatty_acid_length = matching_rows['Fatty acids length'].values[0]
            fatty_acid_chain = matching_rows['Fatty acids chain'].values[0]
            specie = matching_rows['Specie'].values[0]
            phylum = matching_rows['Phylum'].values[0]
            print(f'Protein.accession: {protein_accession} | Fatty acids length: {fatty_acid_length} | Fatty acids chain: {fatty_acid_chain} | Specie: {specie} | Phylum: {phylum} | Cluster: {cluster_label}')
        else:
            print(f'Protein.accession: {protein_accession} | No matching information found | Cluster: {cluster_label}')

# Function to perform BERT embedding on k-mers
def bert_embedding(kmers, model_name):
    tokenizer = BertTokenizer.from_pretrained(model_name)
    model = BertModel.from_pretrained(model_name)

    embeddings = []
    for kmer in kmers:
        inputs = tokenizer(kmer, return_tensors="pt", padding=True, truncation=True)
        with torch.no_grad():
            outputs = model(**inputs)
        embeddings.append(outputs.last_hidden_state.mean(dim=1).numpy())

    return embeddings

def run_mafft(input_file, table_file):
    global embedding, cluster_labels, fasta_ids, table_data
    table_data = pd.read_csv(table_file, delimiter="\t")

    # Obtain relevant columns from the table
    fatty_acids_lengths = table_data["Fatty acids length"]
    species = table_data["Specie"]
    phylum = table_data["Phylum"]
    protein_accession = table_data["Protein.accession"]

    # Execute MAFFT
    mafft_cline = MafftCommandline(input=input_file, thread=40, maxiterate=1000, localpair=True)
    stdout, stderr = mafft_cline()

    with open("alignment.fasta", "w") as handle:
        handle.write(stdout)

    # Read the alignment file
    alignment = AlignIO.read("alignment.fasta", 'fasta')
    
    # Generate k-mers
    
    k = 3
    kmer_groups = {}
    concatenated_kmers = []  # Initialize the list of all k-mers
#    kmer_info = []  # Initialize the list of k-mer information
    
    
    # looping through each record of the sequence alignment, and for each record, generating all possible k-mers of the sequence.
    #These k-mers are then added to the concatenated_kmers list, which contains all k-mers from all sequences
      # Agrupar os k-mers de acordo com o comprimento do Ã¡cido graxo (fatty acid)  
    for record in alignment:
        sequence = str(record.seq)
        seq_len = len(sequence)
        protein_accession_alignment = record.id.split()[0]  # Extract the first column of Protein.accession from the alignment
        matching_rows = table_data['Protein.accession'].str.split().str[0] == protein_accession_alignment
        matching_info = table_data[matching_rows]
        if not matching_info.empty:
            protein_accession_table = matching_info['Protein.accession'].values[0]
            fatty_acid_length = matching_info['Fatty acids length'].values[0]
            if fatty_acid_length not in kmer_groups:
                kmer_groups[fatty_acid_length] = []
            kmers = [sequence[i:i + k] for i in range(seq_len - k + 1)]
           # concatenated_kmers.extend(kmers)  # Add the k-mers to the list of all k-mers
         #   kmer_info.extend([(protein_accession_table, species, fatty_acid_length)] * len(kmers)) #Not needed # Add the sequence information for each k-mer
            kmer_groups[fatty_acid_length].extend(kmers)
            kmer_groups[fatty_acid_length].append((protein_accession_table, kmers))
        else:
            print(f"No matching rows found for Protein.accession: {protein_accession_alignment}")



    # Realizar o embedding BERT para cada grupo de k-mers por comprimento de ácido graxo
    bert_models = {}
    concatenated_embeddings = []
    for fatty_acid_length, kmers in kmer_groups.items():
        model_name = "Rostlab/prot_bert_bfd"
        embeddings = bert_embedding(kmers, model_name)  # Criação dos embeddings BERT
        bert_models[fatty_acid_length] = embeddings
        concatenated_embeddings.extend(embeddings)

        # Print k-mers for each group
    for fatty_acid_length, group in kmer_groups.items():
        print(f"Fatty Acid Length: {fatty_acid_length}")
        for protein_accession, kmers in group:
            print(f"Protein Accession: {protein_accession}")
            print(f"K-mers: {kmers}")
            print()	

   #Perform K-means clustering on the concatenated model of fatty acid groups.
    k = 3
    n_init = 10

    kmeans = KMeans(n_clusters=k, random_state=0)
    cluster_labels = kmeans.fit_predict(concatenated_embeddings)

    # Map cluster labels back to k-mers
    clustered_kmers = {i: [] for i in range(k)}
    for kmer, label in zip(concatenated_kmers, cluster_labels):
        clustered_kmers[label].append(kmer)    

    # Calculate cosine similarity between embeddings
    similarity_scores = {}

    for fatty_acid_length_1, embeddings_1 in bert_embeddings.items():
        similarity_scores[fatty_acid_length_1] = {}
        for fatty_acid_length_2, embeddings_2 in bert_embeddings.items():
            if fatty_acid_length_1 == fatty_acid_length_2:
                continue
            similarity_matrix = cosine_similarity(embeddings_1, embeddings_2)
            similarity_score = similarity_matrix.mean()
            similarity_scores[fatty_acid_length_1][fatty_acid_length_2] = similarity_score

    # Print similarity scores
    for fatty_acid_length_1, similarity_dict in similarity_scores.items():
        for fatty_acid_length_2, similarity_score in similarity_dict.items():
            print(f"Similarity between fatty acid lengths {fatty_acid_length_1} and {fatty_acid_length_2}: {similarity_score}\n")

   
# Crie um mapeamento de cores para os comprimentos de ácidos graxos
    color_map = {
            'Short fatty acids': 'red',
            'Medium fatty acids': 'green',
            'Long fatty acids': 'blue'
            }
# Crie um mapeamento de marcadores para os rótulos de cluster
    marker_map = {0: 'o', 1: 's', 2: '^'}

# Execute o t-SNE e a plotagem
# Perform t-SNE dimensionality reduction
    tsne = TSNE(n_components=3)
    tsne_embeddings = tsne.fit_transform(concatenated_embeddings)

# Plot the results with ellipses for cluster separation
    fig = plt.figure(figsize=(16, 12))
    ax_tsne = fig.add_subplot(111, projection='3d')

    for fatty_acid_length, embedding, cluster_label in zip(fatty_acids_lengths, tsne_embeddings, cluster_labels):
    # Use a different color for each fatty acid length
        color = color_map.get(fatty_acid_length, 'black')  # Use black as the default color

    # Use a different marker for each cluster label
        marker = marker_map.get(cluster_label, 'o')  # Use 'o' as the default marker

        ax_tsne.scatter(
            embedding[0],
            embedding[1],
            embedding[2],
            color=color,
            marker=marker,
            s=50,
            alpha=0.8,
            label=f"Cluster {cluster_label} Fatty Acid Length {fatty_acid_length}",
        )

    ax_tsne.set_xlabel('Component 1', fontsize=12)
    ax_tsne.set_ylabel('Component 2', fontsize=12)
    ax_tsne.set_zlabel('Component 3', fontsize=12)

# Add a legend to the plot
    ax_tsne.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), fancybox=True, shadow=True, ncol=5, fontsize=10)

# Calculate silhouette score
    silhouette_avg = silhouette_score(concatenated_embeddings, cluster_labels)
    print(f'Silhouette Score: {silhouette_avg:.4f}')

# Add the silhouette score to the right of the plot
    plt.annotate(
        f'Silhouette Score: {silhouette_avg:.4f}',  # Limit to 4 decimal places
        xy=(1, 0.5), 
        xycoords='axes fraction', 
        textcoords='offset points', 
        xytext=(30,0),  # Move further to the right
        ha='left', 
        va='center',
        fontsize='large',  # Increase font size
        weight='bold'  # Make the text bold
    )

    plt.savefig("tsne_plot.png", dpi=300)
    plt.show()
    
    
        # Create a dictionary to store cluster counts
    cluster_counts = {i: {} for i in range(k)}
    for fatty_acid_length, cluster_label in zip(fatty_acids_lengths, cluster_labels):
        if fatty_acid_length not in cluster_counts[cluster_label]:
            cluster_counts[cluster_label][fatty_acid_length] = 0
        cluster_counts[cluster_label][fatty_acid_length] += 1
    
    cluster_accuracies = {}
    fatty_acid_lengths_most_common = {}
    for cluster, counts in cluster_counts.items():
        total = sum(counts.values())
        fatty_acid_length_most_common = max(counts, key=counts.get)
        max_count = counts[fatty_acid_length_most_common]
        cluster_accuracies[cluster] = max_count / total if total > 0 else 0
        fatty_acid_lengths_most_common[cluster] = fatty_acid_length_most_common


# Create a pie chart
    fig3 = plt.figure(figsize=(8, 6))
    ax_pie = fig3.add_subplot(111)
    ax_pie.pie(
        cluster_accuracies.values(), 
        labels=[f'Cluster {i}, Fatty Acid Length {fatty_acid_lengths_most_common[i]}' for i in cluster_accuracies.keys()],
        autopct='%1.1f%%'
    )
    ax_pie.set_title('Cluster Accuracies')
    plt.savefig("pie_chart.png", dpi=300)
    plt.show()
    
#counts the number of sequences for each fatty acid length in each cluster.
#Then, it calculates the accuracy of each cluster as the number of sequences with the most common fatty acid length divided by the total number of sequences in the cluster. 

    #Cada linha do heatmap representa uma sequência e cada coluna representa um vetor de embedding de palavra
#A cor em cada célula do heatmap indica o valor do embedding de palavra para uma determinada sequência.
#Padrões semelhantes podem indicar sequências que são semanticamente ou funcionalmente relacionadas.

# Generate the linkage matrix for the heatmap
    linkage_matrix = linkage(concatenated_embeddings, method='ward')
    fig2 = plt.figure(figsize=(18, 9))

# Create the heatmap with color legend
    ax_heatmap = fig2.add_subplot(111)
    sns.heatmap(
        concatenated_embeddings,
        cmap=sns.color_palette(list(color_map.values())),
        cbar=True,
        cbar_ax=fig2.add_axes([0.92, 0.1, 0.02, 0.8]),
        ax=ax_heatmap,
    )
    ax_heatmap.set_xlabel('Word Embeddings')
    ax_heatmap.set_ylabel('Sequences')
    plt.savefig("heatmap.png", dpi=300)

    plt.show()

    # Save sequences and print information
    save_sequences_and_print_info(alignment, cluster_labels, table_data)
    
def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_file> <table_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    table_file = sys.argv[2]
    run_mafft(input_file, table_file)

if __name__ == "__main__":
    main()

