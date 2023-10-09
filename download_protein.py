import os
import subprocess
import tkinter as tk
from tkinter import filedialog, messagebox
from PIL import ImageTk, Image
from io import BytesIO
import requests
from Bio import Entrez, SeqIO

# Configuração das informações de e-mail para a API do NCBI
Entrez.email = "mattoslmp@gmail.com"

window = tk.Tk()
window.title("Protein Downloader")
window.geometry("500x600")

# Função para selecionar o arquivo de IDs de proteínas
def select_fasta_ids_file():
    fasta_ids_file_path = filedialog.askopenfilename(filetypes=[("Text Files", "*.txt")])
    if fasta_ids_file_path:
        fasta_ids_file_label.config(text=fasta_ids_file_path)

# Botão para selecionar o arquivo de IDs de proteínas
select_fasta_ids_file_button = tk.Button(window, text="Select IDs File", command=select_fasta_ids_file)
select_fasta_ids_file_button.pack()

# Label para exibir o arquivo de IDs de proteínas selecionado
fasta_ids_file_label = tk.Label(window, text="", wraplength=400)
fasta_ids_file_label.pack()

# Função para selecionar o diretório de saída
def select_output_directory():
    output_directory = filedialog.askdirectory()
    if output_directory:
        output_directory_label.config(text=output_directory)

# Botão para selecionar o diretório de saída
select_output_directory_button = tk.Button(window, text="Select Output Directory", command=select_output_directory)
select_output_directory_button.pack()

# Label para exibir o diretório de saída selecionado
output_directory_label = tk.Label(window, text="", wraplength=400)
output_directory_label.pack()

# Função para fazer o download das sequências FASTA usando Biopython
def download_sequences():
    fasta_ids_file = fasta_ids_file_label.cget("text")
    output_directory = output_directory_label.cget("text")

    if not fasta_ids_file or not output_directory:
        messagebox.showwarning("Missing Information", "Please select the IDs file and output directory.")
        return

    try:
        with open(fasta_ids_file, 'r') as file:
            fasta_ids = file.read().splitlines()

        all_sequences = []  # List to store all downloaded sequences

        for fasta_id in fasta_ids:
            # Check if the ID corresponds to UniProtKB
            if fasta_id.startswith("UniProtKB:"):
                uniprot_id = fasta_id.split(":")[1]
                # Faz o download da sequência FASTA usando o UniProt website
                url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
                response = requests.get(url)
                if response.ok:
                    fasta_data = response.text
                    record = SeqIO.read(fasta_data.splitlines(), "fasta")
                    # Cria o arquivo de saída com a sequência FASTA
                    output_file = os.path.join(output_directory, f'{fasta_id}.fasta')
                    SeqIO.write(record, output_file, "fasta")
                    all_sequences.append(record)  # Add sequence to the list

            # Check if the ID corresponds to PDB
            elif fasta_id.startswith("PDB:"):
                pdb_id = fasta_id.split(":")[1]
                # Faz o download da sequência FASTA usando o PDB website
                url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
                response = requests.get(url)
                if response.ok:
                    fasta_data = response.text
                    record = SeqIO.read(fasta_data.splitlines(), "fasta")
                    # Cria o arquivo de saída com a sequência FASTA
                    output_file = os.path.join(output_directory, f'{fasta_id}.fasta')
                    SeqIO.write(record, output_file, "fasta")
                    all_sequences.append(record)  # Add sequence to the list

            else:
                # Assume NCBI protein database
                # Faz o download da sequência FASTA usando o Biopython
                handle = Entrez.efetch(db="protein", id=fasta_id, rettype="fasta", retmode="text")
                record = SeqIO.read(handle, "fasta")
                handle.close()
                # Cria o arquivo de saída com a sequência FASTA
                output_file = os.path.join(output_directory, f'{fasta_id}.fasta')
                SeqIO.write(record, output_file, "fasta")
                all_sequences.append(record)  # Add sequence to the list

        # Cria o arquivo FASTA contendo todas as sequências
        all_sequences_file = os.path.join(output_directory, "all_sequences.fasta")
        SeqIO.write(all_sequences, all_sequences_file, "fasta")

        messagebox.showinfo("Download Completed", "The FASTA sequences have been downloaded successfully!")
    except Exception as e:
        messagebox.showerror("Error", str(e))

# Botão para iniciar o download das sequências FASTA
download_button = tk.Button(window, text="Download FASTA Sequences", command=download_sequences)
download_button.pack()

# Função para importar a logomarca do diretório fornecido
def import_logo():
    logo_path = filedialog.askopenfilename(filetypes=[("Image Files", "*.png;*.jpg;*.jpeg")])
    if logo_path:
        try:
            image = Image.open(logo_path)
            logo_image = ImageTk.PhotoImage(image)

            # Cria um label para exibir a logomarca
            logo_label = tk.Label(window, image=logo_image)
            logo_label.pack()

            # Atualiza a imagem da logomarca
            logo_label.image = logo_image
        except Exception as e:
            messagebox.showerror("Error", str(e))

# Botão para importar a logomarca
import_logo_button = tk.Button(window, text="Import Logo", command=import_logo)
import_logo_button.pack()

# Label com mais informações e link
more_info_label = tk.Label(window, text="More Information: ", fg="blue", cursor="hand2")
more_info_label.pack()
more_info_label.bind("<Button-1>", lambda e: os.system("start https://leaolab.wixsite.com/leaolab"))

# Label com os direitos autorais
copyright_label = tk.Label(window, text="© 2023 CNP Team. All rights reserved.")

# Posiciona os elementos da interface gráfica
select_fasta_ids_file_button.pack()
fasta_ids_file_label.pack()
select_output_directory_button.pack()
output_directory_label.pack()
download_button.pack()
import_logo_button.pack()
more_info_label.pack()
copyright_label.pack()

# Inicia o loop principal da janela
window.mainloop()

