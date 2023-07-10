# Authors: Leandro de Mattos Pereira
# Pedro Leao - Team Leader - CNP
import os

# Step 1
with open('Genomes_id') as file:
    for line in file:
        genome_id = line.strip()
        os.mkdir(genome_id)

# Step 2
with open('Genomes_id') as file:
    for line in file:
        genome_id = line.strip()
        os.chdir(genome_id)
        os.system(f'wget -qN "ftp://ftp.patricbrc.org/genomes/{genome_id}/{genome_id}*"')
        os.chdir('..')
