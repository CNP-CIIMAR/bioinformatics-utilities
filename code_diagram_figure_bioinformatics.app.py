### Autor: Leandro de Mattos, CNP Team - Leao Lab.
import matplotlib
matplotlib.use('Agg')  # Define o backend de gravação
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import matplotlib.patches as patches

fig, ax = plt.subplots(figsize=(20, 15))  # Aumentar o tamanho da figura

# Define colors
circle_color = "#1f77b4"
rectangle_color = "#aec7e8"
line_color = "#228B22"  # Verde intenso para destacar as linhas

# Draw center circle
center_circle = patches.Circle((0.5, 0.5), 0.15, fill=True, color=circle_color, edgecolor=line_color)  # Aumentar o tamanho do círculo
ax.add_patch(center_circle)
ax.text(0.5, 0.5, 'Bioinformatics', va='center', ha='center', color="white", fontsize=24, weight='bold')

# Draw rectangles and lines on the left
left_texts = ['Biological databases dev.', 'Molecular Structures\nPredictions',
              'Software dev. (BLAST,\nantiSMASH)', 'Comparative Genomics\nand Phylogenomics',
              'Molecular Evolution and\nPhylogenetics', '“Omics” fields:\nAnalysis and Interpretation']
for i, text in enumerate(left_texts):
    y = 0.9 - i * 0.15  # Calculate y position
    rectangle = patches.Rectangle((0.15, y-0.05), 0.2, 0.1, fill=True, color=rectangle_color)  # Draw rectangle
    ax.add_patch(rectangle)
    ax.text(0.25, y, text, va='center', ha='center', fontsize=16, wrap=True, color="black")
    ax.plot([0.35, 0.4], [y, 0.5], color=line_color, lw=2)  # Linhas mais destacadas com linewidth (lw) igual a 2

# Draw rectangles and lines on the right
right_texts = ['Synthetic Biology: Enzyme, \nMetabolic pathways\n and BGCs mining',
               'Synthetic Biology: Enzyme, \nProtein\n and Metabolic Engineering',
               'Pharmacogenomic', 'Epidemiology and\n Vaccine development',
               'Personalized Medicine and\n Diseases diagnosis', 'Ecological\n and\n Biodiversity studies',
               'Forensic Sciences']
for i, text in enumerate(right_texts):
    y = 0.9 - i * 0.13  # Calculate y position
    rectangle = patches.Rectangle((0.65, y-0.05), 0.2, 0.1, fill=True, color=rectangle_color)  # Draw rectangle
    ax.add_patch(rectangle)
    ax.text(0.75, y, text, va='center', ha='center', fontsize=16, wrap=True, color="black")
    ax.plot([0.65, 0.6], [y, 0.5], color=line_color, lw=2)  # Linhas mais destacadas com linewidth (lw) igual a 2

ax.set_xlim([0, 1])
ax.set_ylim([0, 1])
ax.axis('off')  # Hide axis
plt.savefig("figure_esquema.png")  # save figure to file
