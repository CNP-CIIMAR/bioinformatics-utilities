### Autor Leandro de Mattos Pereira
### CNP - Laboratory 16 August - 2024
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import argparse

def convert_to_number(fa_string):
    """Converte notaÃ§Ã£o de Ã¡cidos graxos para valores numÃ©ricos."""
    try:
        if 'C' in fa_string:
            return int(fa_string.split('C')[1].split(':')[0].strip().split(' ')[0].split('(')[0])
    except (ValueError, IndexError):
        pass
    return None

def convert_range_to_numbers(range_string):
    """Converte uma string de intervalo de Ã¡cidos graxos em uma lista de nÃºmeros."""
    if '-' in range_string:
        start, end = range_string.split('-')
        start = convert_to_number(f'C{start.strip()}')
        end = convert_to_number(f'C{end.strip().split(" ")[0].split("(")[0]}')
        if start is not None and end is not None:
            return list(range(start, end + 1))
    else:
        return [convert_to_number(range_string)]
    return []

def process_data(data):
    """Processa os dados da tabela."""
    processed_data = []
    for _, row in data.iterrows():
        if pd.isna(row['Protein']):
            continue
        protein = row['Protein']
        species = row['Specie']
        faal_pred = row['FAALPred']
        adenyl_pred = row.get('AdenylPred', None)
        exp_val = row['Associated variable']
        
        faal_pred_numbers = [convert_to_number(c) for c in faal_pred.split('-') if convert_to_number(c) is not None]
        adenyl_pred_numbers = [num for range_str in (adenyl_pred or '').split('-') for num in convert_range_to_numbers(range_str)] if pd.notna(adenyl_pred) else None
        exp_val_numbers = [convert_to_number(c.split(':')[0]) for c in exp_val.split(', ')]
        
        processed_data.append({
            'protein': protein,
            'species': species,
            'FAALPred': faal_pred_numbers,
            'AdenylPred': adenyl_pred_numbers,
            'ExperimentalValidation': exp_val_numbers
        })
    return processed_data

def plot_protein_specificity(data, output_dir):
    """Gera um Ãºnico grÃ¡fico para todas as proteÃ­nas."""
    fig, ax = plt.subplots(figsize=(18, len(data) * 0.7))  # Tamanho da figura ajustado para publicaÃ§Ã£o
    
    y_position = 0
    yticks = []
    yticklabels = []
    markersize = 7  # Tamanho das bolinhas ajustado

    for entry in data:
        protein = entry['protein']
        faal_pred = entry['FAALPred']
        adenyl_pred = entry.get('AdenylPred', None)
        exp_val = entry['ExperimentalValidation']
        
        y_position += 1
        yticks.append(y_position)
        yticklabels.append(protein)
        
        if faal_pred:
            ax.plot(faal_pred, [y_position - 0.3] * len(faal_pred), 'bo-', markersize=markersize)
        
        if adenyl_pred:
            ax.plot(adenyl_pred, [y_position] * len(adenyl_pred), 'go--', markersize=markersize)
        
        if exp_val:
            ax.plot(exp_val, [y_position + 0.3] * len(exp_val), 'ro-', markersize=markersize)
    
    # Criar uma legenda simplificada
    handles = [
        plt.Line2D([0], [0], color='b', lw=2, marker='o', markersize=markersize, label='FAALPred'),
        plt.Line2D([0], [0], color='g', lw=2, linestyle='--', marker='o', markersize=markersize, label='AdenylPred'),
        plt.Line2D([0], [0], color='r', lw=2, marker='o', markersize=markersize, label='ExperimentalValidation')
    ]
    ax.legend(handles=handles, loc='center left', bbox_to_anchor=(1, 0.5), prop={'weight': 'bold'})
    
    ax.set_xlabel('Specificity: Length of Fatty Acid Chain (C2 to C18)', fontweight='bold')
    ax.set_ylabel('Proteins', fontweight='bold')

    ax.set_xticks(np.arange(2, 19))
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels, fontweight='bold')
    ax.set_xticklabels(np.arange(2, 19), fontweight='bold')
    ax.grid(True)
    
    # Ajustar margens para nÃ£o cortar o grÃ¡fico
    fig.subplots_adjust(top=0.85, bottom=0.15, left=0.2, right=0.8, hspace=0.2, wspace=0.2)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Salvar em formatos de alta resoluÃ§Ã£o
    plt.savefig(f'{output_dir}/all_proteins_specificity.svg', format='svg')
    plt.savefig(f'{output_dir}/all_proteins_specificity.png', format='png', dpi=600)  # Alta resoluÃ§Ã£o para PNG
    plt.savefig(f'{output_dir}/all_proteins_specificity.jpeg', format='jpeg', dpi=600)  # Alta resoluÃ§Ã£o para JPEG
    
    plt.show()

def main():
    parser = argparse.ArgumentParser(description='Generate protein specificity prediction graphs.')
    parser.add_argument('input_file', type=str, help='Path to the input TSV file.')
    parser.add_argument('output_dir', type=str, help='Directory to save the output figures.')
    
    args = parser.parse_args()
    
    input_file = args.input_file
    output_dir = args.output_dir
    
    try:
        data = pd.read_csv(input_file, sep='\t')
        data.dropna(subset=['Protein'], inplace=True)
        processed_data = process_data(data)
        plot_protein_specificity(processed_data, output_dir)
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == '__main__':
    main()

