import pandas as pd
import sys
import itertools

def assign_color(taxon, color_dict, used_colors, available_colors):
    if taxon in color_dict:
        return color_dict[taxon]
    else:
        color = next(available_colors)
        while color in used_colors:
            color = next(available_colors)
        color_dict[taxon] = color
        used_colors.add(color)
        return color

def get_order_value(row, taxonomic_hierarchy):
    for rank in taxonomic_hierarchy:
        if pd.notna(row[rank]):
            return row[rank]
    return "Unknown order"

def main(file_path):
    color_dict = {
        "Candidatus Termititenacales": "#ff6347",
        "Chroococcales": "#40e0d0",
        "Coleofasciculales": "#ffd700",
        "Gloeobacterales": "#6a5acd",
        "Gloeomargaritales": "#98fb98",
        "Leptolyngbyales": "#dda0dd",
        "Nostocales": "#cd5c5c",
        "Oculatellales": "#b0e0e6",
        "Oscillatoriales": "#ffa07a",
        "Pleurocapsales": "#00ced1",
        "Pseudanabaenales": "#bdb76b",
        "Spirulinales": "#f08080",
        "Synechococcales": "#ffe4b5",
        "Thermostichales": "#fafad2",
        "Vampirovibrionales": "#000000",
        "Não conhecido": "#a9a9a9"
    }

    # Define a paleta de cores, excluindo verdes muito parecidos
    available_colors = itertools.cycle([
        "#7fffd4", "#8a2be2", "#a52a2a", "#deb887", "#5f9ea0", "#d2691e",
        "#ff7f50", "#6495ed", "#fff8dc", "#dc143c", "#00ffff", "#00008b",
        "#008b8b", "#b8860b", "#8b008b", "#556b2f", "#ff8c00", "#9932cc",
        "#8b0000", "#e9967a", "#8fbc8f", "#483d8b", "#2f4f4f", "#9400d3",
        "#ff1493", "#00bfff", "#696969", "#1e90ff", "#b22222", "#fffaf0",
        "#228b22", "#dcdcdc", "#f8f8ff", "#daa520", "#808080", "#008000",
        "#f0fff0", "#ff69b4", "#4b0082", "#f0e68c", "#e6e6fa", "#fff0f5",
        "#7cfc00", "#fffacd", "#add8e6", "#f08080", "#e0ffff", "#d3d3d3",
        "#90ee90", "#ffb6c1", "#ffa07a", "#20b2aa", "#87cefa", "#778899",
        "#b0c4de", "#ffffe0", "#32cd32", "#faf0e6", "#ff00ff", "#800000",
        "#66cdaa", "#0000cd", "#ba55d3", "#9370db", "#3cb371", "#7b68ee",
        "#00fa9a", "#48d1cc", "#c71585", "#191970", "#f5fffa", "#ffe4e1",
        "#ffe4b5", "#ffdead", "#fdf5e6", "#6b8e23", "#ff4500", "#da70d6",
        "#eee8aa", "#afeeee", "#db7093", "#ffefd5", "#ffdab9", "#cd853f",
        "#ffc0cb", "#dda0dd", "#800080", "#663399", "#ff0000", "#bc8f8f",
        "#4169e1", "#8b4513", "#fa8072", "#f4a460", "#2e8b57", "#fff5ee",
        "#a0522d", "#87ceeb", "#6a5acd", "#708090", "#fffafa", "#00ff7f",
        "#4682b4"
    ])

    taxonomic_hierarchy = ["order", "superorder", "subclass", "class", "superclass", "subphylum", "phylum", "superphylum", "kingdom", "superkingdom"]

    data = pd.read_csv(file_path, sep='\t', dtype=str)  # Ajuste o separador se necessário

    # Trata espaços e caracteres invisíveis no cabeçalho das colunas
    data.columns = data.columns.str.strip()

    # Substitui valores ausentes na coluna 'order' pelo valor do ranking precedente ou "Unknown order"
    data['order'] = data.apply(lambda row: get_order_value(row, taxonomic_hierarchy), axis=1)

    used_colors = set()
    data['color'] = data['order'].apply(lambda order: assign_color(order, color_dict, used_colors, available_colors))

    result = data[['Genome ID', 'order', 'color']]
    result.to_csv("output_with_colors.csv", index=False)
    print("Nova tabela criada com sucesso em 'output_with_colors.csv'.")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Por favor, forneça o caminho do arquivo CSV como argumento.")
    else:
        main(sys.argv[1])
