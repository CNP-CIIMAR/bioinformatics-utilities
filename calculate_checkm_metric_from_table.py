import numpy as np
from collections import defaultdict
from pathlib import Path

def read_table_from_csv(file_path):
    data = []
    with open(file_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            data.append(row)
    return data

def safe_float(value):
    try:
        return float(value)
    except (ValueError, TypeError):
        return None

def calculate_metrics(data):
    id_to_values = defaultdict(list)

    for row in data:
        accession = row["Assembly Accession"]
        completeness = safe_float(row["CheckM completeness"])
        completeness_percentile = safe_float(row["CheckM completeness percentile"])

        if completeness is not None and completeness_percentile is not None:
            id_to_values[accession].append((completeness, completeness_percentile))

    max_completeness_values = [max(vals, key=lambda x: x[0])[0] for vals in id_to_values.values()]
    max_completeness_percentile_values = [max(vals, key=lambda x: x[1])[1] for vals in id_to_values.values()]

    completeness_mean = np.mean(max_completeness_values)
    completeness_percentile_mean = np.mean(max_completeness_percentile_values)

    return completeness_mean, completeness_percentile_mean

def main():
    if len(sys.argv) != 2:
        print("Usage: python calculate_checkm_metrics.py <file_path>")
        return

    file_path = sys.argv[1]
    data = read_table_from_csv(file_path)

    completeness_mean, completeness_percentile_mean = calculate_metrics(data)

    # Salvar os resultados em um arquivo de saída
    output_file = Path(file_path).stem + "_resultados.txt"
    with open(output_file, "w") as f:
        f.write(f"Mean CheckM completeness: {completeness_mean:.2f}\n")
        f.write(f"Mean CheckM completeness percentile: {completeness_percentile_mean:.2f}\n")

    print(f"Resultados salvos em {output_file}")

if __name__ == "__main__":
    main()
