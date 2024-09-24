import numpy as np

def read_data(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    data = [line.strip().split('\t') for line in lines[1:]]
    return data

def bin_escores(data, bin_percentage):
    escores = np.array([float(row[2]) for row in data])
    bin_count = int(100 / bin_percentage)
    bin_edges = np.linspace(escores.max(), escores.min(), bin_count + 1)
    binned_escores = np.digitize(escores, bin_edges, right=True)
    return escores, binned_escores, bin_edges

def calculate_stats(escores, binned_escores, bin_edges):
    stats = []
    for i in range(len(bin_edges)-1):
        indices = np.where(binned_escores == i + 1)[0]
        bin_escores = escores[indices]
        avg = np.mean(bin_escores) if bin_escores.size > 0 else 0
        std = np.std(bin_escores) if bin_escores.size > 0 else 0
        stats.append((bin_edges[i], bin_edges[i+1], avg, std))
    return stats

def main():
    file_path = 'GATA4/GATA4_anti-GST/GATA4_anti-GST_8mers_top_enrichment.txt'
    for i in [10,20,50]:
        bin_percentage = i;
        data = read_data(file_path)
        escores, binned_escores, bin_edges = bin_escores(data, bin_percentage)
        stats = calculate_stats(escores, binned_escores, bin_edges)
        
        print(f"Binning {i}%")
        metrics = []
        for i, (bin_start, bin_end, avg, std) in enumerate(stats):
            metric = avg/std
            metrics.append(metric)
        for i, (bin_start, bin_end, avg, std) in enumerate(stats):
            weight = metrics[i]/sum(metrics)
            print(f"Bin {i+1}: {bin_start:.4f} to {bin_end:.4f} - Average: {avg:.4f}, Std Dev: {std:.4f}, Metric: {metrics[i]:.4f}, Weight: {weight:.4f}")
        
        print("***************")
        print("\n")

if __name__ == "__main__":
    main()
