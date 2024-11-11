import numpy as np
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram

def plot_emst(x, y, labels, edges, save_path="emst.png"):
    plt.figure(figsize=(10, 8))
    plt.plot(x, y, "o")
    
    for xx, yy, ll in zip(x, y, labels):
        plt.text(xx, yy, ll, fontsize=8)

    for ed in edges:
        plt.plot(ed[0], ed[1], 'k-', linewidth=0.5)

    plt.title("Euclidean Minimum Spanning Tree (EMST)")
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    plt.grid(True)
    plt.savefig(save_path, dpi=300)
    plt.close()

def plot_dendrogram(matrix, labels, save_path="dendrogram.png"):
    plt.figure(figsize=(10, 8))
    dendrogram(matrix, color_threshold=0, labels=labels, leaf_rotation=90, leaf_font_size=8)
    plt.title("Dendrogram")
    plt.xlabel("Cluster Labels")
    plt.ylabel("Distance")
    plt.savefig(save_path, dpi=300)
    plt.close()

if __name__ == "__main__":
    x, y, labels, edges = [], [], [], []
    lett = ord('A')
    
    with open("emst.txt", "r") as file:
        for line in file:
            datum = line.split()
            if int(datum[0]) == 0:
                x.append(float(datum[1]))
                y.append(float(datum[2]))
                labels.append(chr(lett))
                lett += 1
            else:
                edges.append([[x[int(datum[1])], x[int(datum[2])]], 
                              [y[int(datum[1])], y[int(datum[2])]]])

    plot_emst(x, y, labels, edges, "emst.png")

    matrix = []
    with open("dendro.txt", "r") as file:
        for line in file:
            matrix.append([float(x) for x in line.split()])
    matrix = np.array(matrix)

    plot_dendrogram(matrix, labels, "dendrogram.png")
