import numpy as np
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram

if __name__ == "__main__":
    x = []
    y = []
    l = []
    edges = []
    lett = ord('A')
    with open("emst.txt", "r+") as file:
        for line in file:
            datum = line.split()
            if (int(datum[0]) == 0):
                x.append(np.double(datum[1]))
                y.append(np.double(datum[2]))
                l.append(chr(lett))
                lett += 1
            else:
                edges.append([[x[int(datum[1])], x[int(datum[2])]], [y[int(datum[1])], y[int(datum[2])]]])

    fig = plt.plot(x, y, "o")
    for xx, yy, ll in zip(x, y, l):
        plt.text(xx, yy, ll)
    for ed in edges:
        plt.plot(ed[0], ed[1], 'k-')
    plt.savefig("emst.png")

    arr = []
    with open("dendro.txt", "r+") as file:
        for line in file:
            arr.append(np.array([np.double(x) for x in line.split()]));

    arr = np.array(arr);
    
    fig = plt.figure(figsize=(3,2.3))
    dn = dendrogram(arr, color_threshold=0, labels=l)
    fig.savefig("dendrogram.png")
