import matplotlib.pyplot as plt
import numpy as np


def showData(file_name):
    for name in file_name:
        fname = ".\\SpectralClustering\\data\\" + name
        X = []
        y = []
        labels = {}
        data = []
        with open(fname, 'r') as fp:
            lines = fp.readlines()
            for line in lines:
                point = [float(x) for x in line.split()]
                X.append(point[0])
                y.append(point[1])
                if point[2] in labels.keys():
                    labels[point[2]] += 1
                else:
                    labels.update({point[2]: 1})
                data.append(point)

        sorted_dict = sorted(labels.items(), key=lambda item: item[1], reverse=True)
        # print(data)
        value1, value2 = sorted_dict[:2]
        # print(f"value1  {value1[0]}")
        # print(f"value2  {value2[0]}")
        mean = value1[0] + value2[0] / 2
        # print(mean)
        label_2 = []

        for d in data:
            if abs(d[2] - value1[0]) > abs(d[2] - value2[0]):
                label_2.append(1)
            else:
                label_2.append(2)
        # print(label_2)
        plt.scatter(X, y, c=label_2)
        plt.show()


if __name__ == "__main__":
    fname = ["circle__result.txt","moon__result.txt"]
    showData(fname)
