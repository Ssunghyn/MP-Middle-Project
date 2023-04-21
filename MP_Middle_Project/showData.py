import matplotlib.pyplot as plt
import sys
import os

path = os.getcwd()


def showData(file_name):
    fig, axes = plt.subplots(1, 2)
    i = 0
    for name in file_name[1:]:
        fname = path + "\\data\\" + name
        X = []
        y = []
        labels = []
        with open(fname, 'r') as fp:
            lines = fp.readlines()
            for line in lines:
                point = [float(x) for x in line.split()]
                X.append(point[0])
                y.append(point[1])
                labels.append(point[2])

        axes[i].scatter(X, y, c=labels)
        axes[i].set_title(name[:-4])
        i += 1

    plt.show()


if __name__ == "__main__":
    fname = sys.argv
    showData(fname)
