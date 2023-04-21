import matplotlib.pyplot as plt
from sklearn.datasets import make_circles, make_moons, make_blobs
from numpy import random
import os
import sys

random_state = 8
random.seed(11)
n_sample = 2000
path = os.getcwd()


def get_circle_data():
    X, y = make_circles(n_samples=n_sample, factor=0.3, noise=0.03, random_state=random_state)
    plt.scatter(X[:, 0], X[:, 1])
    plt.show()
    x_to_str = [f"{p[0]} {p[1]}" for p in X]
    point_str = "\n".join(x_to_str)
    # print(f"{X[0][0]} {X[0][1]}")
    # print(point_str.split("\n")[0])
    with open(path + "\\data\\circle_data.txt", "w") as f:
        f.write(point_str)


def get_moon_data():
    X, y = make_moons(n_samples=n_sample, noise=0.03, random_state=random_state)
    plt.scatter(X[:, 0], X[:, 1])
    plt.show()
    x_to_str = [f"{p[0]} {p[1]}" for p in X]
    point_str = "\n".join(x_to_str)
    with open(path + "\\data\\moon_data.txt", "w") as f:
        f.write(point_str)


def get_blob_data():
    X, y = make_blobs(n_samples=n_sample, random_state=random_state)
    plt.scatter(X[:, 0], X[:, 1])
    plt.show()
    x_to_str = [f"{p[0]} {p[1]}" for p in X]
    point_str = "\n".join(x_to_str)
    with open(path + "\\data\\blob_data.txt", "w") as f:
        f.write(point_str)


if __name__ == "__main__":
    n_sample = int(sys.argv[1])
    directory = "data"
    if not os.path.exists(directory):
        os.makedirs(directory)
    get_circle_data()
    get_moon_data()
    get_blob_data()
