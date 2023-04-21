import matplotlib.pyplot as plt
from sklearn import datasets
import numpy as np

np.random.seed(0)
random_state = 8
n_sample = 2000

def get_circle_data():
    X, y = datasets.make_circles(n_samples=n_sample, factor=0.3, noise=0.03)
    plt.scatter(X[:, 0], X[:, 1])
    plt.show()
    x_to_str = [f"{p[0]} {p[1]}" for p in X]
    point_str = "\n".join(x_to_str)
    # print(f"{X[0][0]} {X[0][1]}")
    # print(point_str.split("\n")[0])
    with open("SpectralClustering/data/circle_data.txt", "w") as f:
        f.write(point_str)


def get_moon_data():
    X, y = datasets.make_moons(n_samples=n_sample, noise=0.03, random_state=random_state)
    plt.scatter(X[:, 0], X[:, 1])
    plt.show()
    x_to_str = [f"{p[0]} {p[1]}" for p in X]
    point_str = "\n".join(x_to_str)
    with open("SpectralClustering/data/moon_data.txt", "w") as f:
        f.write(point_str)


def get_blob_data():
    X, y = datasets.make_blobs(n_samples=n_sample, random_state=random_state)
    plt.scatter(X[:, 0], X[:, 1])
    plt.show()
    x_to_str = [f"{p[0]} {p[1]}" for p in X]
    point_str = "\n".join(x_to_str)
    with open("SpectralClustering/data/blob_data.txt", "w") as f:
        f.write(point_str)


if __name__ == "__main__":
    get_circle_data()
    get_moon_data()
    get_blob_data()
