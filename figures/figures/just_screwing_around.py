import numpy as np
from matplotlib import pyplot as plt

if __name__ == '__main__':
    x = [*range(1000)]
    y = [*range(1000)]

    plt.plot(x, y)
    plt.yscale('log')
    plt.show()
