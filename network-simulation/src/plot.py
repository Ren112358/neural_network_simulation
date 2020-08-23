import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(style='darkgrid')


def main():
    time = np.genfromtxt('../data/time.txt')
    I_ext = np.genfromtxt('../data/input_currents.txt', delimiter=',').T
    V = np.genfromtxt('../data/potentials.txt', delimiter=',').T

    N = I_ext.shape[0]
    ymin = np.min(I_ext)
    ymax = np.max(I_ext)
    mergin = 0.1 * (ymax - ymin)
    for i in range(N):
        plt.subplot(2, N, i + 1)
        plt.ylim(ymin - mergin, ymax + mergin)
        plt.plot(time, I_ext[i])
    for i in range(N):
        plt.subplot(2, N, i + N + 1)
        plt.plot(time, V[i])
    plt.show()


if __name__ == '__main__':
    main()
