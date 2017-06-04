from __future__ import print_function

import matplotlib
import matplotlib.pyplot as plt
import numpy as np


def maxwell(vel, T):
    m = 3.9529263e-25  # Mass of U-238 in kg
    k_b = 1.38065e-23  # Boltzmann constant in J / K

    return (np.sqrt( (m / 2 / np.pi / k_b / T)**3 ) * 4 * np.pi * vel**2
            * np.exp(-m * vel**2 / 2 / k_b / T))


def cxs(vel, T):
    return vel * maxwell(vel, T)


if __name__ == '__main__':
    with open('vel.txt') as fh:
        lines = fh.readlines()
    dat = [float(l.strip()) for l in lines]
    plt.hist(dat, bins=100, normed=True)

    v_grid = np.linspace(0, 1000, 200)
    ref = cxs(v_grid, 300)
    ref /= np.trapz(ref, v_grid)
    plt.plot(v_grid, ref)

    plt.show()
