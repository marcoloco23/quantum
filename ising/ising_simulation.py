import numpy as np
import random


def initial_lattice(N):
    return np.random.choice([-1, 1], size=(N, N, N))


def local_energy(spin_config, i, j, k, J, h, N):
    return (
        -J
        * spin_config[i, j, k]
        * (
            spin_config[i, (j + 1) % N, k]
            + spin_config[(i + 1) % N, j, k]
            + spin_config[i, (j - 1) % N, k]
            + spin_config[(i - 1) % N, j, k]
            + spin_config[i, j, (k + 1) % N]
            + spin_config[i, j, (k - 1) % N]
        )
        - h * spin_config[i, j, k]
    )


def total_energy(spin_config, J, h):
    N = len(spin_config)
    E = 0
    for i in range(N):
        for j in range(N):
            for k in range(N):
                E += local_energy(spin_config, i, j, k, J, h, N)
    return E / 2


def metropolis(spin_config, J, h, T, N):
    i, j, k = (
        random.randint(0, N - 1),
        random.randint(0, N - 1),
        random.randint(0, N - 1),
    )
    dE = -2 * local_energy(spin_config, i, j, k, J, h, N)
    if dE <= 0 or random.random() < np.exp(-dE / T):
        spin_config[i, j, k] *= -1
        return dE
    return 0


def simulate(args):
    N, J, h, T, steps, equilibration_steps = args

    spin_config = initial_lattice(N)

    # Equilibration
    for _ in range(equilibration_steps):
        metropolis(spin_config, J, h, T, N)

    # Measurement
    M = 0
    E = total_energy(spin_config, J, h)
    M2 = 0
    E2 = 0
    for _ in range(steps):
        dE = metropolis(spin_config, J, h, T, N)
        E += dE
        M += np.sum(spin_config)
        E2 += E * E
        M2 += np.sum(spin_config) * np.sum(spin_config)

    M_avg = M / steps
    E_avg = E / steps
    M2_avg = M2 / steps
    E2_avg = E2 / steps

    magnetization = M_avg / (N * N * N)
    energy = E_avg / (N * N * N)
    susceptibility = (M2_avg / (N * N * N) - M_avg * M_avg) / (N * N * N * T)
    heat_capacity = (E2_avg - E_avg * E_avg) / (N * N * N * T * T)

    return magnetization, energy, susceptibility, heat_capacity
