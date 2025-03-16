import numpy as np
import itertools
from scipy.linalg import eigh
from scipy.linalg import expm
from scipy.sparse.linalg import eigs


def generate_interaction_coefficients(N):
    coefficients = np.random.normal(0, 1, size=(N, N, N, N))
    return coefficients


def construct_hamiltonian(N, J):
    dim = 2 ** (N // 2)
    hamiltonian = np.zeros((dim, dim))

    for i, j, k, l in itertools.product(range(N), repeat=4):
        if i < j and k < l and len(set([i, j, k, l])) == 4:
            term = np.zeros((dim, dim))
            for p, q in itertools.product(range(dim), repeat=2):
                p_binary = format(p, f"0{N//2}b")
                q_binary = format(q, f"0{N//2}b")
                new_q_binary = list(q_binary)
                new_q_binary[i // 2] = q_binary[j // 2]
                new_q_binary[j // 2] = q_binary[i // 2]
                new_q_binary[k // 2] = q_binary[l // 2]
                new_q_binary[l // 2] = q_binary[k // 2]
                sign = (-1) ** (
                    sum([int(p_binary[x]) for x in [i // 2, j // 2, k // 2, l // 2]])
                )
                term[p, q] = sign * int("".join(new_q_binary), 2) == p
            hamiltonian += J[i, j, k, l] * term

    return hamiltonian


def two_point_correlation(N, eigenvalues, eigenvectors, tau):
    G = np.zeros((N, N))
    Z = np.sum(np.exp(-eigenvalues * tau))
    for i, j in itertools.product(range(N), repeat=2):
        for m, n in itertools.product(range(len(eigenvalues)), repeat=2):
            G[i, j] += (
                np.conj(eigenvectors[m, i])
                * eigenvectors[m, j]
                * eigenvectors[n, i]
                * np.conj(eigenvectors[n, j])
                * np.exp(-tau * (eigenvalues[m] - eigenvalues[n]))
                / Z
            )
    return G


def calculate_observables(N, eigenvalues, eigenvectors, tau_values):
    two_point_correlations = np.zeros((len(tau_values), N, N))

    for t_idx, tau in enumerate(tau_values):
        two_point_correlations[t_idx] = two_point_correlation(
            N, eigenvalues, eigenvectors, tau
        )

    return two_point_correlations


def boltzmann_distribution(eigenvalues, T):
    Z = np.sum(np.exp(-eigenvalues / T))
    return np.exp(-eigenvalues / T) / Z


def entropy(eigenvalues, T):
    p = boltzmann_distribution(eigenvalues, T)
    return -np.sum(p * np.log(p))


def calculate_entropy(eigenvalues, T_values):
    entropies = np.zeros_like(T_values)
    for i, T in enumerate(T_values):
        entropies[i] = entropy(eigenvalues, T)
    return entropies


def specific_heat(eigenvalues, T):
    p = boltzmann_distribution(eigenvalues, T)
    energy = np.sum(p * eigenvalues)
    energy_squared = np.sum(p * (eigenvalues**2))
    return (energy_squared - energy**2) / (T**2)


def calculate_specific_heat(eigenvalues, T_values):
    specific_heats = np.zeros_like(T_values)
    for i, T in enumerate(T_values):
        specific_heats[i] = specific_heat(eigenvalues, T)
    return specific_heats


def calculate_spectral_function(eigenvalues, eigenvectors, beta):
    N = len(eigenvalues)
    rho = np.zeros(N)
    for i in range(N):
        for j in range(N):
            energy_diff = eigenvalues[i] - eigenvalues[j]
            rho[i] += np.abs(eigenvectors[:, i] @ eigenvectors[:, j]) ** 2 * np.exp(
                -beta * energy_diff
            )
    return rho


def run_syk_simulation(N):
    J = generate_interaction_coefficients(N)
    hamiltonian = construct_hamiltonian(N, J)
    eigenvalues, eigenvectors = eigh(hamiltonian)
    return eigenvalues, eigenvectors


def run_multiple_syk_simulations(N, num_simulations):
    eigenvalues_list = []
    eigenvectors_list = []

    for i in range(num_simulations):
        np.random.seed(
            i
        )  # Set the random seed to a different value for each simulation
        eigenvalues, eigenvectors = run_syk_simulation(N)
        eigenvalues_list.append(eigenvalues)
        eigenvectors_list.append(eigenvectors)

    return eigenvalues_list, eigenvectors_list


# ----------------- #


def simplified_syk_hamiltonian(N, J):
    """
    Create a simplified SYK Hamiltonian.
    """
    H = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if i != j:
                H[i, j] = J * np.random.normal()
    return H


def calculate_otoc(N, J, num_samples, time_steps):
    # Initialize arrays for time and OTOC values
    times = np.linspace(0, 10, time_steps)
    oto_correlations = np.zeros(time_steps)

    # Create initial operators A and B
    A = np.zeros((N, N))
    B = np.zeros((N, N))

    # Fill diagonals with random values
    np.fill_diagonal(A, np.random.rand(N))
    np.fill_diagonal(B, np.random.rand(N))

    for sample in range(num_samples):
        # Generate simplified SYK Hamiltonian
        H = simplified_syk_hamiltonian(N, J)

        # Compute the time evolution of the operators A and B
        for t_idx, t in enumerate(times):
            U_t = expm(-1j * H * t)
            U_t_dag = expm(1j * H * t)

            # Compute the commutator [A(t), B]
            A_t = U_t @ A @ U_t_dag
            commutator = A_t @ B - B @ A_t

            # Calculate the OTOC
            oto_correlations[t_idx] += np.abs(np.trace(commutator @ commutator))

    # Normalize the OTOC by the number of samples
    oto_correlations /= num_samples

    return times, oto_correlations
