import numpy as np
from numpy import kron
from scipy.linalg import expm
import matplotlib.pyplot as plt

###############################################################################
# 1. Define basic Pauli matrices and helpers
###############################################################################
I2 = np.array([[1, 0], [0, 1]], dtype=complex)

sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)

sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)


def partial_trace(rho, subsystem=0):
    """
    Compute the partial trace of a 2-qubit density matrix over one qubit.

    Args:
        rho (np.ndarray): 4x4 density matrix (two qubits).
        subsystem (int): Which qubit to trace out.
                        0 -> trace out the first qubit,
                        1 -> trace out the second qubit.

    Returns:
        np.ndarray: 2x2 density matrix of the remaining qubit.
    """
    # Reshape rho into (2,2,2,2)
    rho_reshaped = np.reshape(rho, (2, 2, 2, 2))
    if subsystem == 0:
        # Trace out the first qubit
        return np.einsum("ijik->jk", rho_reshaped)
    else:
        # Trace out the second qubit
        return np.einsum("ijkj->ik", rho_reshaped)


def density_matrix(state):
    """Given a state vector, return the density matrix |psi><psi|."""
    return np.outer(state, state.conjugate())


###############################################################################
# 2. Build the total Hamiltonian for a 2-qubit system
###############################################################################
# Parameters:
omega_x = 1.0  # "energy scale" for the system qubit
omega_y = 0.0  # "energy scale" for the environment qubit (set to 0 for demo)
g = 1.0  # interaction strength

# H_x acts on qubit x only
H_x = omega_x * sigma_z
# H_y acts on qubit y only
H_y = omega_y * sigma_z
# H_int couples x and y; e.g., H_int = g * sigma_x (x) tensor sigma_x (y)
H_int = g * kron(sigma_x, sigma_x)

# Convert H_x, H_y to 4x4 operators by tensoring with identity
H_x_total = kron(H_x, I2)
H_y_total = kron(I2, H_y)

# Total Hamiltonian
H_total = H_x_total + H_y_total + H_int

###############################################################################
# 3. Define time evolution parameters
###############################################################################
t_final = 5.0
dt = 0.05
timesteps = int(t_final / dt)
times = np.linspace(0, t_final, timesteps)


# Time evolution operator for one small step dt:
# U(dt) = exp(-i * H_total * dt / hbar). We set hbar=1 for simplicity.
def time_evolution_operator(H, dt):
    return expm(-1j * H * dt)


U_dt = time_evolution_operator(H_total, dt)

###############################################################################
# 4. Choose an initial state and set up the simulation
###############################################################################
# Let's say system qubit x starts in |0>, environment qubit y in |0>,
# but we give them a superposition or phase.

# Basic computational basis states for a single qubit:
zero = np.array([1, 0], dtype=complex)
one = np.array([0, 1], dtype=complex)

# Build 2-qubit basis states:
zero_zero = kron(zero, zero)  # |00>
zero_one = kron(zero, one)  # |01>
one_zero = kron(one, zero)  # |10>
one_one = kron(one, one)  # |11>

# Example initial state: a superposition for qubit x, environment in |0>
# psi(0) = (|0> + |1>) / sqrt(2) for x, and |0> for y
psi_0 = (1.0 / np.sqrt(2)) * one + (0.0) * zero  # This is a superposition for x alone
# Combine it with environment y in |0>:
psi_init = kron(psi_0, zero)

# Normalize just in case (should already be normalized):
psi_init = psi_init / np.linalg.norm(psi_init)

rho_init = density_matrix(psi_init)

###############################################################################
# 5. Simulation loop: evolve the state, compute partial trace over environment
###############################################################################
rho_x_history = []
pop_x_0 = []  # Probability to find x in |0>
pop_x_1 = []  # Probability to find x in |1>

rho = rho_init.copy()
for t in times:
    # Partial trace to get rho_x from the full 2-qubit density matrix
    rho_x = partial_trace(rho, subsystem=1)  # trace out qubit y
    rho_x_history.append(rho_x)

    # Probabilities for x in |0> or |1>:
    p0 = np.real_if_close(np.dot(zero.conjugate(), rho_x.dot(zero)))
    p1 = np.real_if_close(np.dot(one.conjugate(), rho_x.dot(one)))
    pop_x_0.append(p0)
    pop_x_1.append(p1)

    # Evolve the full density matrix by dt
    rho = U_dt @ rho @ U_dt.conjugate().T

rho_x_history = np.array(rho_x_history)

###############################################################################
# 6. Plot results
###############################################################################
plt.figure()
plt.plot(times, pop_x_0, label="P(x=0)")
plt.plot(times, pop_x_1, label="P(x=1)")
plt.xlabel("Time")
plt.ylabel("Population of qubit x")
plt.legend()
plt.title("System Qubit Populations Over Time")
plt.show()
