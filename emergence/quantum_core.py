import numpy as np
from numpy import kron
from scipy.linalg import expm
from functools import reduce

"""
Quantum core module providing fundamental quantum mechanics utilities.
"""

###############################################################################
# Basic quantum operators and states
###############################################################################
# Pauli matrices and identity
sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)
I2 = np.eye(2, dtype=complex)

# Common quantum states
zero = np.array([1, 0], dtype=complex)
one = np.array([0, 1], dtype=complex)
plus = (1.0 / np.sqrt(2)) * np.array([1, 1], dtype=complex)
minus = (1.0 / np.sqrt(2)) * np.array([1, -1], dtype=complex)

# Bell states
bell_plus = (1.0 / np.sqrt(2)) * np.array(
    [1, 0, 0, 1], dtype=complex
)  # (|00⟩ + |11⟩)/√2
bell_minus = (1.0 / np.sqrt(2)) * np.array(
    [1, 0, 0, -1], dtype=complex
)  # (|00⟩ - |11⟩)/√2


###############################################################################
# Quantum operator utilities
###############################################################################
def multi_qubit_op(op_list):
    """Create a multi-qubit operator from a list of single-qubit operators.

    Args:
        op_list: List of operators [op_1, op_2, ..., op_n]
                 where each op_i acts on the i-th qubit

    Returns:
        Tensor product of all operators
    """
    return reduce(np.kron, op_list)


def four_qubit_op(op1, op2, op3, op4):
    """Create a 4-qubit operator from single-qubit operators.

    Args:
        op1: Operator for the first qubit
        op2: Operator for the second qubit
        op3: Operator for the third qubit
        op4: Operator for the fourth qubit

    Returns:
        4-qubit operator as tensor product
    """
    return kron(kron(kron(op1, op2), op3), op4)


###############################################################################
# Density matrix and measurement utilities
###############################################################################
def density_matrix(state):
    """Convert a pure state vector to density matrix.

    Args:
        state: Pure quantum state as a complex vector

    Returns:
        Density matrix representation of the state
    """
    return np.outer(state, state.conjugate())


def partial_trace_to_qubit(rho_full, keep_qubit, n_qubits=4):
    """Extract single-qubit density matrix by tracing out others.

    Args:
        rho_full: Full density matrix for n_qubits
        keep_qubit: Index of the qubit to keep (0-based)
        n_qubits: Total number of qubits in the system

    Returns:
        2x2 density matrix for the kept qubit
    """
    if n_qubits == 4:
        # Specialized implementation for 4-qubit system
        # Reshape the 16x16 matrix into 2x2x2x2 x 2x2x2x2 form
        rho_reshaped = rho_full.reshape([2, 2, 2, 2, 2, 2, 2, 2])

        if keep_qubit == 0:  # Keep first qubit
            return np.trace(
                np.trace(np.trace(rho_reshaped, axis1=1, axis2=5), axis1=1, axis2=3),
                axis1=1,
                axis2=2,
            )
        elif keep_qubit == 1:  # Keep second qubit
            return np.trace(
                np.trace(np.trace(rho_reshaped, axis1=0, axis2=4), axis1=1, axis2=3),
                axis1=1,
                axis2=2,
            )
        else:
            raise ValueError(
                "For 4-qubit system, only implemented for system qubits 0 and 1"
            )
    else:
        # Generic implementation for arbitrary qubit systems
        # Only implemented for first two qubits (system qubits)
        if keep_qubit > 1:
            raise ValueError("Generic implementation only supports qubits 0 and 1")

        # Dimensions
        dim_rest = 2 ** (n_qubits - 1)

        # Initialize the reduced density matrix
        rho_qubit = np.zeros((2, 2), dtype=complex)

        # Perform partial trace by direct summation
        for i in range(2):
            for j in range(2):
                for k in range(dim_rest):
                    if keep_qubit == 0:
                        # First qubit: stride is 2*dim_rest
                        row = i * dim_rest + k
                        col = j * dim_rest + k
                    else:  # keep_qubit == 1
                        # Second qubit: stride is dim_rest/2
                        block_size = dim_rest // 2
                        block = k // block_size
                        offset = k % block_size
                        row = block * 2 * block_size + i * block_size + offset
                        col = block * 2 * block_size + j * block_size + offset

                    rho_qubit[i, j] += rho_full[row, col]

        return rho_qubit


def get_system_density_matrix(full_state, n_qubits=4, n_system_qubits=2):
    """Get reduced density matrix for the system qubits.

    Args:
        full_state: State vector of the full system
        n_qubits: Total number of qubits
        n_system_qubits: Number of system qubits (assumed to be first qubits)

    Returns:
        Reduced density matrix for the system qubits
    """
    if n_qubits == 4 and n_system_qubits == 2:
        # Specialized implementation for 2 system qubits in a 4-qubit system
        rho_full = density_matrix(full_state)
        # Reshape to 2x2x2x2 x 2x2x2x2 form
        rho_reshaped = rho_full.reshape([2, 2, 2, 2, 2, 2, 2, 2])
        # Trace out environment qubits (last two)
        # First trace out the 4th qubit (indices 3 and 7)
        rho_traced_1 = np.trace(rho_reshaped, axis1=3, axis2=7)
        # Then trace out the 3rd qubit (now indices 2 and 5 after first trace)
        rho_traced_2 = np.trace(rho_traced_1, axis1=2, axis2=5)
        return rho_traced_2
    else:
        # Generic implementation for arbitrary qubit systems
        rho_full = density_matrix(full_state)

        # Dimensions
        dim_sys = 2**n_system_qubits
        dim_env = 2 ** (n_qubits - n_system_qubits)

        # Initialize the reduced density matrix
        rho_sys = np.zeros((dim_sys, dim_sys), dtype=complex)

        # Perform partial trace by direct summation
        for i in range(dim_sys):
            for j in range(dim_sys):
                for k in range(dim_env):
                    row = i * dim_env + k
                    col = j * dim_env + k
                    rho_sys[i, j] += rho_full[row, col]

        return rho_sys


def entanglement_entropy(rho):
    """Compute von Neumann entropy S = -Tr(ρ ln ρ).

    Args:
        rho: Density matrix

    Returns:
        Entanglement entropy
    """
    # Get eigenvalues of density matrix
    eigenvalues = np.linalg.eigvalsh(rho)
    # Remove very small eigenvalues (numerical stability)
    eigenvalues = eigenvalues[eigenvalues > 1e-10]

    # Handle case where there might be no valid eigenvalues
    if len(eigenvalues) == 0:
        return 0.0

    # Compute entropy
    entropy = -np.sum(eigenvalues * np.log(eigenvalues))

    # Return real part (should be real anyway, but ensure no numerical artifacts)
    return np.real_if_close(entropy)


def purity(rho):
    """Calculate purity of a quantum state: Tr(ρ²).

    Args:
        rho: Density matrix

    Returns:
        Purity value
    """
    return np.real_if_close(np.trace(rho @ rho))


def coherence(rho):
    """Calculate the coherence |ρ01| of a qubit.

    Args:
        rho: 2x2 density matrix

    Returns:
        Coherence value
    """
    return np.abs(rho[0, 1])


def expectation_value(rho, operator):
    """Calculate the expectation value of an operator.

    Args:
        rho: Density matrix
        operator: Quantum operator

    Returns:
        Expectation value
    """
    return np.real_if_close(np.trace(rho @ operator))


def time_evolution_operator(hamiltonian, time_step):
    """Generate time evolution operator U = exp(-i*H*dt).

    Args:
        hamiltonian: Hamiltonian matrix
        time_step: Time step for evolution

    Returns:
        Time evolution operator
    """
    return expm(-1j * hamiltonian * time_step)
