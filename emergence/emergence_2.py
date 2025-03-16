import numpy as np
from numpy import kron
from scipy.linalg import expm
import matplotlib.pyplot as plt

###############################################################################
# 1. Define a larger discrete clock register for finer time resolution
###############################################################################
D_clock = 20  # Increased clock states for smoother evolution
clock_basis = np.eye(D_clock, dtype=complex)

###############################################################################
# 2. Define system + environment: three qubits total
###############################################################################
# Now dim=8 for 3 qubits (2 system + 1 environment)
dim_system = 8

# Basic operators
sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)
sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
I2 = np.eye(2, dtype=complex)

# Build expanded Hamiltonian including environment
omega = 1.0    # system energy scale
g = 0.5        # system-system interaction
g_env = 0.3    # system-environment interaction

# System terms (first two qubits)
H_1 = kron(kron(sigma_z, I2), I2)  # σz on first qubit
H_2 = kron(kron(I2, sigma_z), I2)  # σz on second qubit
H_int = g * kron(kron(sigma_x, sigma_x), I2)  # interaction between system qubits

# Environment coupling (third qubit)
H_env = kron(kron(sigma_z, I2), sigma_z)  # coupling first qubit to environment
H_env += kron(kron(I2, sigma_z), sigma_z)  # coupling second qubit to environment
H_env *= g_env

# Total Hamiltonian
H_sys = omega * (H_1 + H_2) + H_int + H_env

# Time evolution operator (smaller dt for finer resolution)
dt = 0.2
U = expm(-1j * H_sys * dt)

###############################################################################
# 3. Initial state: system in superposition, environment in |0>
###############################################################################
plus = (1.0/np.sqrt(2))*np.array([1, 1], dtype=complex)
zero = np.array([1, 0], dtype=complex)

# Two system qubits in |+>, environment in |0>
psi0_sys = kron(plus, plus)  # system qubits
psi0_total = kron(psi0_sys, zero)  # add environment
psi0_total /= np.linalg.norm(psi0_total)

###############################################################################
# 4. Build global state with clock correlations
###############################################################################
dim_total = D_clock * dim_system
Psi_global = np.zeros(dim_total, dtype=complex)

for T in range(D_clock):
    clock_state = clock_basis[T,:]
    sys_state = np.linalg.matrix_power(U, T).dot(psi0_total)
    kron_state = np.kron(clock_state, sys_state)
    Psi_global += kron_state

Psi_global /= np.linalg.norm(Psi_global)

###############################################################################
# 5. Analysis with reduced density matrices
###############################################################################
Psi_global_reshaped = np.reshape(Psi_global, (D_clock, dim_system))

clock_probs = []
exp_values_Z1 = []  # First system qubit
exp_values_Z2 = []  # Second system qubit
coherence = []      # Off-diagonal element of reduced density matrix

# Observables for first two qubits
Z1 = kron(kron(sigma_z, I2), I2)
Z2 = kron(kron(I2, sigma_z), I2)

def partial_trace_to_qubit(rho_full, keep_qubit):
    """Extract single-qubit density matrix by tracing out others.
    
    Args:
        rho_full: 8x8 density matrix for 3 qubits
        keep_qubit: Which qubit to keep (0, 1, or 2)
    
    Returns:
        2x2 reduced density matrix for the specified qubit
    """
    # Reshape the 8x8 matrix into 2x2x2 x 2x2x2 form
    rho_reshaped = rho_full.reshape([2,2,2, 2,2,2])
    
    if keep_qubit == 0:  # Keep first qubit
        return np.trace(np.trace(rho_reshaped, axis1=1, axis2=4), axis1=1, axis2=2)
    elif keep_qubit == 1:  # Keep second qubit
        return np.trace(np.trace(rho_reshaped, axis1=0, axis2=3), axis1=1, axis2=2)
    else:  # Keep third qubit
        return np.trace(np.trace(rho_reshaped, axis1=0, axis2=3), axis1=0, axis2=1)

for T in range(D_clock):
    # Get system state for this clock time
    block = Psi_global_reshaped[T,:]
    prob_T = np.vdot(block, block).real
    clock_probs.append(prob_T)
    
    if prob_T > 1e-12:
        sys_state = block / np.sqrt(prob_T)
    else:
        sys_state = block
        
    # Full density matrix
    rho_full = np.outer(sys_state, sys_state.conjugate())
    
    # Measure observables
    exp_values_Z1.append(np.real_if_close(np.trace(rho_full @ Z1)))
    exp_values_Z2.append(np.real_if_close(np.trace(rho_full @ Z2)))
    
    # Get reduced density matrix of first qubit and measure coherence
    rho_1 = partial_trace_to_qubit(rho_full, 0)
    coherence.append(np.abs(rho_1[0,1]))

###############################################################################
# 6. Enhanced visualization
###############################################################################
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 12))

# Clock probabilities
ax1.bar(range(D_clock), clock_probs)
ax1.set_xlabel('Clock State T')
ax1.set_ylabel('P(T)')
ax1.set_title('Clock State Probabilities')

# System observables
T_values = np.arange(D_clock)
ax2.plot(T_values, exp_values_Z1, 'b-o', label='⟨σz⊗I⊗I⟩', markersize=4)
ax2.plot(T_values, exp_values_Z2, 'r-o', label='⟨I⊗σz⊗I⟩', markersize=4)
ax2.set_xlabel('Clock Time T')
ax2.set_ylabel('Expectation Value')
ax2.set_title('System Observables vs Clock Time')
ax2.legend()
ax2.grid(True)

# Coherence decay
ax3.plot(T_values, coherence, 'g-o', label='|ρ01|', markersize=4)
ax3.set_xlabel('Clock Time T')
ax3.set_ylabel('Coherence')
ax3.set_title('First Qubit Coherence vs Clock Time')
ax3.grid(True)
ax3.legend()

plt.tight_layout()
plt.show()

# Print key numerical results
print("\nSelected measurements at different clock times:")
for T in [0, D_clock//4, D_clock//2, 3*D_clock//4, D_clock-1]:
    print(f"\nT = {T}:")
    print(f"P(T) = {clock_probs[T]:.4f}")
    print(f"⟨σz⊗I⊗I⟩ = {exp_values_Z1[T]:.4f}")
    print(f"⟨I⊗σz⊗I⟩ = {exp_values_Z2[T]:.4f}")
    print(f"Coherence = {coherence[T]:.4f}")
