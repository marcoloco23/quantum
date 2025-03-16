import numpy as np
from numpy import kron
from scipy.linalg import expm
import matplotlib.pyplot as plt

###############################################################################
# 1. Define a shorter discrete clock register
###############################################################################
D_clock = 10  # Reduced clock states to see early-time dynamics
clock_basis = np.eye(D_clock, dtype=complex)

###############################################################################
# 2. Define system + environment: four qubits total (2 system + 2 environment)
###############################################################################
# Now dim=16 for 4 qubits (2 system + 2 environment)
dim_system = 16

# Basic operators
sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)
sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
I2 = np.eye(2, dtype=complex)

# Build expanded Hamiltonian including environment
omega = 1.0    # system energy scale
g = 0.5        # system-system interaction
g_env1 = 0.02  # first qubit-environment coupling (reduced)
g_env2 = 0.01  # second qubit-environment coupling (reduced)

# Helper function for 4-qubit operators
def four_qubit_op(op1, op2, op3, op4):
    """Create a 4-qubit operator from single-qubit operators."""
    return kron(kron(kron(op1, op2), op3), op4)

# System terms (first two qubits)
H_1 = four_qubit_op(sigma_z, I2, I2, I2)  # σz on first qubit
H_2 = four_qubit_op(I2, sigma_z, I2, I2)  # σz on second qubit
H_int = g * four_qubit_op(sigma_x, sigma_x, I2, I2)  # interaction between system qubits

# Environment coupling (σx-type instead of σz-type)
# First system qubit couples to first environment qubit via σx⊗σx
H_env1 = g_env1 * four_qubit_op(sigma_x, I2, sigma_x, I2)
# Second system qubit couples to second environment qubit via σx⊗σx
H_env2 = g_env2 * four_qubit_op(I2, sigma_x, I2, sigma_x)

# Total Hamiltonian
H_sys = omega * (H_1 + H_2) + H_int + H_env1 + H_env2

# Time evolution operator (smaller dt for finer resolution)
dt = 0.2
U = expm(-1j * H_sys * dt)

###############################################################################
# 3. Initial state: system in superposition, environments in |0>
###############################################################################
plus = (1.0/np.sqrt(2))*np.array([1, 1], dtype=complex)
zero = np.array([1, 0], dtype=complex)

# Two system qubits in |+>, two environment qubits in |0>
psi0_sys = kron(plus, plus)  # system qubits
psi0_total = kron(kron(psi0_sys, zero), zero)  # add two environment qubits
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
coherence1 = []     # First qubit coherence
coherence2 = []     # Second qubit coherence
purity1 = []        # Purity of first qubit (Tr(ρ²))
purity2 = []        # Purity of second qubit (Tr(ρ²))

# Observables for system qubits
Z1 = four_qubit_op(sigma_z, I2, I2, I2)
Z2 = four_qubit_op(I2, sigma_z, I2, I2)

def partial_trace_to_qubit(rho_full, keep_qubit):
    """Extract single-qubit density matrix by tracing out others.
    
    Args:
        rho_full: 16x16 density matrix for 4 qubits
        keep_qubit: Which qubit to keep (0-3)
    
    Returns:
        2x2 reduced density matrix for the specified qubit
    """
    # Reshape the 16x16 matrix into 2x2x2x2 x 2x2x2x2 form
    rho_reshaped = rho_full.reshape([2,2,2,2, 2,2,2,2])
    
    if keep_qubit == 0:  # Keep first qubit
        return np.trace(np.trace(np.trace(rho_reshaped, 
            axis1=1, axis2=5), axis1=1, axis2=3), axis1=1, axis2=2)
    elif keep_qubit == 1:  # Keep second qubit
        return np.trace(np.trace(np.trace(rho_reshaped, 
            axis1=0, axis2=4), axis1=1, axis2=3), axis1=1, axis2=2)
    else:
        raise ValueError("Only implemented for system qubits 0 and 1")

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
    
    # Get reduced density matrices and measure coherence
    rho_1 = partial_trace_to_qubit(rho_full, 0)
    rho_2 = partial_trace_to_qubit(rho_full, 1)
    coherence1.append(np.abs(rho_1[0,1]))
    coherence2.append(np.abs(rho_2[0,1]))
    
    # Calculate purity (Tr(ρ²)) for each qubit
    purity1.append(np.real_if_close(np.trace(rho_1 @ rho_1)))
    purity2.append(np.real_if_close(np.trace(rho_2 @ rho_2)))

###############################################################################
# 6. Enhanced visualization
###############################################################################
fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(10, 16))

# Clock probabilities
ax1.bar(range(D_clock), clock_probs)
ax1.set_xlabel('Clock State T')
ax1.set_ylabel('P(T)')
ax1.set_title('Clock State Probabilities')

# System observables
T_values = np.arange(D_clock)
ax2.plot(T_values, exp_values_Z1, 'b-o', label='⟨σz⊗I⊗I⊗I⟩', markersize=4)
ax2.plot(T_values, exp_values_Z2, 'r-o', label='⟨I⊗σz⊗I⊗I⟩', markersize=4)
ax2.set_xlabel('Clock Time T')
ax2.set_ylabel('Expectation Value')
ax2.set_title('System Observables vs Clock Time')
ax2.legend()
ax2.grid(True)

# Coherence decay
ax3.plot(T_values, coherence1, 'b-o', label='|ρ01| (qubit 1)', markersize=4)
ax3.plot(T_values, coherence2, 'r-o', label='|ρ01| (qubit 2)', markersize=4)
ax3.set_xlabel('Clock Time T')
ax3.set_ylabel('Coherence')
ax3.set_title('Qubit Coherences vs Clock Time')
ax3.grid(True)
ax3.legend()

# Purity
ax4.plot(T_values, purity1, 'b-o', label='Tr(ρ²) (qubit 1)', markersize=4)
ax4.plot(T_values, purity2, 'r-o', label='Tr(ρ²) (qubit 2)', markersize=4)
ax4.set_xlabel('Clock Time T')
ax4.set_ylabel('Purity')
ax4.set_title('Qubit Purity vs Clock Time')
ax4.grid(True)
ax4.legend()

plt.tight_layout()
plt.show()

# Print key numerical results
print("\nSelected measurements at different clock times:")
for T in [0, D_clock//4, D_clock//2, 3*D_clock//4, D_clock-1]:
    print(f"\nT = {T}:")
    print(f"P(T) = {clock_probs[T]:.4f}")
    print(f"⟨σz⊗I⊗I⊗I⟩ = {exp_values_Z1[T]:.4f}")
    print(f"⟨I⊗σz⊗I⊗I⟩ = {exp_values_Z2[T]:.4f}")
    print(f"Coherence 1 = {coherence1[T]:.4f}")
    print(f"Coherence 2 = {coherence2[T]:.4f}")
    print(f"Purity 1 = {purity1[T]:.4f}")
    print(f"Purity 2 = {purity2[T]:.4f}")
