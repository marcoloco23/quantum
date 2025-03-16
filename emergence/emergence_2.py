import numpy as np
from numpy import kron
from scipy.linalg import expm
import matplotlib.pyplot as plt

###############################################################################
# 1. Define a discrete clock register of dimension D_clock
###############################################################################
D_clock = 4  # "time" states: T=0,1,2,3
# We'll label them as basis vectors |0>, |1>, |2>, |3> in clock space.
clock_basis = np.eye(D_clock, dtype=complex)

###############################################################################
# 2. Define the system as two qubits (dim=4)
###############################################################################
dim_system = 4  # 2 qubits => 2^2 = 4

# Basic Pauli matrices and identity
sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)
sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
I2 = np.eye(2, dtype=complex)

# Build system Hamiltonian (similar to emergence.py but simplified)
omega = 1.0  # energy scale
g = 1.0      # interaction strength

H_x = kron(sigma_z, I2)  # acts on first qubit
H_y = kron(I2, sigma_z)  # acts on second qubit
H_int = g * kron(sigma_x, sigma_x)  # interaction term
H_sys = omega * (H_x + H_y) + H_int

# Define single time-step operator U = exp(-i * H_sys * dt)
dt = 0.5
U = expm(-1j * H_sys * dt)

###############################################################################
# 3. Define the initial system state |psi_0> (2-qubit state)
###############################################################################
# Let's start with both qubits in |+> = (|0> + |1>)/sqrt(2)
plus = (1.0/np.sqrt(2))*np.array([1, 1], dtype=complex)
psi0_2qubit = kron(plus, plus)  # 2-qubit product: |+>_1 x |+>_2
psi0_2qubit /= np.linalg.norm(psi0_2qubit)

###############################################################################
# 4. Build the global wavefunction correlating clock states with system states
###############################################################################
# |Psi> = sum_{T=0}^{D_clock-1} |T>_clock ⊗ [U^T |psi0_2qubit>]_system
dim_total = D_clock * dim_system
Psi_global = np.zeros(dim_total, dtype=complex)

for T in range(D_clock):
    # Clock basis vector for time T
    clock_state = clock_basis[T,:]
    # System state after T steps: U^T|psi0>
    sys_state = np.linalg.matrix_power(U, T).dot(psi0_2qubit)
    # Build tensor product |T>_clock ⊗ (U^T|psi0>)
    kron_state = np.kron(clock_state, sys_state)
    Psi_global += kron_state

# Normalize the global state
Psi_global /= np.linalg.norm(Psi_global)

###############################################################################
# 5. Analysis: Measure clock and system observables
###############################################################################
# Reshape global state for analysis
Psi_global_reshaped = np.reshape(Psi_global, (D_clock, dim_system))

# Store results for each clock time T
clock_probs = []  # P(T) - probability of measuring clock = T
system_states = []  # Conditional system state given clock = T
exp_values_Z1 = []  # <σz⊗I> for each T
exp_values_Z2 = []  # <I⊗σz> for each T

# Observables for measurements
Z1 = kron(sigma_z, I2)  # σz on first qubit
Z2 = kron(I2, sigma_z)  # σz on second qubit

# Analyze state for each possible clock reading T
for T in range(D_clock):
    # Extract system state conditioned on clock = T
    block = Psi_global_reshaped[T,:]
    prob_T = np.vdot(block, block).real
    clock_probs.append(prob_T)
    
    # Normalize the conditional system state
    if prob_T > 1e-12:
        sys_state = block / np.sqrt(prob_T)
    else:
        sys_state = block
    system_states.append(sys_state)
    
    # Measure observables in this state
    rho_sys = np.outer(sys_state, sys_state.conjugate())
    exp_values_Z1.append(np.real_if_close(np.trace(rho_sys @ Z1)))
    exp_values_Z2.append(np.real_if_close(np.trace(rho_sys @ Z2)))

###############################################################################
# 6. Visualization
###############################################################################
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

# Plot clock state probabilities
ax1.bar(range(D_clock), clock_probs)
ax1.set_xlabel('Clock State T')
ax1.set_ylabel('P(T)')
ax1.set_title('Clock State Probabilities')

# Plot expectation values of observables
T_values = np.arange(D_clock)
ax2.plot(T_values, exp_values_Z1, 'b-o', label='⟨σz⊗I⟩')
ax2.plot(T_values, exp_values_Z2, 'r-o', label='⟨I⊗σz⟩')
ax2.set_xlabel('Clock Time T')
ax2.set_ylabel('Expectation Value')
ax2.set_title('System Observables vs Clock Time')
ax2.legend()
ax2.grid(True)

plt.tight_layout()
plt.show()

# Print numerical results
print("\nClock probabilities and system measurements:")
for T in range(D_clock):
    print(f"\nT = {T}:")
    print(f"P(T) = {clock_probs[T]:.4f}")
    print(f"⟨σz⊗I⟩ = {exp_values_Z1[T]:.4f}")
    print(f"⟨I⊗σz⟩ = {exp_values_Z2[T]:.4f}")
