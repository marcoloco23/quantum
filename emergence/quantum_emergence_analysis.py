import numpy as np
from numpy import kron
from scipy.linalg import expm
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

###############################################################################
# 1. Define system parameters and helper functions
###############################################################################
# Clock parameters
D_clock = 10  # Number of clock states

# Basic operators
sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)
sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
I2 = np.eye(2, dtype=complex)

# Helper function for 4-qubit operators
def four_qubit_op(op1, op2, op3, op4):
    """Create a 4-qubit operator from single-qubit operators."""
    return kron(kron(kron(op1, op2), op3), op4)

# System Hamiltonian components (fixed)
omega = 1.0    # system energy scale
g = 0.5        # system-system interaction
dt = 0.2       # time step

# System terms (first two qubits)
H_1 = four_qubit_op(sigma_z, I2, I2, I2)  # σz on first qubit
H_2 = four_qubit_op(I2, sigma_z, I2, I2)  # σz on second qubit
H_int = g * four_qubit_op(sigma_x, sigma_x, I2, I2)  # interaction between system qubits

# Initial state preparation
plus = (1.0/np.sqrt(2))*np.array([1, 1], dtype=complex)
zero = np.array([1, 0], dtype=complex)
psi0_sys = kron(plus, plus)  # system qubits in |+>

# Partial trace function
def partial_trace_to_qubit(rho_full, keep_qubit):
    """Extract single-qubit density matrix by tracing out others."""
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

# Function to compute entanglement entropy
def entanglement_entropy(rho):
    """Compute von Neumann entropy S = -Tr(ρ ln ρ)."""
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

# Function to trace out environment and get system density matrix
def get_system_density_matrix(full_state):
    """Get 4x4 density matrix for the two system qubits."""
    # Full state is 16-dimensional
    rho_full = np.outer(full_state, full_state.conjugate())
    # Reshape to 2x2x2x2 x 2x2x2x2 form
    rho_reshaped = rho_full.reshape([2,2,2,2, 2,2,2,2])
    # Trace out environment qubits (last two)
    # First trace out the 4th qubit (indices 3 and 7)
    rho_traced_1 = np.trace(rho_reshaped, axis1=3, axis2=7)
    # Then trace out the 3rd qubit (now indices 2 and 5 after first trace)
    rho_traced_2 = np.trace(rho_traced_1, axis1=2, axis2=5)
    return rho_traced_2

###############################################################################
# 2. Simulation function for a single parameter set
###############################################################################
def run_simulation(g_env1, g_env2, coupling_type="xx"):
    """Run simulation with given environment coupling strengths and type."""
    # Special case: no environment coupling
    if g_env1 == 0.0 and g_env2 == 0.0:
        # Just system Hamiltonian
        H_sys = omega * (H_1 + H_2) + H_int
    else:
        # Environment coupling
        if coupling_type == "xx":
            # σx⊗σx coupling
            H_env1 = g_env1 * four_qubit_op(sigma_x, I2, sigma_x, I2)
            H_env2 = g_env2 * four_qubit_op(I2, sigma_x, I2, sigma_x)
        elif coupling_type == "zz":
            # σz⊗σz coupling
            H_env1 = g_env1 * four_qubit_op(sigma_z, I2, sigma_z, I2)
            H_env2 = g_env2 * four_qubit_op(I2, sigma_z, I2, sigma_z)
        elif coupling_type == "mixed":
            # Mixed σx⊗σx and σz⊗σz coupling
            H_env1 = g_env1 * (0.5 * four_qubit_op(sigma_x, I2, sigma_x, I2) + 
                            0.5 * four_qubit_op(sigma_z, I2, sigma_z, I2))
            H_env2 = g_env2 * (0.5 * four_qubit_op(I2, sigma_x, I2, sigma_x) + 
                            0.5 * four_qubit_op(I2, sigma_z, I2, sigma_z))
        else:
            raise ValueError(f"Unknown coupling type: {coupling_type}")
        
        # Total Hamiltonian
        H_sys = omega * (H_1 + H_2) + H_int + H_env1 + H_env2
    
    # Time evolution operator
    U = expm(-1j * H_sys * dt)
    
    # Initial state: system in |+>⊗|+>, environment in |0>⊗|0>
    psi0_total = kron(kron(psi0_sys, zero), zero)
    psi0_total /= np.linalg.norm(psi0_total)
    
    # Clock basis
    clock_basis = np.eye(D_clock, dtype=complex)
    
    # Build global state
    dim_system = 16  # 4 qubits
    dim_total = D_clock * dim_system
    Psi_global = np.zeros(dim_total, dtype=complex)
    
    for T in range(D_clock):
        clock_state = clock_basis[T,:]
        sys_state = np.linalg.matrix_power(U, T).dot(psi0_total)
        kron_state = np.kron(clock_state, sys_state)
        Psi_global += kron_state
    
    Psi_global /= np.linalg.norm(Psi_global)
    
    # Reshape for analysis
    Psi_global_reshaped = np.reshape(Psi_global, (D_clock, dim_system))
    
    # Initialize result arrays
    clock_probs = []
    exp_values_Z1 = []
    exp_values_Z2 = []
    coherence1 = []
    coherence2 = []
    purity1 = []
    purity2 = []
    entropy = []
    
    # Observables
    Z1 = four_qubit_op(sigma_z, I2, I2, I2)
    Z2 = four_qubit_op(I2, sigma_z, I2, I2)
    
    # Analyze each clock time
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
        
        # Calculate purity for each qubit
        purity1.append(np.real_if_close(np.trace(rho_1 @ rho_1)))
        purity2.append(np.real_if_close(np.trace(rho_2 @ rho_2)))
        
        # Calculate system-environment entanglement entropy
        rho_sys = get_system_density_matrix(sys_state)
        try:
            entropy_val = entanglement_entropy(rho_sys)
            entropy.append(entropy_val)
        except Exception as e:
            print(f"Warning: Error calculating entropy at T={T}: {e}")
            entropy.append(0.0)
    
    # Return all results
    return {
        'clock_probs': np.array(clock_probs),
        'exp_values_Z1': np.array(exp_values_Z1),
        'exp_values_Z2': np.array(exp_values_Z2),
        'coherence1': np.array(coherence1),
        'coherence2': np.array(coherence2),
        'purity1': np.array(purity1),
        'purity2': np.array(purity2),
        'entropy': np.array(entropy)
    }

###############################################################################
# 3. Parameter sweep and analysis
###############################################################################
def coupling_strength_sweep():
    """Sweep environment coupling strengths and analyze results."""
    # Define parameter ranges
    g_env1_values = np.linspace(0.0, 0.25, 6)  # 0.0, 0.05, 0.1, ..., 0.25
    g_env2_values = np.linspace(0.0, 0.25, 6)  # 0.0, 0.05, 0.1, ..., 0.25
    
    # Initialize result arrays
    avg_coherence1 = np.zeros((len(g_env1_values), len(g_env2_values)))
    avg_coherence2 = np.zeros((len(g_env1_values), len(g_env2_values)))
    avg_purity1 = np.zeros((len(g_env1_values), len(g_env2_values)))
    avg_purity2 = np.zeros((len(g_env1_values), len(g_env2_values)))
    avg_entropy = np.zeros((len(g_env1_values), len(g_env2_values)))
    
    # Run simulations for each parameter combination
    for i, g1 in enumerate(g_env1_values):
        for j, g2 in enumerate(g_env2_values):
            print(f"Running simulation with g_env1={g1:.2f}, g_env2={g2:.2f}")
            results = run_simulation(g1, g2, coupling_type="xx")
            
            # Store average values
            avg_coherence1[i, j] = np.mean(results['coherence1'])
            avg_coherence2[i, j] = np.mean(results['coherence2'])
            avg_purity1[i, j] = np.mean(results['purity1'])
            avg_purity2[i, j] = np.mean(results['purity2'])
            avg_entropy[i, j] = np.mean(results['entropy'])
    
    # Create plots
    fig = plt.figure(figsize=(18, 12))
    
    # 3D surface plot for average coherence
    ax1 = fig.add_subplot(231, projection='3d')
    X, Y = np.meshgrid(g_env1_values, g_env2_values)
    surf = ax1.plot_surface(X, Y, avg_coherence1.T, cmap=cm.viridis)
    ax1.set_xlabel('g_env1')
    ax1.set_ylabel('g_env2')
    ax1.set_zlabel('Avg Coherence (Qubit 1)')
    ax1.set_title('Average Coherence vs Coupling Strengths')
    fig.colorbar(surf, ax=ax1, shrink=0.5)
    
    # 3D surface plot for average purity
    ax2 = fig.add_subplot(232, projection='3d')
    surf = ax2.plot_surface(X, Y, avg_purity1.T, cmap=cm.viridis)
    ax2.set_xlabel('g_env1')
    ax2.set_ylabel('g_env2')
    ax2.set_zlabel('Avg Purity (Qubit 1)')
    ax2.set_title('Average Purity vs Coupling Strengths')
    fig.colorbar(surf, ax=ax2, shrink=0.5)
    
    # 3D surface plot for average entropy
    ax3 = fig.add_subplot(233, projection='3d')
    surf = ax3.plot_surface(X, Y, avg_entropy.T, cmap=cm.viridis)
    ax3.set_xlabel('g_env1')
    ax3.set_ylabel('g_env2')
    ax3.set_zlabel('Avg Entropy')
    ax3.set_title('Average Entanglement Entropy vs Coupling Strengths')
    fig.colorbar(surf, ax=ax3, shrink=0.5)
    
    # 2D heatmap for average coherence
    ax4 = fig.add_subplot(234)
    im = ax4.imshow(avg_coherence1, origin='lower', extent=[g_env1_values[0], g_env1_values[-1], 
                                                          g_env2_values[0], g_env2_values[-1]])
    ax4.set_xlabel('g_env1')
    ax4.set_ylabel('g_env2')
    ax4.set_title('Average Coherence Heatmap (Qubit 1)')
    fig.colorbar(im, ax=ax4)
    
    # 2D heatmap for average purity
    ax5 = fig.add_subplot(235)
    im = ax5.imshow(avg_purity1, origin='lower', extent=[g_env1_values[0], g_env1_values[-1], 
                                                       g_env2_values[0], g_env2_values[-1]])
    ax5.set_xlabel('g_env1')
    ax5.set_ylabel('g_env2')
    ax5.set_title('Average Purity Heatmap (Qubit 1)')
    fig.colorbar(im, ax=ax5)
    
    # 2D heatmap for average entropy
    ax6 = fig.add_subplot(236)
    im = ax6.imshow(avg_entropy, origin='lower', extent=[g_env1_values[0], g_env1_values[-1], 
                                                       g_env2_values[0], g_env2_values[-1]])
    ax6.set_xlabel('g_env1')
    ax6.set_ylabel('g_env2')
    ax6.set_title('Average Entropy Heatmap')
    fig.colorbar(im, ax=ax6)
    
    plt.tight_layout()
    plt.savefig('images/coupling_sweep_results.png', dpi=300)
    plt.show()
    
    return {
        'g_env1_values': g_env1_values,
        'g_env2_values': g_env2_values,
        'avg_coherence1': avg_coherence1,
        'avg_coherence2': avg_coherence2,
        'avg_purity1': avg_purity1,
        'avg_purity2': avg_purity2,
        'avg_entropy': avg_entropy
    }

###############################################################################
# 4. Compare different coupling types
###############################################################################
def compare_coupling_types():
    """Compare xx, zz, and mixed coupling types."""
    # Fixed coupling strengths
    g_env1 = 0.05
    g_env2 = 0.05
    
    # Run simulations for each coupling type
    results_xx = run_simulation(g_env1, g_env2, coupling_type="xx")
    results_zz = run_simulation(g_env1, g_env2, coupling_type="zz")
    results_mixed = run_simulation(g_env1, g_env2, coupling_type="mixed")
    
    # Create plots
    fig, axes = plt.subplots(3, 1, figsize=(10, 12))
    
    # Time values
    T_values = np.arange(D_clock)
    
    # Plot coherence
    ax1 = axes[0]
    ax1.plot(T_values, results_xx['coherence1'], 'b-o', label='XX Coupling', markersize=4)
    ax1.plot(T_values, results_zz['coherence1'], 'r-o', label='ZZ Coupling', markersize=4)
    ax1.plot(T_values, results_mixed['coherence1'], 'g-o', label='Mixed Coupling', markersize=4)
    ax1.set_xlabel('Clock Time T')
    ax1.set_ylabel('Coherence |ρ01|')
    ax1.set_title(f'Qubit 1 Coherence vs Clock Time (g_env={g_env1:.2f})')
    ax1.grid(True)
    ax1.legend()
    
    # Plot purity
    ax2 = axes[1]
    ax2.plot(T_values, results_xx['purity1'], 'b-o', label='XX Coupling', markersize=4)
    ax2.plot(T_values, results_zz['purity1'], 'r-o', label='ZZ Coupling', markersize=4)
    ax2.plot(T_values, results_mixed['purity1'], 'g-o', label='Mixed Coupling', markersize=4)
    ax2.set_xlabel('Clock Time T')
    ax2.set_ylabel('Purity Tr(ρ²)')
    ax2.set_title(f'Qubit 1 Purity vs Clock Time (g_env={g_env1:.2f})')
    ax2.grid(True)
    ax2.legend()
    
    # Plot entropy
    ax3 = axes[2]
    ax3.plot(T_values, results_xx['entropy'], 'b-o', label='XX Coupling', markersize=4)
    ax3.plot(T_values, results_zz['entropy'], 'r-o', label='ZZ Coupling', markersize=4)
    ax3.plot(T_values, results_mixed['entropy'], 'g-o', label='Mixed Coupling', markersize=4)
    ax3.set_xlabel('Clock Time T')
    ax3.set_ylabel('Entanglement Entropy')
    ax3.set_title(f'System-Environment Entanglement Entropy (g_env={g_env1:.2f})')
    ax3.grid(True)
    ax3.legend()
    
    plt.tight_layout()
    plt.savefig('images/coupling_type_comparison.png', dpi=300)
    plt.show()

###############################################################################
# 5. Main execution
###############################################################################
if __name__ == "__main__":
    print("Emergent Time with Environment Coupling Analysis")
    print("1. Coupling Strength Sweep")
    print("2. Compare Coupling Types")
    
    choice = input("Enter your choice (1 or 2): ")
    
    if choice == "1":
        print("Running coupling strength sweep (this may take a while)...")
        results = coupling_strength_sweep()
        print("Sweep complete. Results saved to 'images/coupling_sweep_results.png'")
    elif choice == "2":
        print("Comparing different coupling types...")
        compare_coupling_types()
        print("Comparison complete. Results saved to 'images/coupling_type_comparison.png'")
    else:
        print("Invalid choice. Please run again and select 1 or 2.")
