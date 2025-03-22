import numpy as np
from numpy import kron
from scipy.linalg import expm
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from scipy.fft import fft, fftfreq  # Add FFT imports
from scipy.signal import correlate  # Add for autocorrelation
from functools import reduce  # For larger environment simulation

###############################################################################
# 1. Define system parameters and helper functions
###############################################################################
# Clock parameters
D_clock = 100  # Number of clock states

# Basic operators
sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)
sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
I2 = np.eye(2, dtype=complex)


# Helper function for 4-qubit operators
def four_qubit_op(op1, op2, op3, op4):
    """Create a 4-qubit operator from single-qubit operators."""
    return kron(kron(kron(op1, op2), op3), op4)


# System Hamiltonian components (fixed)
omega = 1.0  # system energy scale
g = 0.5  # system-system interaction
dt = 0.2  # time step

# System terms (first two qubits)
H_1 = four_qubit_op(sigma_z, I2, I2, I2)  # σz on first qubit
H_2 = four_qubit_op(I2, sigma_z, I2, I2)  # σz on second qubit
H_int = g * four_qubit_op(sigma_x, sigma_x, I2, I2)  # interaction between system qubits

# Initial state preparation
plus = (1.0 / np.sqrt(2)) * np.array([1, 1], dtype=complex)
zero = np.array([1, 0], dtype=complex)
psi0_sys = kron(plus, plus)  # system qubits in |+>


# Partial trace function
def partial_trace_to_qubit(rho_full, keep_qubit):
    """Extract single-qubit density matrix by tracing out others."""
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
    rho_reshaped = rho_full.reshape([2, 2, 2, 2, 2, 2, 2, 2])
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
            H_env1 = g_env1 * (
                0.5 * four_qubit_op(sigma_x, I2, sigma_x, I2)
                + 0.5 * four_qubit_op(sigma_z, I2, sigma_z, I2)
            )
            H_env2 = g_env2 * (
                0.5 * four_qubit_op(I2, sigma_x, I2, sigma_x)
                + 0.5 * four_qubit_op(I2, sigma_z, I2, sigma_z)
            )
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
        clock_state = clock_basis[T, :]
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
        block = Psi_global_reshaped[T, :]
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
        coherence1.append(np.abs(rho_1[0, 1]))
        coherence2.append(np.abs(rho_2[0, 1]))

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
        "clock_probs": np.array(clock_probs),
        "exp_values_Z1": np.array(exp_values_Z1),
        "exp_values_Z2": np.array(exp_values_Z2),
        "coherence1": np.array(coherence1),
        "coherence2": np.array(coherence2),
        "purity1": np.array(purity1),
        "purity2": np.array(purity2),
        "entropy": np.array(entropy),
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
            avg_coherence1[i, j] = np.mean(results["coherence1"])
            avg_coherence2[i, j] = np.mean(results["coherence2"])
            avg_purity1[i, j] = np.mean(results["purity1"])
            avg_purity2[i, j] = np.mean(results["purity2"])
            avg_entropy[i, j] = np.mean(results["entropy"])

    # Create plots
    fig = plt.figure(figsize=(18, 12))

    # 3D surface plot for average coherence
    ax1 = fig.add_subplot(231, projection="3d")
    X, Y = np.meshgrid(g_env1_values, g_env2_values)
    surf = ax1.plot_surface(X, Y, avg_coherence1.T, cmap=cm.viridis)
    ax1.set_xlabel("g_env1")
    ax1.set_ylabel("g_env2")
    ax1.set_zlabel("Avg Coherence (Qubit 1)")
    ax1.set_title("Average Coherence vs Coupling Strengths")
    fig.colorbar(surf, ax=ax1, shrink=0.5)

    # 3D surface plot for average purity
    ax2 = fig.add_subplot(232, projection="3d")
    surf = ax2.plot_surface(X, Y, avg_purity1.T, cmap=cm.viridis)
    ax2.set_xlabel("g_env1")
    ax2.set_ylabel("g_env2")
    ax2.set_zlabel("Avg Purity (Qubit 1)")
    ax2.set_title("Average Purity vs Coupling Strengths")
    fig.colorbar(surf, ax=ax2, shrink=0.5)

    # 3D surface plot for average entropy
    ax3 = fig.add_subplot(233, projection="3d")
    surf = ax3.plot_surface(X, Y, avg_entropy.T, cmap=cm.viridis)
    ax3.set_xlabel("g_env1")
    ax3.set_ylabel("g_env2")
    ax3.set_zlabel("Avg Entropy")
    ax3.set_title("Average Entanglement Entropy vs Coupling Strengths")
    fig.colorbar(surf, ax=ax3, shrink=0.5)

    # 2D heatmap for average coherence
    ax4 = fig.add_subplot(234)
    im = ax4.imshow(
        avg_coherence1,
        origin="lower",
        extent=[
            g_env1_values[0],
            g_env1_values[-1],
            g_env2_values[0],
            g_env2_values[-1],
        ],
    )
    ax4.set_xlabel("g_env1")
    ax4.set_ylabel("g_env2")
    ax4.set_title("Average Coherence Heatmap (Qubit 1)")
    fig.colorbar(im, ax=ax4)

    # 2D heatmap for average purity
    ax5 = fig.add_subplot(235)
    im = ax5.imshow(
        avg_purity1,
        origin="lower",
        extent=[
            g_env1_values[0],
            g_env1_values[-1],
            g_env2_values[0],
            g_env2_values[-1],
        ],
    )
    ax5.set_xlabel("g_env1")
    ax5.set_ylabel("g_env2")
    ax5.set_title("Average Purity Heatmap (Qubit 1)")
    fig.colorbar(im, ax=ax5)

    # 2D heatmap for average entropy
    ax6 = fig.add_subplot(236)
    im = ax6.imshow(
        avg_entropy,
        origin="lower",
        extent=[
            g_env1_values[0],
            g_env1_values[-1],
            g_env2_values[0],
            g_env2_values[-1],
        ],
    )
    ax6.set_xlabel("g_env1")
    ax6.set_ylabel("g_env2")
    ax6.set_title("Average Entropy Heatmap")
    fig.colorbar(im, ax=ax6)

    plt.tight_layout()
    plt.savefig("images/coupling_sweep_results.png", dpi=300)
    plt.show()

    return {
        "g_env1_values": g_env1_values,
        "g_env2_values": g_env2_values,
        "avg_coherence1": avg_coherence1,
        "avg_coherence2": avg_coherence2,
        "avg_purity1": avg_purity1,
        "avg_purity2": avg_purity2,
        "avg_entropy": avg_entropy,
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
    ax1.plot(
        T_values, results_xx["coherence1"], "b-o", label="XX Coupling", markersize=4
    )
    ax1.plot(
        T_values, results_zz["coherence1"], "r-o", label="ZZ Coupling", markersize=4
    )
    ax1.plot(
        T_values,
        results_mixed["coherence1"],
        "g-o",
        label="Mixed Coupling",
        markersize=4,
    )
    ax1.set_xlabel("Clock Time T")
    ax1.set_ylabel("Coherence |ρ01|")
    ax1.set_title(f"Qubit 1 Coherence vs Clock Time (g_env={g_env1:.2f})")
    ax1.grid(True)
    ax1.legend()

    # Plot purity
    ax2 = axes[1]
    ax2.plot(T_values, results_xx["purity1"], "b-o", label="XX Coupling", markersize=4)
    ax2.plot(T_values, results_zz["purity1"], "r-o", label="ZZ Coupling", markersize=4)
    ax2.plot(
        T_values, results_mixed["purity1"], "g-o", label="Mixed Coupling", markersize=4
    )
    ax2.set_xlabel("Clock Time T")
    ax2.set_ylabel("Purity Tr(ρ²)")
    ax2.set_title(f"Qubit 1 Purity vs Clock Time (g_env={g_env1:.2f})")
    ax2.grid(True)
    ax2.legend()

    # Plot entropy
    ax3 = axes[2]
    ax3.plot(T_values, results_xx["entropy"], "b-o", label="XX Coupling", markersize=4)
    ax3.plot(T_values, results_zz["entropy"], "r-o", label="ZZ Coupling", markersize=4)
    ax3.plot(
        T_values, results_mixed["entropy"], "g-o", label="Mixed Coupling", markersize=4
    )
    ax3.set_xlabel("Clock Time T")
    ax3.set_ylabel("Entanglement Entropy")
    ax3.set_title(f"System-Environment Entanglement Entropy (g_env={g_env1:.2f})")
    ax3.grid(True)
    ax3.legend()

    plt.tight_layout()
    plt.savefig("images/coupling_type_comparison.png", dpi=300)
    plt.show()

    # Frequency analysis of coherence and purity signals
    if D_clock >= 50:  # Only perform if we have enough data points
        analyze_frequency_components(results_xx, results_zz, results_mixed)


###############################################################################
# 5. Frequency domain analysis
###############################################################################
def analyze_frequency_components(results_xx, results_zz, results_mixed):
    """Analyze frequency components of coherence and purity signals."""
    # Sample spacing (assuming uniform time steps)
    sample_spacing = dt

    # Time domain signals
    coherence_xx = results_xx["coherence1"]
    coherence_mixed = results_mixed["coherence1"]
    purity_xx = results_xx["purity1"]
    purity_zz = results_zz["purity1"]
    purity_mixed = results_mixed["purity1"]

    # Compute FFT
    fft_coherence_xx = fft(coherence_xx)
    fft_coherence_mixed = fft(coherence_mixed)
    fft_purity_xx = fft(purity_xx)
    fft_purity_zz = fft(purity_zz)
    fft_purity_mixed = fft(purity_mixed)

    # Compute frequency axis
    n = len(coherence_xx)
    freq = fftfreq(n, d=sample_spacing)

    # Only plot positive frequencies up to Nyquist frequency
    positive_freq_idx = np.arange(1, n // 2)

    # Create plots
    fig, axes = plt.subplots(2, 1, figsize=(12, 10))

    # Plot coherence frequency components
    ax1 = axes[0]
    ax1.plot(
        freq[positive_freq_idx],
        2.0 / n * np.abs(fft_coherence_xx[positive_freq_idx]),
        "b-",
        label="XX Coupling",
    )
    ax1.plot(
        freq[positive_freq_idx],
        2.0 / n * np.abs(fft_coherence_mixed[positive_freq_idx]),
        "g-",
        label="Mixed Coupling",
    )
    ax1.set_xlabel("Frequency (1/time)")
    ax1.set_ylabel("Amplitude")
    ax1.set_title("Frequency Components of Coherence Signal")
    ax1.grid(True)
    ax1.legend()

    # Plot purity frequency components
    ax2 = axes[1]
    ax2.plot(
        freq[positive_freq_idx],
        2.0 / n * np.abs(fft_purity_xx[positive_freq_idx]),
        "b-",
        label="XX Coupling",
    )
    ax2.plot(
        freq[positive_freq_idx],
        2.0 / n * np.abs(fft_purity_zz[positive_freq_idx]),
        "r-",
        label="ZZ Coupling",
    )
    ax2.plot(
        freq[positive_freq_idx],
        2.0 / n * np.abs(fft_purity_mixed[positive_freq_idx]),
        "g-",
        label="Mixed Coupling",
    )
    ax2.set_xlabel("Frequency (1/time)")
    ax2.set_ylabel("Amplitude")
    ax2.set_title("Frequency Components of Purity Signal")
    ax2.grid(True)
    ax2.legend()

    plt.tight_layout()
    plt.savefig("images/frequency_analysis.png", dpi=300)
    plt.show()


###############################################################################
# 6. Recurrence and correlation analysis
###############################################################################
def analyze_recurrence_dynamics(results_xx, results_zz, results_mixed):
    """Analyze recurrence dynamics and temporal correlations."""
    # Time axis
    T_values = np.arange(D_clock)

    # Create figure
    fig, axes = plt.subplots(3, 2, figsize=(15, 12))

    # Extract time series data
    coherence_xx = results_xx["coherence1"]
    coherence_zz = results_zz["coherence1"]
    coherence_mixed = results_mixed["coherence1"]

    purity_xx = results_xx["purity1"]
    purity_zz = results_zz["purity1"]
    purity_mixed = results_mixed["purity1"]

    entropy_xx = results_xx["entropy"]
    entropy_zz = results_zz["entropy"]
    entropy_mixed = results_mixed["entropy"]

    # Compute auto-correlation functions
    # Normalize to have auto-correlation=1 at lag=0
    def normalized_autocorr(x):
        result = correlate(x - np.mean(x), x - np.mean(x), mode="full")
        return result[len(result) // 2 :] / result[len(result) // 2]

    # Coherence auto-correlations
    ac_coherence_xx = normalized_autocorr(coherence_xx)
    ac_coherence_mixed = normalized_autocorr(coherence_mixed)

    # Purity auto-correlations
    ac_purity_xx = normalized_autocorr(purity_xx)
    ac_purity_zz = normalized_autocorr(purity_zz)
    ac_purity_mixed = normalized_autocorr(purity_mixed)

    # Entropy auto-correlations
    ac_entropy_xx = normalized_autocorr(entropy_xx)
    ac_entropy_zz = normalized_autocorr(entropy_zz)
    ac_entropy_mixed = normalized_autocorr(entropy_mixed)

    # Plot auto-correlation functions
    lag_values = np.arange(len(ac_coherence_xx))

    # Coherence auto-correlation
    ax1 = axes[0, 0]
    ax1.plot(lag_values, ac_coherence_xx, "b-", label="XX Coupling")
    ax1.plot(lag_values, ac_coherence_mixed, "g-", label="Mixed Coupling")
    ax1.set_xlabel("Lag (Clock Steps)")
    ax1.set_ylabel("Auto-correlation")
    ax1.set_title("Coherence Auto-correlation")
    ax1.grid(True)
    ax1.legend()

    # Purity auto-correlation
    ax2 = axes[1, 0]
    ax2.plot(lag_values, ac_purity_xx, "b-", label="XX Coupling")
    ax2.plot(lag_values, ac_purity_zz, "r-", label="ZZ Coupling")
    ax2.plot(lag_values, ac_purity_mixed, "g-", label="Mixed Coupling")
    ax2.set_xlabel("Lag (Clock Steps)")
    ax2.set_ylabel("Auto-correlation")
    ax2.set_title("Purity Auto-correlation")
    ax2.grid(True)
    ax2.legend()

    # Entropy auto-correlation
    ax3 = axes[2, 0]
    ax3.plot(lag_values, ac_entropy_xx, "b-", label="XX Coupling")
    ax3.plot(lag_values, ac_entropy_zz, "r-", label="ZZ Coupling")
    ax3.plot(lag_values, ac_entropy_mixed, "g-", label="Mixed Coupling")
    ax3.set_xlabel("Lag (Clock Steps)")
    ax3.set_ylabel("Auto-correlation")
    ax3.set_title("Entropy Auto-correlation")
    ax3.grid(True)
    ax3.legend()

    # Compute recurrence quantifiers
    # Distance matrices
    def compute_distance_matrix(signal):
        n = len(signal)
        dist_matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                dist_matrix[i, j] = abs(signal[i] - signal[j])
        return dist_matrix

    # Compute recurrence plots
    # Coherence recurrence (XX coupling)
    dist_coherence_xx = compute_distance_matrix(coherence_xx)

    # Purity recurrence (XX coupling)
    dist_purity_xx = compute_distance_matrix(purity_xx)

    # Entropy recurrence (XX coupling)
    dist_entropy_xx = compute_distance_matrix(entropy_xx)

    # Plot recurrence plots
    ax4 = axes[0, 1]
    im4 = ax4.imshow(
        dist_coherence_xx,
        origin="lower",
        cmap="viridis",
        extent=[0, D_clock, 0, D_clock],
    )
    ax4.set_xlabel("Clock Time T")
    ax4.set_ylabel("Clock Time T")
    ax4.set_title("Coherence Distance Matrix (XX Coupling)")
    plt.colorbar(im4, ax=ax4)

    ax5 = axes[1, 1]
    im5 = ax5.imshow(
        dist_purity_xx, origin="lower", cmap="viridis", extent=[0, D_clock, 0, D_clock]
    )
    ax5.set_xlabel("Clock Time T")
    ax5.set_ylabel("Clock Time T")
    ax5.set_title("Purity Distance Matrix (XX Coupling)")
    plt.colorbar(im5, ax=ax5)

    ax6 = axes[2, 1]
    im6 = ax6.imshow(
        dist_entropy_xx, origin="lower", cmap="viridis", extent=[0, D_clock, 0, D_clock]
    )
    ax6.set_xlabel("Clock Time T")
    ax6.set_ylabel("Clock Time T")
    ax6.set_title("Entropy Distance Matrix (XX Coupling)")
    plt.colorbar(im6, ax=ax6)

    plt.tight_layout()
    plt.savefig("images/recurrence_analysis.png", dpi=300)
    plt.show()

    # Calculate recurrence times
    # Find peaks in auto-correlation
    def find_recurrence_times(ac, threshold=0.5):
        peaks = []
        for i in range(2, len(ac) - 1):
            if ac[i] > threshold and ac[i] > ac[i - 1] and ac[i] > ac[i + 1]:
                peaks.append(i)
        return peaks

    # Find recurrence times
    rec_times_coherence_xx = find_recurrence_times(ac_coherence_xx)
    rec_times_purity_xx = find_recurrence_times(ac_purity_xx)
    rec_times_entropy_xx = find_recurrence_times(ac_entropy_xx)

    # Print recurrence times
    print("Recurrence analysis results:")
    print(f"Coherence recurrence times (XX coupling): {rec_times_coherence_xx}")
    print(f"Purity recurrence times (XX coupling): {rec_times_purity_xx}")
    print(f"Entropy recurrence times (XX coupling): {rec_times_entropy_xx}")

    # Return recurrence information
    return {
        "recurrence_times_coherence": rec_times_coherence_xx,
        "recurrence_times_purity": rec_times_purity_xx,
        "recurrence_times_entropy": rec_times_entropy_xx,
        "autocorr_coherence": ac_coherence_xx,
        "autocorr_purity": ac_purity_xx,
        "autocorr_entropy": ac_entropy_xx,
    }


###############################################################################
# 7. Initial state analysis
###############################################################################
def analyze_initial_states():
    """Analyze how different initial system states affect evolution."""
    # Define different initial states for comparison
    # State 1: |+⟩⊗|+⟩ (superposition - the default)
    plus = (1.0 / np.sqrt(2)) * np.array([1, 1], dtype=complex)
    psi0_plus_plus = kron(plus, plus)

    # State 2: |0⟩⊗|0⟩ (computational basis state)
    zero = np.array([1, 0], dtype=complex)
    psi0_zero_zero = kron(zero, zero)

    # State 3: |+⟩⊗|0⟩ (mixed)
    psi0_plus_zero = kron(plus, zero)

    # State 4: Bell state (|00⟩ + |11⟩)/√2 (entangled)
    psi0_bell = (1.0 / np.sqrt(2)) * np.array([1, 0, 0, 1], dtype=complex)

    # Set environment coupling parameters
    g_env1 = 0.05
    g_env2 = 0.05

    # Modified simulation function to use custom initial state
    def run_sim_with_state(init_state, coupling_type="xx"):
        """Run simulation with a given initial system state."""
        # Environment coupling
        if coupling_type == "xx":
            # σx⊗σx coupling
            H_env1 = g_env1 * four_qubit_op(sigma_x, I2, sigma_x, I2)
            H_env2 = g_env2 * four_qubit_op(I2, sigma_x, I2, sigma_x)
        else:
            # σz⊗σz coupling
            H_env1 = g_env1 * four_qubit_op(sigma_z, I2, sigma_z, I2)
            H_env2 = g_env2 * four_qubit_op(I2, sigma_z, I2, sigma_z)

        # Total Hamiltonian
        H_sys = omega * (H_1 + H_2) + H_int + H_env1 + H_env2

        # Time evolution operator
        U = expm(-1j * H_sys * dt)

        # Initial state: custom system state, environment in |0⟩⊗|0⟩
        psi0_total = kron(kron(init_state, zero), zero)
        psi0_total /= np.linalg.norm(psi0_total)

        # Clock basis
        clock_basis = np.eye(D_clock, dtype=complex)

        # Build global state
        dim_system = 16  # 4 qubits
        dim_total = D_clock * dim_system
        Psi_global = np.zeros(dim_total, dtype=complex)

        for T in range(D_clock):
            clock_state = clock_basis[T, :]
            sys_state = np.linalg.matrix_power(U, T).dot(psi0_total)
            kron_state = np.kron(clock_state, sys_state)
            Psi_global += kron_state

        Psi_global /= np.linalg.norm(Psi_global)

        # Reshape for analysis
        Psi_global_reshaped = np.reshape(Psi_global, (D_clock, dim_system))

        # Initialize result arrays
        coherence1 = []
        coherence2 = []
        purity1 = []
        purity2 = []
        entropy = []

        # Analyze each clock time
        for T in range(D_clock):
            # Get system state for this clock time
            block = Psi_global_reshaped[T, :]
            prob_T = np.vdot(block, block).real

            if prob_T > 1e-12:
                sys_state = block / np.sqrt(prob_T)
            else:
                sys_state = block

            # Full density matrix
            rho_full = np.outer(sys_state, sys_state.conjugate())

            # Get reduced density matrices and measure coherence
            rho_1 = partial_trace_to_qubit(rho_full, 0)
            rho_2 = partial_trace_to_qubit(rho_full, 1)
            coherence1.append(np.abs(rho_1[0, 1]))
            coherence2.append(np.abs(rho_2[0, 1]))

            # Calculate purity for each qubit
            purity1.append(np.real_if_close(np.trace(rho_1 @ rho_1)))
            purity2.append(np.real_if_close(np.trace(rho_2 @ rho_2)))

            # Calculate system-environment entanglement entropy
            rho_sys = get_system_density_matrix(sys_state)
            try:
                entropy_val = entanglement_entropy(rho_sys)
                entropy.append(entropy_val)
            except Exception as e:
                entropy.append(0.0)

        # Return results
        return {
            "coherence1": np.array(coherence1),
            "coherence2": np.array(coherence2),
            "purity1": np.array(purity1),
            "purity2": np.array(purity2),
            "entropy": np.array(entropy),
        }

    # Run simulations for each initial state with XX coupling
    results_plus_plus = run_sim_with_state(psi0_plus_plus, "xx")
    results_zero_zero = run_sim_with_state(psi0_zero_zero, "xx")
    results_plus_zero = run_sim_with_state(psi0_plus_zero, "xx")
    results_bell = run_sim_with_state(psi0_bell, "xx")

    # Create plots
    fig, axes = plt.subplots(3, 1, figsize=(12, 15))

    # Time values
    T_values = np.arange(D_clock)

    # Plot coherence
    ax1 = axes[0]
    ax1.plot(
        T_values, results_plus_plus["coherence1"], "b-", label="|+⟩⊗|+⟩", linewidth=1.5
    )
    ax1.plot(
        T_values, results_zero_zero["coherence1"], "r-", label="|0⟩⊗|0⟩", linewidth=1.5
    )
    ax1.plot(
        T_values, results_plus_zero["coherence1"], "g-", label="|+⟩⊗|0⟩", linewidth=1.5
    )
    ax1.plot(
        T_values, results_bell["coherence1"], "k-", label="Bell state", linewidth=1.5
    )
    ax1.set_xlabel("Clock Time T")
    ax1.set_ylabel("Coherence |ρ01|")
    ax1.set_title(
        "Qubit 1 Coherence vs Clock Time for Different Initial States (XX Coupling)"
    )
    ax1.grid(True)
    ax1.legend()

    # Plot purity
    ax2 = axes[1]
    ax2.plot(
        T_values, results_plus_plus["purity1"], "b-", label="|+⟩⊗|+⟩", linewidth=1.5
    )
    ax2.plot(
        T_values, results_zero_zero["purity1"], "r-", label="|0⟩⊗|0⟩", linewidth=1.5
    )
    ax2.plot(
        T_values, results_plus_zero["purity1"], "g-", label="|+⟩⊗|0⟩", linewidth=1.5
    )
    ax2.plot(T_values, results_bell["purity1"], "k-", label="Bell state", linewidth=1.5)
    ax2.set_xlabel("Clock Time T")
    ax2.set_ylabel("Purity Tr(ρ²)")
    ax2.set_title(
        "Qubit 1 Purity vs Clock Time for Different Initial States (XX Coupling)"
    )
    ax2.grid(True)
    ax2.legend()

    # Plot entropy
    ax3 = axes[2]
    ax3.plot(
        T_values, results_plus_plus["entropy"], "b-", label="|+⟩⊗|+⟩", linewidth=1.5
    )
    ax3.plot(
        T_values, results_zero_zero["entropy"], "r-", label="|0⟩⊗|0⟩", linewidth=1.5
    )
    ax3.plot(
        T_values, results_plus_zero["entropy"], "g-", label="|+⟩⊗|0⟩", linewidth=1.5
    )
    ax3.plot(T_values, results_bell["entropy"], "k-", label="Bell state", linewidth=1.5)
    ax3.set_xlabel("Clock Time T")
    ax3.set_ylabel("Entanglement Entropy")
    ax3.set_title(
        "System-Environment Entanglement Entropy for Different Initial States (XX Coupling)"
    )
    ax3.grid(True)
    ax3.legend()

    plt.tight_layout()
    plt.savefig("images/initial_state_comparison.png", dpi=300)
    plt.show()

    # Analyze how initial states affect recurrence dynamics
    # Compute auto-correlation for coherence of each initial state
    def normalized_autocorr(x):
        result = correlate(x - np.mean(x), x - np.mean(x), mode="full")
        return result[len(result) // 2 :] / result[len(result) // 2]

    ac_plus_plus = normalized_autocorr(results_plus_plus["coherence1"])
    ac_zero_zero = normalized_autocorr(results_zero_zero["coherence1"])
    ac_plus_zero = normalized_autocorr(results_plus_zero["coherence1"])
    ac_bell = normalized_autocorr(results_bell["coherence1"])

    # Plot autocorrelation functions
    fig, ax = plt.subplots(figsize=(10, 6))
    lag_values = np.arange(len(ac_plus_plus))

    ax.plot(lag_values, ac_plus_plus, "b-", label="|+⟩⊗|+⟩", linewidth=1.5)
    ax.plot(lag_values, ac_zero_zero, "r-", label="|0⟩⊗|0⟩", linewidth=1.5)
    ax.plot(lag_values, ac_plus_zero, "g-", label="|+⟩⊗|0⟩", linewidth=1.5)
    ax.plot(lag_values, ac_bell, "k-", label="Bell state", linewidth=1.5)
    ax.set_xlabel("Lag (Clock Steps)")
    ax.set_ylabel("Auto-correlation")
    ax.set_title(
        "Coherence Auto-correlation for Different Initial States (XX Coupling)"
    )
    ax.grid(True)
    ax.legend()

    plt.tight_layout()
    plt.savefig("images/initial_state_autocorrelation.png", dpi=300)
    plt.show()

    # Return results for potential further analysis
    return {
        "results_plus_plus": results_plus_plus,
        "results_zero_zero": results_zero_zero,
        "results_plus_zero": results_plus_zero,
        "results_bell": results_bell,
    }


###############################################################################
# 8. System size scaling analysis
###############################################################################
def analyze_environment_scaling():
    """Analyze how environment size affects decoherence and recurrence patterns."""
    # Define parameters for the fixed system (2 qubits)
    plus = (1.0 / np.sqrt(2)) * np.array([1, 1], dtype=complex)
    zero = np.array([1, 0], dtype=complex)
    psi0_sys = kron(plus, plus)  # system state |+⟩⊗|+⟩

    # Fixed system energy scale and coupling
    omega = 1.0  # system energy scale
    g = 0.5  # system-system interaction
    dt = 0.2  # time step

    # Environment coupling strength
    g_env = 0.05  # fixed for all environment qubits

    # Number of clock states (reduced for computational efficiency)
    D_clock_reduced = 50  # use smaller clock dimension for larger environments

    # Define operator function for arbitrary qubit systems
    def system_operator(op_list):
        """Create a multi-qubit operator from a list of single-qubit operators.

        Args:
            op_list: List of operators [op_1, op_2, ..., op_n]
                     where each op_i acts on the i-th qubit

        Returns:
            Tensor product of all operators
        """
        return reduce(np.kron, op_list)

    def run_simulation_with_env_size(n_env, coupling_type="xx"):
        """Run simulation with a variable number of environment qubits.

        Args:
            n_env: Number of environment qubits (1-4)
            coupling_type: Type of system-environment coupling ("xx" or "zz")

        Returns:
            Simulation results
        """
        print(f"Running simulation with {n_env} environment qubits...")

        # Total number of qubits (2 system + n_env environment)
        n_total = 2 + n_env

        # System Hamiltonian terms (fixed for all simulations)
        # First system qubit σz
        H1_ops = [sigma_z] + [I2] * (n_total - 1)
        H_1 = system_operator(H1_ops)

        # Second system qubit σz
        H2_ops = [I2, sigma_z] + [I2] * (n_total - 2)
        H_2 = system_operator(H2_ops)

        # System-system interaction
        Hint_ops = [sigma_x, sigma_x] + [I2] * (n_total - 2)
        H_int = g * system_operator(Hint_ops)

        # Base system Hamiltonian
        H_sys = omega * (H_1 + H_2) + H_int

        # Environment coupling terms
        H_env_terms = []

        for i in range(n_env):
            # System qubit 1 coupled to environment qubit i
            if coupling_type == "xx":
                # σx⊗σx coupling
                ops1 = [I2] * n_total
                ops1[0] = sigma_x  # First system qubit
                ops1[2 + i] = sigma_x  # Environment qubit i
                H_env1 = g_env * system_operator(ops1)

                # System qubit 2 coupled to environment qubit i
                ops2 = [I2] * n_total
                ops2[1] = sigma_x  # Second system qubit
                ops2[2 + i] = sigma_x  # Environment qubit i
                H_env2 = g_env * system_operator(ops2)
            else:  # "zz" coupling
                # σz⊗σz coupling
                ops1 = [I2] * n_total
                ops1[0] = sigma_z  # First system qubit
                ops1[2 + i] = sigma_z  # Environment qubit i
                H_env1 = g_env * system_operator(ops1)

                # System qubit 2 coupled to environment qubit i
                ops2 = [I2] * n_total
                ops2[1] = sigma_z  # Second system qubit
                ops2[2 + i] = sigma_z  # Environment qubit i
                H_env2 = g_env * system_operator(ops2)

            H_env_terms.append(H_env1)
            H_env_terms.append(H_env2)

        # Full Hamiltonian
        H_total = H_sys + sum(H_env_terms)

        # Time evolution operator
        U = expm(-1j * H_total * dt)

        # Initial state: system in |+⟩⊗|+⟩, environment in |0⟩⊗|0⟩⊗...
        psi0_env = zero
        for _ in range(n_env - 1):
            psi0_env = kron(psi0_env, zero)

        psi0_total = kron(psi0_sys, psi0_env)
        psi0_total /= np.linalg.norm(psi0_total)

        # Clock basis
        clock_basis = np.eye(D_clock_reduced, dtype=complex)

        # Build global state
        dim_system = 2 ** (n_total)  # 2^n qubits
        dim_total = D_clock_reduced * dim_system
        Psi_global = np.zeros(dim_total, dtype=complex)

        for T in range(D_clock_reduced):
            clock_state = clock_basis[T, :]
            sys_state = np.linalg.matrix_power(U, T).dot(psi0_total)
            kron_state = np.kron(clock_state, sys_state)
            Psi_global += kron_state

        Psi_global /= np.linalg.norm(Psi_global)

        # Reshape for analysis
        Psi_global_reshaped = np.reshape(Psi_global, (D_clock_reduced, dim_system))

        # Initialize result arrays
        coherence1 = []
        coherence2 = []
        purity1 = []
        purity2 = []
        system_purity = []

        # Analyze each clock time
        for T in range(D_clock_reduced):
            # Get system state for this clock time
            block = Psi_global_reshaped[T, :]
            prob_T = np.vdot(block, block).real

            if prob_T > 1e-12:
                sys_state = block / np.sqrt(prob_T)
            else:
                sys_state = block

            # Full density matrix
            rho_full = np.outer(sys_state, sys_state.conjugate())

            # Extract reduced density matrix for the system (first 2 qubits)
            # For large systems, we use a direct approach to avoid reshaping large arrays
            rho_system = extract_system_density_matrix(rho_full, n_total)
            system_purity.append(np.real_if_close(np.trace(rho_system @ rho_system)))

            # Extract single-qubit density matrices
            rho_1 = extract_qubit_density_matrix(rho_full, 0, n_total)
            rho_2 = extract_qubit_density_matrix(rho_full, 1, n_total)

            # Measure coherence
            coherence1.append(np.abs(rho_1[0, 1]))
            coherence2.append(np.abs(rho_2[0, 1]))

            # Calculate purity for each qubit
            purity1.append(np.real_if_close(np.trace(rho_1 @ rho_1)))
            purity2.append(np.real_if_close(np.trace(rho_2 @ rho_2)))

        # Return results
        return {
            "coherence1": np.array(coherence1),
            "coherence2": np.array(coherence2),
            "purity1": np.array(purity1),
            "purity2": np.array(purity2),
            "system_purity": np.array(system_purity),
            "n_env": n_env,
        }

    # Helper function to extract system density matrix for arbitrary size
    def extract_system_density_matrix(rho_full, n_total):
        """Extract the 4x4 density matrix for the 2 system qubits.

        Args:
            rho_full: Full density matrix (2^n x 2^n)
            n_total: Total number of qubits

        Returns:
            2-qubit system density matrix (4x4)
        """
        # Dimensions
        dim_full = 2**n_total
        dim_env = 2 ** (n_total - 2)
        dim_sys = 4

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

    # Helper function to extract single-qubit density matrix
    def extract_qubit_density_matrix(rho_full, qubit_index, n_total):
        """Extract single-qubit density matrix by tracing out others.

        Args:
            rho_full: Full density matrix (2^n x 2^n)
            qubit_index: Index of the qubit to keep (0-based)
            n_total: Total number of qubits

        Returns:
            Single-qubit density matrix (2x2)
        """
        # Only implemented for first two qubits (system qubits)
        if qubit_index > 1:
            raise ValueError("Only implemented for system qubits 0 and 1")

        # Dimensions
        dim_full = 2**n_total
        dim_rest = 2 ** (n_total - 1)

        # Initialize the reduced density matrix
        rho_qubit = np.zeros((2, 2), dtype=complex)

        # Perform partial trace by direct summation
        for i in range(2):
            for j in range(2):
                for k in range(dim_rest):
                    if qubit_index == 0:
                        # First qubit: stride is 2*dim_rest
                        row = i * dim_rest + k
                        col = j * dim_rest + k
                    else:  # qubit_index == 1
                        # Second qubit: stride is dim_rest/2
                        block_size = dim_rest // 2
                        block = k // block_size
                        offset = k % block_size
                        row = block * 2 * block_size + i * block_size + offset
                        col = block * 2 * block_size + j * block_size + offset

                    rho_qubit[i, j] += rho_full[row, col]

        return rho_qubit

    # Run simulations with different environment sizes
    # Note: Each additional qubit doubles the memory requirement
    env_sizes = [1, 2, 3, 4]  # 1 to 4 environment qubits

    results_xx = {}
    results_zz = {}

    for n_env in env_sizes:
        results_xx[n_env] = run_simulation_with_env_size(n_env, "xx")
        results_zz[n_env] = run_simulation_with_env_size(n_env, "zz")

    # Time values
    T_values = np.arange(D_clock_reduced)

    # Create plots for XX coupling
    fig, axes = plt.subplots(2, 1, figsize=(12, 10))

    # Plot coherence vs environment size
    ax1 = axes[0]
    colors = ["b", "g", "r", "k"]

    for i, n_env in enumerate(env_sizes):
        ax1.plot(
            T_values,
            results_xx[n_env]["coherence1"],
            color=colors[i],
            linestyle="-",
            label=f"{n_env} Env. Qubits",
            linewidth=1.5,
        )

    ax1.set_xlabel("Clock Time T")
    ax1.set_ylabel("Coherence |ρ01|")
    ax1.set_title("Qubit 1 Coherence vs Environment Size (XX Coupling)")
    ax1.grid(True)
    ax1.legend()

    # Plot system purity vs environment size
    ax2 = axes[1]

    for i, n_env in enumerate(env_sizes):
        ax2.plot(
            T_values,
            results_xx[n_env]["system_purity"],
            color=colors[i],
            linestyle="-",
            label=f"{n_env} Env. Qubits",
            linewidth=1.5,
        )

    ax2.set_xlabel("Clock Time T")
    ax2.set_ylabel("System Purity Tr(ρ²)")
    ax2.set_title("System Purity vs Environment Size (XX Coupling)")
    ax2.grid(True)
    ax2.legend()

    plt.tight_layout()
    plt.savefig("images/environment_scaling_xx.png", dpi=300)
    plt.show()

    # Create plots for ZZ coupling
    fig, axes = plt.subplots(2, 1, figsize=(12, 10))

    # Plot coherence vs environment size
    ax1 = axes[0]

    for i, n_env in enumerate(env_sizes):
        ax1.plot(
            T_values,
            results_zz[n_env]["coherence1"],
            color=colors[i],
            linestyle="-",
            label=f"{n_env} Env. Qubits",
            linewidth=1.5,
        )

    ax1.set_xlabel("Clock Time T")
    ax1.set_ylabel("Coherence |ρ01|")
    ax1.set_title("Qubit 1 Coherence vs Environment Size (ZZ Coupling)")
    ax1.grid(True)
    ax1.legend()

    # Plot system purity vs environment size
    ax2 = axes[1]

    for i, n_env in enumerate(env_sizes):
        ax2.plot(
            T_values,
            results_zz[n_env]["system_purity"],
            color=colors[i],
            linestyle="-",
            label=f"{n_env} Env. Qubits",
            linewidth=1.5,
        )

    ax2.set_xlabel("Clock Time T")
    ax2.set_ylabel("System Purity Tr(ρ²)")
    ax2.set_title("System Purity vs Environment Size (ZZ Coupling)")
    ax2.grid(True)
    ax2.legend()

    plt.tight_layout()
    plt.savefig("images/environment_scaling_zz.png", dpi=300)
    plt.show()

    # Analyze recurrence properties vs environment size
    fig, axes = plt.subplots(2, 1, figsize=(12, 10))

    # Compute autocorrelation functions
    def normalized_autocorr(x):
        result = correlate(x - np.mean(x), x - np.mean(x), mode="full")
        return result[len(result) // 2 :] / result[len(result) // 2]

    # Plot autocorrelation for XX coupling
    ax1 = axes[0]

    for i, n_env in enumerate(env_sizes):
        ac = normalized_autocorr(results_xx[n_env]["coherence1"])
        ax1.plot(
            np.arange(len(ac)),
            ac,
            color=colors[i],
            linestyle="-",
            label=f"{n_env} Env. Qubits",
            linewidth=1.5,
        )

    ax1.set_xlabel("Lag (Clock Steps)")
    ax1.set_ylabel("Autocorrelation")
    ax1.set_title("Coherence Autocorrelation vs Environment Size (XX Coupling)")
    ax1.grid(True)
    ax1.legend()

    # Plot autocorrelation for ZZ coupling
    ax2 = axes[1]

    for i, n_env in enumerate(env_sizes):
        if np.any(
            results_zz[n_env]["coherence1"] > 1e-10
        ):  # Only if meaningful coherence exists
            ac = normalized_autocorr(results_zz[n_env]["coherence1"])
            ax2.plot(
                np.arange(len(ac)),
                ac,
                color=colors[i],
                linestyle="-",
                label=f"{n_env} Env. Qubits",
                linewidth=1.5,
            )

    ax2.set_xlabel("Lag (Clock Steps)")
    ax2.set_ylabel("Autocorrelation")
    ax2.set_title("Coherence Autocorrelation vs Environment Size (ZZ Coupling)")
    ax2.grid(True)
    ax2.legend()

    plt.tight_layout()
    plt.savefig("images/environment_scaling_autocorr.png", dpi=300)
    plt.show()

    # Return results for potential further analysis
    return {"results_xx": results_xx, "results_zz": results_zz}


###############################################################################
# 9. Main execution
###############################################################################
if __name__ == "__main__":
    print("Emergent Time with Environment Coupling Analysis")
    print("1. Coupling Strength Sweep")
    print("2. Compare Coupling Types")
    print("3. Advanced Analysis (Frequency and Recurrence)")
    print("4. Initial State Analysis")
    print("5. Environment Size Scaling Analysis")

    choice = input("Enter your choice (1-5): ")

    if choice == "1":
        print("Running coupling strength sweep (this may take a while)...")
        results = coupling_strength_sweep()
        print("Sweep complete. Results saved to 'images/coupling_sweep_results.png'")
    elif choice == "2":
        print("Comparing different coupling types...")
        compare_coupling_types()
        print(
            "Comparison complete. Results saved to 'images/coupling_type_comparison.png'"
        )
    elif choice == "3":
        print("Running advanced analysis...")
        # Run simulations for three coupling types
        g_env1 = 0.05
        g_env2 = 0.05
        results_xx = run_simulation(g_env1, g_env2, coupling_type="xx")
        results_zz = run_simulation(g_env1, g_env2, coupling_type="zz")
        results_mixed = run_simulation(g_env1, g_env2, coupling_type="mixed")

        # Run frequency analysis
        print("Analyzing frequency components...")
        analyze_frequency_components(results_xx, results_zz, results_mixed)

        # Run recurrence analysis
        print("Analyzing recurrence dynamics...")
        recurrence_info = analyze_recurrence_dynamics(
            results_xx, results_zz, results_mixed
        )

        print(
            "Advanced analysis complete. Results saved to 'images/frequency_analysis.png' and 'images/recurrence_analysis.png'"
        )
    elif choice == "4":
        print("Analyzing different initial states...")
        analyze_initial_states()
        print(
            "Initial state analysis complete. Results saved to 'images/initial_state_comparison.png' and 'images/initial_state_autocorrelation.png'"
        )
    elif choice == "5":
        print(
            "Running environment size scaling analysis (this may take a long time)..."
        )
        print("WARNING: This analysis requires significant computational resources.")
        print("The simulation time increases exponentially with environment size.")
        confirm = input("Are you sure you want to continue? (y/n): ")
        if confirm.lower() == "y":
            scaling_results = analyze_environment_scaling()
            print(
                "Environment scaling analysis complete. Results saved to 'images/environment_scaling_*.png'"
            )
        else:
            print("Analysis cancelled.")
    else:
        print("Invalid choice. Please run again and select options 1-5.")
