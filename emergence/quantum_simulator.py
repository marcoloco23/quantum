import numpy as np
from numpy import kron

from .quantum_core import (
    sigma_x,
    sigma_z,
    I2,
    zero,
    plus,
    four_qubit_op,
    multi_qubit_op,
    density_matrix,
    partial_trace_to_qubit,
    get_system_density_matrix,
    entanglement_entropy,
    purity,
    coherence,
    expectation_value,
    time_evolution_operator,
)

"""
Quantum simulator module for emergent time dynamics.
"""


class SystemConfiguration:
    """Configuration parameters for quantum simulations."""

    def __init__(
        self,
        num_clock_states: int = 100,
        omega: float = 1.0,
        system_coupling: float = 0.5,
        env_coupling1: float = 0.05,
        env_coupling2: float = 0.05,
        time_step: float = 0.2,
        coupling_type: str = "xx",
        num_env_qubits: int = 2,
    ):
        """Initialize system configuration.

        Args:
            num_clock_states: Number of clock states
            omega: System energy scale
            system_coupling: Interaction strength between system qubits
            env_coupling1: Coupling strength to first environment qubit
            env_coupling2: Coupling strength to second environment qubit
            time_step: Time step for evolution
            coupling_type: Type of system-environment coupling ("xx", "zz", or "mixed")
            num_env_qubits: Number of environment qubits
        """
        self.D_clock = num_clock_states
        self.omega = omega
        self.g = system_coupling
        self.g_env1 = env_coupling1
        self.g_env2 = env_coupling2
        self.dt = time_step
        self.coupling_type = coupling_type
        self.num_env_qubits = num_env_qubits

        # Total number of qubits
        self.num_system_qubits = 2
        self.num_total_qubits = self.num_system_qubits + self.num_env_qubits


class QuantumSimulator:
    """Simulator for quantum dynamics with emergent time."""

    def __init__(self, config: SystemConfiguration):
        """Initialize the simulator with configuration.

        Args:
            config: System configuration parameters
        """
        self.config = config

        # Initialize clock basis
        self.clock_basis = np.eye(config.D_clock, dtype=complex)

        # Create default system Hamiltonian
        self._build_hamiltonian()

        # Create default initial state
        self._prepare_initial_state()

    def _build_hamiltonian(self):
        """Build the system Hamiltonian based on configuration."""
        # Default 4-qubit system
        if self.config.num_total_qubits == 4:
            # System terms (first two qubits)
            self.H_1 = four_qubit_op(sigma_z, I2, I2, I2)  # σz on first qubit
            self.H_2 = four_qubit_op(I2, sigma_z, I2, I2)  # σz on second qubit
            self.H_int = self.config.g * four_qubit_op(
                sigma_x, sigma_x, I2, I2
            )  # σxσx interaction

            # Environment coupling
            if self.config.g_env1 == 0.0 and self.config.g_env2 == 0.0:
                # No environment coupling
                self.H_total = self.config.omega * (self.H_1 + self.H_2) + self.H_int
            else:
                # Environment coupling
                if self.config.coupling_type == "xx":
                    # σx⊗σx coupling
                    self.H_env1 = self.config.g_env1 * four_qubit_op(
                        sigma_x, I2, sigma_x, I2
                    )
                    self.H_env2 = self.config.g_env2 * four_qubit_op(
                        I2, sigma_x, I2, sigma_x
                    )
                elif self.config.coupling_type == "zz":
                    # σz⊗σz coupling
                    self.H_env1 = self.config.g_env1 * four_qubit_op(
                        sigma_z, I2, sigma_z, I2
                    )
                    self.H_env2 = self.config.g_env2 * four_qubit_op(
                        I2, sigma_z, I2, sigma_z
                    )
                elif self.config.coupling_type == "mixed":
                    # Mixed σx⊗σx and σz⊗σz coupling
                    self.H_env1 = self.config.g_env1 * (
                        0.5 * four_qubit_op(sigma_x, I2, sigma_x, I2)
                        + 0.5 * four_qubit_op(sigma_z, I2, sigma_z, I2)
                    )
                    self.H_env2 = self.config.g_env2 * (
                        0.5 * four_qubit_op(I2, sigma_x, I2, sigma_x)
                        + 0.5 * four_qubit_op(I2, sigma_z, I2, sigma_z)
                    )
                else:
                    raise ValueError(
                        f"Unknown coupling type: {self.config.coupling_type}"
                    )

                self.H_total = (
                    self.config.omega * (self.H_1 + self.H_2)
                    + self.H_int
                    + self.H_env1
                    + self.H_env2
                )
        else:
            # Generic multi-qubit system implementation
            n_total = self.config.num_total_qubits

            # System Hamiltonian terms
            # First system qubit σz
            H1_ops = [sigma_z] + [I2] * (n_total - 1)
            self.H_1 = multi_qubit_op(H1_ops)

            # Second system qubit σz
            H2_ops = [I2, sigma_z] + [I2] * (n_total - 2)
            self.H_2 = multi_qubit_op(H2_ops)

            # System-system interaction σxσx
            Hint_ops = [sigma_x, sigma_x] + [I2] * (n_total - 2)
            self.H_int = self.config.g * multi_qubit_op(Hint_ops)

            # Base system Hamiltonian
            self.H_total = self.config.omega * (self.H_1 + self.H_2) + self.H_int

            # Environment coupling
            if not (self.config.g_env1 == 0.0 and self.config.g_env2 == 0.0):
                H_env_terms = []

                # Build coupling terms
                for i in range(self.config.num_env_qubits):
                    # Environment index starts at 2 (after system qubits)
                    env_idx = 2 + i

                    # Coupling strength - distribute among environment qubits
                    # Use g_env1 for first half, g_env2 for second half
                    if i < self.config.num_env_qubits // 2:
                        g_env = self.config.g_env1 / max(
                            1, self.config.num_env_qubits // 2
                        )
                    else:
                        g_env = self.config.g_env2 / max(
                            1,
                            self.config.num_env_qubits
                            - self.config.num_env_qubits // 2,
                        )

                    if (
                        self.config.coupling_type == "xx"
                        or self.config.coupling_type == "mixed"
                    ):
                        # System qubit 1 coupled to environment qubit i (σx⊗σx)
                        ops1 = [I2] * n_total
                        ops1[0] = sigma_x  # First system qubit
                        ops1[env_idx] = sigma_x  # Environment qubit i
                        H_env1 = g_env * multi_qubit_op(ops1)
                        H_env_terms.append(H_env1)

                        # System qubit 2 coupled to environment qubit i (σx⊗σx)
                        ops2 = [I2] * n_total
                        ops2[1] = sigma_x  # Second system qubit
                        ops2[env_idx] = sigma_x  # Environment qubit i
                        H_env2 = g_env * multi_qubit_op(ops2)
                        H_env_terms.append(H_env2)

                    if (
                        self.config.coupling_type == "zz"
                        or self.config.coupling_type == "mixed"
                    ):
                        # System qubit 1 coupled to environment qubit i (σz⊗σz)
                        ops1 = [I2] * n_total
                        ops1[0] = sigma_z  # First system qubit
                        ops1[env_idx] = sigma_z  # Environment qubit i
                        H_env1 = g_env * multi_qubit_op(ops1)
                        H_env_terms.append(H_env1)

                        # System qubit 2 coupled to environment qubit i (σz⊗σz)
                        ops2 = [I2] * n_total
                        ops2[1] = sigma_z  # Second system qubit
                        ops2[env_idx] = sigma_z  # Environment qubit i
                        H_env2 = g_env * multi_qubit_op(ops2)
                        H_env_terms.append(H_env2)

                # Add all environment coupling terms to total Hamiltonian
                self.H_total += sum(H_env_terms)

        # Generate time evolution operator
        self.U = time_evolution_operator(self.H_total, self.config.dt)

    def _prepare_initial_state(self, system_state=None):
        """Prepare the initial state of the system.

        Args:
            system_state: Initial state of the system qubits (default: |+>⊗|+>)
        """
        # Default system state: |+>⊗|+>
        if system_state is None:
            self.psi0_sys = kron(plus, plus)
        else:
            self.psi0_sys = system_state

        # Environment state: |0>⊗...⊗|0>
        psi0_env = zero
        for _ in range(1, self.config.num_env_qubits):
            psi0_env = kron(psi0_env, zero)

        # Total initial state
        self.psi0_total = kron(self.psi0_sys, psi0_env)
        self.psi0_total /= np.linalg.norm(self.psi0_total)

    def set_initial_state(self, system_state):
        """Set a custom initial state for the system qubits.

        Args:
            system_state: Initial state of the system qubits
        """
        self._prepare_initial_state(system_state)

    def run_simulation(self):
        """Run the quantum simulation.

        Returns:
            Dictionary containing simulation results
        """
        # System dimensions
        dim_system = 2**self.config.num_total_qubits
        dim_total = self.config.D_clock * dim_system

        # Build global state
        Psi_global = np.zeros(dim_total, dtype=complex)

        for T in range(self.config.D_clock):
            clock_state = self.clock_basis[T, :]
            sys_state = np.linalg.matrix_power(self.U, T).dot(self.psi0_total)
            kron_state = np.kron(clock_state, sys_state)
            Psi_global += kron_state

        Psi_global /= np.linalg.norm(Psi_global)

        # Reshape for analysis
        Psi_global_reshaped = np.reshape(Psi_global, (self.config.D_clock, dim_system))

        # Initialize result arrays
        clock_probs = []
        exp_values_Z1 = []
        exp_values_Z2 = []
        coherence1 = []
        coherence2 = []
        purity1 = []
        purity2 = []
        system_purity = []
        entropy = []

        # Observables
        Z1_ops = [sigma_z] + [I2] * (self.config.num_total_qubits - 1)
        Z2_ops = [I2, sigma_z] + [I2] * (self.config.num_total_qubits - 2)

        Z1 = multi_qubit_op(Z1_ops)
        Z2 = multi_qubit_op(Z2_ops)

        # Analyze each clock time
        for T in range(self.config.D_clock):
            # Get system state for this clock time
            block = Psi_global_reshaped[T, :]
            prob_T = np.vdot(block, block).real
            clock_probs.append(prob_T)

            if prob_T > 1e-12:
                sys_state = block / np.sqrt(prob_T)
            else:
                sys_state = block

            # Full density matrix
            rho_full = density_matrix(sys_state)

            # Measure observables
            exp_values_Z1.append(expectation_value(rho_full, Z1))
            exp_values_Z2.append(expectation_value(rho_full, Z2))

            # Get reduced density matrices for individual qubits
            try:
                rho_1 = partial_trace_to_qubit(
                    rho_full, 0, self.config.num_total_qubits
                )
                rho_2 = partial_trace_to_qubit(
                    rho_full, 1, self.config.num_total_qubits
                )

                # Calculate coherence
                coherence1.append(coherence(rho_1))
                coherence2.append(coherence(rho_2))

                # Calculate purity for each qubit
                purity1.append(purity(rho_1))
                purity2.append(purity(rho_2))
            except ValueError as e:
                print(f"Warning at T={T}: {e}")
                # Use fallback values if trace operation fails
                coherence1.append(0.0)
                coherence2.append(0.0)
                purity1.append(1.0)  # Pure state as fallback
                purity2.append(1.0)

            # Get system density matrix (first 2 qubits)
            rho_sys = get_system_density_matrix(
                sys_state,
                n_qubits=self.config.num_total_qubits,
                n_system_qubits=self.config.num_system_qubits,
            )

            # Calculate system purity
            system_purity.append(purity(rho_sys))

            # Calculate system-environment entanglement entropy
            try:
                entropy_val = entanglement_entropy(rho_sys)
                entropy.append(entropy_val)
            except Exception as e:
                print(f"Warning: Error calculating entropy at T={T}: {e}")
                entropy.append(0.0)

        # Return all results
        return {
            "global_state": Psi_global,
            "global_reshaped": Psi_global_reshaped,
            "clock_probs": np.array(clock_probs),
            "exp_values_Z1": np.array(exp_values_Z1),
            "exp_values_Z2": np.array(exp_values_Z2),
            "coherence1": np.array(coherence1),
            "coherence2": np.array(coherence2),
            "purity1": np.array(purity1),
            "purity2": np.array(purity2),
            "system_purity": np.array(system_purity),
            "entropy": np.array(entropy),
        }


def run_simulation_with_custom_state(
    system_state,
    env_coupling1=0.05,
    env_coupling2=0.05,
    coupling_type="xx",
    num_clock_states=100,
):
    """Run simulation with a given initial system state.

    Args:
        system_state: Initial state of the system qubits
        env_coupling1: Coupling strength to first environment qubit
        env_coupling2: Coupling strength to second environment qubit
        coupling_type: Type of system-environment coupling
        num_clock_states: Number of clock states

    Returns:
        Dictionary containing simulation results
    """
    # Create configuration
    config = SystemConfiguration(
        num_clock_states=num_clock_states,
        env_coupling1=env_coupling1,
        env_coupling2=env_coupling2,
        coupling_type=coupling_type,
    )

    # Create simulator
    simulator = QuantumSimulator(config)

    # Set initial state and run
    simulator.set_initial_state(system_state)
    results = simulator.run_simulation()

    return results


def run_simulation_with_env_size(
    n_env, coupling_type="xx", g_env=0.05, num_clock_states=50
):
    """Run simulation with a variable number of environment qubits.

    Args:
        n_env: Number of environment qubits
        coupling_type: Type of system-environment coupling
        g_env: Environment coupling strength (distributed across qubits)
        num_clock_states: Number of clock states (reduced for computational efficiency)

    Returns:
        Dictionary containing simulation results
    """
    # Split coupling strength evenly
    g_env1 = g_env
    g_env2 = g_env

    # Create configuration with specified number of environment qubits
    config = SystemConfiguration(
        num_clock_states=num_clock_states,
        env_coupling1=g_env1,
        env_coupling2=g_env2,
        coupling_type=coupling_type,
        num_env_qubits=n_env,
    )

    # Create simulator and run
    simulator = QuantumSimulator(config)
    results = simulator.run_simulation()

    # Add environment size to results
    results["n_env"] = n_env

    return results
