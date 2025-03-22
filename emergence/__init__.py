"""
Emergence: Quantum Dynamics Simulator for Emergent Time Analysis

This package provides tools for simulating and analyzing quantum systems that
demonstrate emergent time behavior according to the Page-Wootters mechanism.
"""

__version__ = "1.0.0"
__author__ = "Quantum Matter Research Group"

# Import core modules to make them available at package level
from .quantum_core import (
    sigma_x,
    sigma_y,
    sigma_z,
    I2,
    zero,
    one,
    plus,
    minus,
    bell_plus,
    bell_minus,
    entanglement_entropy,
    coherence,
    purity,
)

from .quantum_simulator import (
    SystemConfiguration,
    QuantumSimulator,
    run_simulation_with_custom_state,
    run_simulation_with_env_size,
)

# Import functions for direct usage
from .analysis import (
    normalized_autocorr,
    compute_fft,
    compute_distance_matrix,
    find_recurrence_times,
    analyze_frequency_components,
    analyze_recurrence_dynamics,
)

from .visualization import (
    plot_simulation_results,
    plot_coupling_strength_sweep,
    plot_coupling_type_comparison,
    plot_frequency_analysis,
)
