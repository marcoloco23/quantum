# Quantum Emergence

A Python package for simulating and analyzing quantum systems that demonstrate emergent time behavior according to the Page-Wootters mechanism.

## Overview

This package provides tools for exploring quantum dynamics in systems where time emerges from entanglement between a "clock" subsystem and the rest of the quantum state. The framework is based on the Page-Wootters mechanism, where a global timeless state gives rise to perceived dynamics through quantum correlations.

## Installation

```bash
# Clone the repository
git clone https://github.com/marcoloco23/quantum.git
cd emergence

# Install dependencies
pip install -r requirements.txt

# Install the package
pip install -e .
```

## Package Structure

The package consists of several modules:

1. **quantum_core.py**: Core quantum mechanics utilities and constants
   - Quantum operators (Pauli matrices)
   - Common quantum states
   - Density matrix utilities
   - Measurement functions
   - Entanglement and coherence metrics

2. **quantum_simulator.py**: Simulation classes and functions
   - `SystemConfiguration` class for simulation parameters
   - `QuantumSimulator` class for running simulations
   - Functions for custom simulations with different initial states and environment sizes

3. **analysis.py**: Analysis utilities
   - Signal processing functions (FFT, autocorrelation)
   - Recurrence dynamics analysis
   - Parameter sweep analysis
   - Statistical functions for comparing metrics across simulations

4. **visualization.py**: Plotting functions
   - Time series plots
   - Heatmaps and 3D surface plots
   - Specialized plots for simulation results, coupling sweeps, and frequency analysis

5. **main.py**: Command-line interface for running analyses
   - Subcommands for different types of analyses
   - Options for customizing simulations and visualizations

## Usage Examples

### Basic Simulation

```python
from emergence import SystemConfiguration, QuantumSimulator

# Configure the simulation
config = SystemConfiguration(
    num_clock_states=20,
    omega=1.0,
    system_coupling=0.5,
    env_coupling1=0.2,
    env_coupling2=0.2,
    time_step=0.1,
    coupling_type="xx",  # Case-insensitive: "xx", "zz", or "mixed"
    num_env_qubits=2
)

# Create simulator and run
simulator = QuantumSimulator(config)
results = simulator.run_simulation()

# Access results
print(f"Final coherence: {results['coherence1'][-1]}")
print(f"Final purity: {results['purity1'][-1]}")
```

### Running Analyses from Command Line

```bash
# Sweep coupling strengths (coupling type is case-insensitive)
python -m emergence.main sweep --coupling-type xx --min-coupling 0.0 --max-coupling 1.0 --num-points 10

# Compare different coupling types
python -m emergence.main compare --coupling-strength 0.3

# Run advanced analysis (frequency and recurrence)
python -m emergence.main advanced --coupling-strength 0.05

# Analyze different initial states
python -m emergence.main initial --coupling-strength 0.05

# Analyze environment scaling effects
python -m emergence.main env-scaling --max-size 4 --coupling-type xx
```

## Key Concepts

### Page-Wootters Mechanism

The Page-Wootters mechanism proposes that time emerges from entanglement between a "clock" degree of freedom and the rest of a quantum system. In this framework:

- The global state is timeless and satisfies the Wheeler-DeWitt equation (Hψ = 0)
- Perceived dynamics arise by conditioning on different "clock" states
- The framework naturally explains how time emerges in quantum systems

### Quantum Decoherence

Decoherence is a quantum phenomenon where a system loses its quantum coherence due to interaction with its environment:

- Coherence measures the quantum superposition in a system
- Purity quantifies how mixed (or pure) a quantum state is
- Entanglement entropy measures the entanglement between system and environment

This package allows exploration of how these metrics evolve under different coupling conditions.

## Related Components

This package is part of a larger quantum simulation framework that includes:
- **Ising model simulations**: Classical and quantum Ising models for studying phase transitions
- **SYK model simulations**: Sachdev-Ye-Kitaev model for exploring quantum chaos
- **Web interface**: Interactive visualization and parameter exploration

See the main repository README for more information on these components.

## License

MIT License

