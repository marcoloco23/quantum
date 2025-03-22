# Sachdev-Ye-Kitaev (SYK) Model Simulations

A module for simulating and analyzing the Sachdev-Ye-Kitaev model, which exhibits properties of quantum chaos and holographic duality.

## Overview

The SYK model is a quantum many-body system of randomly interacting Majorana fermions that has received significant attention in quantum gravity research due to its connections to black holes via the AdS/CFT correspondence. This module provides tools to simulate the SYK model and calculate various properties related to quantum chaos and thermodynamics.

## Features

- Random matrix implementation of the SYK model
- Eigenvalue spectrum calculation and analysis
- Entropy and specific heat calculations
- Out-of-time-order correlators (OTOCs) for scrambling analysis
- Averages over disorder realizations for robust statistics
- Temperature-dependent observables

## Usage

### Direct Module Import

```python
from syk.syk_simulation import (
    calculate_entropy,
    calculate_otoc,
    calculate_specific_heat,
    run_multiple_syk_simulations
)

# Run N=4 SYK model with 10 disorder realizations
N = 4  # Number of Majorana fermions (must be even)
num_samples = 10
eigenvalues_list, eigenvectors_list = run_multiple_syk_simulations(N, num_samples)

# Calculate entropy at different temperatures
T_values = np.linspace(0.1, 10, 100)
entropies = [calculate_entropy(eigenvalues, T_values) for eigenvalues in eigenvalues_list]
avg_entropy = np.mean(entropies, axis=0)

# Calculate OTOC to measure quantum chaos
times, otoc_values = calculate_otoc(N, J=1.0, num_samples=10, time_steps=20)
```

### Using the Web Interface

The module can also be accessed through the Streamlit web interface:

```bash
streamlit run app.py
```

Then select "SYK Model" from the sidebar and configure the simulation parameters.

## Theoretical Background

The SYK model Hamiltonian is given by:

H = i^(q/2) ∑_{i,j,k,l} J_{ijkl} ψ_i ψ_j ψ_k ψ_l

Where:
- ψ_i are Majorana fermion operators
- J_{ijkl} are random couplings drawn from a Gaussian distribution
- The model exhibits:
  - Maximal chaos (Lyapunov exponent = 2π/β)
  - Approximate conformal symmetry in the IR
  - Connections to black hole physics through holography

## Implementation Details

This implementation uses numerical diagonalization of the Hamiltonian matrix. Key techniques include:

1. Construction of the SYK Hamiltonian using random couplings
2. Diagonalization to obtain eigenvalues and eigenvectors
3. Calculation of thermal observables using the canonical ensemble
4. Computation of OTOCs to quantify quantum chaos
5. Averaging over disorder realizations for statistical robustness

## Integration with Other Modules

This module is part of a larger quantum simulation framework that includes:
- **Emergent time simulations**: Quantum systems with Page-Wootters mechanism
- **Ising model simulations**: Classical and quantum Ising models
- **Web interface**: Interactive visualization and parameter exploration

See the main repository README for more information on these components.

## References

1. Sachdev, S., & Ye, J. (1993). *Gapless spin-fluid ground state in a random quantum Heisenberg magnet*. Physical Review Letters, 70(21), 3339.
2. Kitaev, A. (2015). *A simple model of quantum holography*. KITP Program: Entanglement in Strongly-Correlated Quantum Matter.
3. Maldacena, J., & Stanford, D. (2016). *Remarks on the Sachdev-Ye-Kitaev model*. Physical Review D, 94(10), 106002.

## License

MIT License 