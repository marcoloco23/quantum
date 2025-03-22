# Ising Model Simulations

A module for simulating and analyzing the 3D Ising model using Monte Carlo methods.

## Overview

The Ising model is one of the most fundamental models in statistical physics, used to describe phase transitions and critical phenomena. This module implements efficient Monte Carlo simulations for the 3D Ising model with configurable parameters.

## Features

- 3D lattice simulations with periodic boundary conditions
- Metropolis-Hastings Monte Carlo algorithm
- Configurable interaction strength and external field
- Calculation of key thermodynamic observables:
  - Magnetization
  - Energy
  - Magnetic susceptibility
  - Specific heat capacity
- Temperature sweep for phase transition analysis

## Usage

### Direct Module Import

```python
from ising.ising_simulation import simulate

# Parameters:
# N: Lattice size
# J: Interaction strength
# h: External magnetic field
# T: Temperature
# steps: Number of Monte Carlo steps
# equilibration_steps: Number of equilibration steps
result = simulate((N, J, h, T, steps, equilibration_steps))

# Unpack results
magnetization, energy, susceptibility, heat_capacity = result

print(f"Magnetization: {magnetization}")
print(f"Energy: {energy}")
print(f"Susceptibility: {susceptibility}")
print(f"Heat Capacity: {heat_capacity}")
```

### Using the Web Interface

The module can also be accessed through the Streamlit web interface:

```bash
streamlit run app.py
```

Then select "3D Ising Model" from the sidebar and configure simulation parameters.

## Implementation Details

The implementation uses a 3D lattice of spins (±1), with the Hamiltonian:

H = -J ∑<i,j> s_i s_j - h ∑_i s_i

Where:
- J is the coupling strength
- h is the external magnetic field
- <i,j> denotes nearest-neighbor pairs

The Monte Carlo algorithm proceeds as follows:
1. Initialize a random spin configuration
2. Perform equilibration steps to reach thermal equilibrium
3. Perform measurement steps, proposing spin flips according to the Metropolis criterion
4. Calculate observables from the equilibrated configurations

## Integration with Other Modules

This module is part of a larger quantum simulation framework that includes:
- **Emergent time simulations**: Quantum systems with Page-Wootters mechanism
- **SYK model simulations**: Models of quantum chaos
- **Web interface**: Interactive visualization and parameter exploration

See the main repository README for more information on these components.

## License

MIT License 