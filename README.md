# Quantum Matter Simulation

A comprehensive Python package for simulating and analyzing quantum systems, focusing on emergent time dynamics, Ising models, and SYK models.

## Project Components

This repository contains several interconnected quantum physics simulation modules:

1. **Emergence**: Simulates emergent time dynamics in quantum systems with environment coupling based on the Page-Wootters mechanism
2. **Ising**: Implements 3D Ising model simulations for studying phase transitions
3. **SYK**: Provides Sachdev-Ye-Kitaev model simulations for exploring quantum chaos
4. **Web Application**: Streamlit-based interface for running and visualizing simulations

## Features

- **Quantum System Simulation**:
  - Configurable environment coupling
  - Multiple coupling types (XX, ZZ, Mixed)
  - Parameter sweep analysis
  - Initial state analysis
  - Environment size scaling

- **Advanced Analysis Tools**:
  - Frequency analysis
  - Recurrence dynamics
  - Eigenvalue analysis
  - Entropy calculations
  - Out-of-time-order correlators (OTOCs)

- **Interactive Visualizations**:
  - Command-line simulation outputs
  - Web interface for interactive analysis
  - Publication-quality plots

## Installation

1. Clone the repository:
```bash
git clone https://github.com/marcoloco23/quantum.git
cd quantum
```

2. Install dependencies:
```bash
pip install -r requirements.txt
pip install -e .
```

## Usage

### Command-Line Simulations

#### Emergent Time Dynamics
```bash
# Sweep coupling strengths
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

### Web Application

Launch the interactive web application:
```bash
streamlit run app.py
```

This provides a user-friendly interface for:
- Running 3D Ising model simulations
- Exploring SYK model properties
- Visualizing results dynamically

## Project Structure

- `emergence/`: Page-Wootters mechanism and emergent time simulations
  - `main.py`: Command-line interface
  - `quantum_core.py`: Core quantum mechanics functions
  - `quantum_simulator.py`: Quantum system simulator
  - `analysis.py`: Analysis tools
  - `visualization.py`: Plotting functions

- `ising/`: Ising model implementation
  - `ising_simulation.py`: Monte Carlo simulations of the 3D Ising model

- `syk/`: Sachdev-Ye-Kitaev model implementation
  - `syk_simulation.py`: Random matrix theory approach to quantum chaos

- `notebooks/`: Jupyter notebooks for research and exploration
  - `emergent_time.ipynb`: Explorations of emergent time dynamics
  - `EntropyInequality.ipynb`: Investigations of entropy inequalities
  - `distance_matrices.ipynb`: Analysis of quantum state distances

- `app.py`: Streamlit web application for interactive simulations

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details. 