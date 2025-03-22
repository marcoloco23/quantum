# Quantum Matter Simulation

A Python package for simulating and analyzing emergent time dynamics in quantum systems with environment coupling.

## Features

- Quantum system simulation with configurable environment coupling
- Multiple coupling types (XX, ZZ, Mixed)
- Parameter sweep analysis
- Advanced analysis tools:
  - Frequency analysis
  - Recurrence dynamics
  - Initial state analysis
  - Environment size scaling

## Installation

1. Clone the repository:
```bash
git clone https://github.com/yourusername/quantum-matter.git
cd quantum-matter
```

2. Install dependencies:
```bash
pip install -e .
```

## Usage

The package provides several analysis modes that can be run through the command line:

### Coupling Strength Sweep
```bash
python -m emergence.main sweep --coupling-type xx --min-coupling 0.0 --max-coupling 1.0 --num-points 10
```

### Compare Coupling Types
```bash
python -m emergence.main compare --coupling-strength 0.3
```

### Advanced Analysis
```bash
python -m emergence.main advanced --coupling-strength 0.05
```

### Initial State Analysis
```bash
python -m emergence.main initial --coupling-strength 0.05
```

### Environment Size Scaling
```bash
python -m emergence.main env-scaling --max-size 4 --coupling-type xx
```

## Project Structure

- `emergence/`
  - `main.py`: Main script with command-line interface
  - `quantum_core.py`: Core quantum mechanics functions
  - `quantum_simulator.py`: Quantum system simulator
  - `analysis.py`: Analysis tools
  - `visualization.py`: Plotting functions

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details. 