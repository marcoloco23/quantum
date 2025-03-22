# Research Notebooks

This directory contains Jupyter notebooks for exploring and visualizing various aspects of quantum systems. These notebooks provide interactive demonstrations, theoretical investigations, and visualization tools for research purposes.

## Notebook Descriptions

### 1. `emergent_time.ipynb`

Explores the concept of emergent time in quantum systems using the Page-Wootters mechanism. This notebook demonstrates:
- Construction of timeless global states
- Conditioning on clock states to extract dynamics
- Visualization of coherence and purity evolution
- Effect of environment coupling on quantum dynamics

### 2. `EntropyInequality.ipynb`

Investigates entropy inequalities in quantum systems, focusing on:
- Strong subadditivity of entropy
- Entanglement entropy calculations
- Mutual information between subsystems
- Entropy bounds in different coupling regimes

### 3. `distance_matrices.ipynb`

Analyzes quantum state distance metrics to quantify:
- Recurrence phenomena in quantum dynamics
- Trace distance between states at different time points
- Fidelity decay and quantum Loschmidt echoes
- Visualization of distance matrices for detecting patterns

### 4. `quantum.ipynb`

Basic introduction to quantum simulation techniques:
- Simple quantum state manipulations
- Tensor product structure
- Basic operator actions
- Getting started with the simulation framework

## Usage

To run these notebooks:

```bash
# Make sure Jupyter is installed
pip install jupyter

# Launch Jupyter Lab/Notebook
jupyter lab
# or
jupyter notebook
```

Navigate to the notebook of interest and execute the cells to see the demonstrations and visualizations.

## Integration with Other Modules

These notebooks leverage functionality from the main repository modules:
- **emergence**: Used for Page-Wootters mechanism simulations
- **ising**: Used for Ising model simulations
- **syk**: Used for SYK model chaos analysis

See the main repository README for more information on these components.

## License

MIT License 