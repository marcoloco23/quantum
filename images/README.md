# Simulation Output Images

This directory contains visualization outputs from various simulations in the Quantum Matter project.

## Directory Purpose

This directory is automatically populated with plot images when running simulations with the `--no-plots` option or when saving figures programmatically. The directory structure is preserved in Git, but the actual image files are gitignored to keep the repository size manageable.

## Common Output Files

Typical outputs include:

- `coupling_sweep_results_xx.png`: Parameter sweep for XX coupling
- `coupling_type_comparison.png`: Comparison between XX, ZZ, and mixed coupling
- `frequency_analysis.png`: Frequency component analysis
- `recurrence_analysis.png`: Recurrence dynamics analysis
- `initial_state_comparison.png`: Comparison of different initial states
- `environment_scaling_xx.png`: Analysis of environment size scaling

## Regenerating Images

To regenerate these images, run the appropriate simulation commands from the repository root:

```bash
# Examples
python -m emergence.main sweep --coupling-type xx --min-coupling 0.0 --max-coupling 1.0 --num-points 10
python -m emergence.main advanced --coupling-strength 0.05
```

## Usage in Documentation

When preparing documentation or papers, you may reference these images or copy them to permanent locations. Consider using a consistent naming scheme to help identify which parameters produced which images. 