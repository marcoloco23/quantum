#!/usr/bin/env python3
import argparse
import numpy as np
from numpy import kron
import matplotlib.pyplot as plt
import os
import sys
from typing import Dict, List, Tuple, Optional, Any

from quantum_core import sigma_x, sigma_z, I2, zero, plus, bell_plus
from quantum_simulator import (
    SystemConfiguration,
    QuantumSimulator,
    run_simulation_with_custom_state,
    run_simulation_with_env_size,
)
from analysis import (
    analyze_frequency_components,
    analyze_recurrence_dynamics,
    analyze_parameter_sweep,
)
from visualization import (
    setup_image_path,
    plot_simulation_results,
    plot_coupling_strength_sweep,
    plot_coupling_type_comparison,
    plot_frequency_analysis,
)

"""
Main script for emergent time quantum simulations.
"""


def create_parser() -> argparse.ArgumentParser:
    """Create command-line argument parser.

    Returns:
        Argument parser
    """
    parser = argparse.ArgumentParser(
        description="Emergent Time with Environment Coupling Analysis"
    )

    # Add subparsers for different analyses
    subparsers = parser.add_subparsers(
        dest="analysis_type", help="Type of analysis to run"
    )

    # 1. Coupling strength sweep
    sweep_parser = subparsers.add_parser(
        "sweep", help="Sweep environment coupling strengths"
    )
    sweep_parser.add_argument(
        "--coupling-type",
        type=str,
        default="xx",
        choices=["xx", "zz", "mixed"],
        help="Type of system-environment coupling",
    )
    sweep_parser.add_argument(
        "--min-g", type=float, default=0.0, help="Minimum coupling strength"
    )
    sweep_parser.add_argument(
        "--max-g", type=float, default=0.25, help="Maximum coupling strength"
    )
    sweep_parser.add_argument(
        "--n-points", type=int, default=6, help="Number of coupling strength points"
    )

    # 2. Compare coupling types
    compare_parser = subparsers.add_parser(
        "compare", help="Compare different coupling types"
    )
    compare_parser.add_argument(
        "--g-env", type=float, default=0.05, help="Environment coupling strength"
    )

    # 3. Advanced analysis
    advanced_parser = subparsers.add_parser(
        "advanced", help="Advanced analysis (frequency/recurrence)"
    )
    advanced_parser.add_argument(
        "--g-env", type=float, default=0.05, help="Environment coupling strength"
    )

    # 4. Initial state analysis
    initial_parser = subparsers.add_parser("initial", help="Initial state analysis")
    initial_parser.add_argument(
        "--g-env", type=float, default=0.05, help="Environment coupling strength"
    )

    # 5. Environment size scaling
    env_parser = subparsers.add_parser(
        "env-scaling", help="Environment size scaling analysis"
    )
    env_parser.add_argument(
        "--max-size", type=int, default=4, help="Maximum number of environment qubits"
    )
    env_parser.add_argument(
        "--coupling-type",
        type=str,
        default="xx",
        choices=["xx", "zz"],
        help="Type of system-environment coupling",
    )

    # General options
    parser.add_argument(
        "--clock-states", type=int, default=100, help="Number of clock states"
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Do not display plots (save to file only)",
    )
    parser.add_argument(
        "--output-dir", type=str, default="images", help="Output directory for plots"
    )

    return parser


def coupling_strength_sweep(args):
    """Run coupling strength sweep analysis.

    Args:
        args: Command-line arguments
    """
    print(f"Running coupling strength sweep (coupling type: {args.coupling_type})...")

    # Define parameter ranges
    g_env1_values = np.linspace(args.min_g, args.max_g, args.n_points)
    g_env2_values = np.linspace(args.min_g, args.max_g, args.n_points)

    # Initialize result arrays
    avg_coherence1 = np.zeros((len(g_env1_values), len(g_env2_values)))
    avg_coherence2 = np.zeros((len(g_env1_values), len(g_env2_values)))
    avg_purity1 = np.zeros((len(g_env1_values), len(g_env2_values)))
    avg_purity2 = np.zeros((len(g_env1_values), len(g_env2_values)))
    avg_entropy = np.zeros((len(g_env1_values), len(g_env2_values)))

    # Run simulations for each parameter combination
    results_grid = np.empty((len(g_env1_values), len(g_env2_values)), dtype=object)

    for i, g1 in enumerate(g_env1_values):
        for j, g2 in enumerate(g_env2_values):
            print(f"  Running simulation with g_env1={g1:.2f}, g_env2={g2:.2f}")

            # Create configuration
            config = SystemConfiguration(
                num_clock_states=args.clock_states,
                env_coupling1=g1,
                env_coupling2=g2,
                coupling_type=args.coupling_type,
            )

            # Create simulator and run
            simulator = QuantumSimulator(config)
            results = simulator.run_simulation()
            results_grid[i, j] = results

            # Store average values
            avg_coherence1[i, j] = np.mean(results["coherence1"])
            avg_coherence2[i, j] = np.mean(results["coherence2"])
            avg_purity1[i, j] = np.mean(results["purity1"])
            avg_purity2[i, j] = np.mean(results["purity2"])
            avg_entropy[i, j] = np.mean(results["entropy"])

    # Find optimal parameters for coherence
    _, (opt_i, opt_j) = analyze_parameter_sweep(
        results_grid, g_env1_values, g_env2_values, "coherence1"
    )
    print(
        f"Optimal coupling strengths for coherence: g1={g_env1_values[opt_i]:.2f}, g2={g_env2_values[opt_j]:.2f}"
    )

    # Prepare results for plotting
    sweep_results = {
        "g_env1_values": g_env1_values,
        "g_env2_values": g_env2_values,
        "avg_coherence1": avg_coherence1,
        "avg_coherence2": avg_coherence2,
        "avg_purity1": avg_purity1,
        "avg_purity2": avg_purity2,
        "avg_entropy": avg_entropy,
    }

    # Plot results
    plot_coupling_strength_sweep(
        sweep_results,
        save_path=f"coupling_sweep_results_{args.coupling_type}.png",
        show=not args.no_plots,
    )

    print(
        f"Sweep complete. Results saved to '{setup_image_path(f'coupling_sweep_results_{args.coupling_type}.png')}'"
    )

    return sweep_results


def compare_coupling_types(args):
    """Compare different coupling types.

    Args:
        args: Command-line arguments
    """
    print(f"Comparing coupling types (g_env={args.g_env:.2f})...")

    # Run simulations for each coupling type
    g_env1 = args.g_env
    g_env2 = args.g_env

    # XX coupling
    print("  Running XX coupling simulation...")
    config_xx = SystemConfiguration(
        num_clock_states=args.clock_states,
        env_coupling1=g_env1,
        env_coupling2=g_env2,
        coupling_type="xx",
    )
    simulator_xx = QuantumSimulator(config_xx)
    results_xx = simulator_xx.run_simulation()

    # ZZ coupling
    print("  Running ZZ coupling simulation...")
    config_zz = SystemConfiguration(
        num_clock_states=args.clock_states,
        env_coupling1=g_env1,
        env_coupling2=g_env2,
        coupling_type="zz",
    )
    simulator_zz = QuantumSimulator(config_zz)
    results_zz = simulator_zz.run_simulation()

    # Mixed coupling
    print("  Running mixed coupling simulation...")
    config_mixed = SystemConfiguration(
        num_clock_states=args.clock_states,
        env_coupling1=g_env1,
        env_coupling2=g_env2,
        coupling_type="mixed",
    )
    simulator_mixed = QuantumSimulator(config_mixed)
    results_mixed = simulator_mixed.run_simulation()

    # Plot results
    plot_coupling_type_comparison(
        results_xx,
        results_zz,
        results_mixed,
        g_env=g_env1,
        save_path="coupling_type_comparison.png",
        show=not args.no_plots,
    )

    print(
        f"Comparison complete. Results saved to '{setup_image_path('coupling_type_comparison.png')}'"
    )

    return results_xx, results_zz, results_mixed


def run_advanced_analysis(args):
    """Run advanced analysis (frequency and recurrence dynamics).

    Args:
        args: Command-line arguments
    """
    print(f"Running advanced analysis (g_env={args.g_env:.2f})...")

    # Get results from coupling type comparison
    results_xx, results_zz, results_mixed = compare_coupling_types(args)

    # Frequency analysis
    print("  Analyzing frequency components...")
    freq_results = analyze_frequency_components(
        results_xx, results_zz, results_mixed, dt=0.2
    )

    # Plot frequency results
    plot_frequency_analysis(
        freq_results, save_path="frequency_analysis.png", show=not args.no_plots
    )

    # Recurrence analysis
    print("  Analyzing recurrence dynamics...")
    recurrence_info = analyze_recurrence_dynamics(results_xx, results_zz, results_mixed)

    # Print recurrence times
    print("Recurrence analysis results:")
    print(
        f"Coherence recurrence times (XX coupling): {recurrence_info['recurrence_times']['coherence_xx']}"
    )
    print(
        f"Purity recurrence times (XX coupling): {recurrence_info['recurrence_times']['purity_xx']}"
    )
    print(
        f"Entropy recurrence times (XX coupling): {recurrence_info['recurrence_times']['entropy_xx']}"
    )

    print(
        f"Advanced analysis complete. Results saved to '{setup_image_path('frequency_analysis.png')}'"
    )

    return freq_results, recurrence_info


def run_initial_state_analysis(args):
    """Analyze different initial states.

    Args:
        args: Command-line arguments
    """
    print(f"Analyzing different initial states (g_env={args.g_env:.2f})...")

    # Define different initial states
    # State 1: |+>⊗|+> (superposition)
    psi0_plus_plus = kron(plus, plus)

    # State 2: |0>⊗|0> (computational basis)
    psi0_zero_zero = kron(zero, zero)

    # State 3: |+>⊗|0> (mixed)
    psi0_plus_zero = kron(plus, zero)

    # State 4: Bell state (|00> + |11>)/√2 (entangled)
    psi0_bell = bell_plus

    # Run simulations for each initial state
    print("  Running simulation with |+>⊗|+> initial state...")
    results_plus_plus = run_simulation_with_custom_state(
        psi0_plus_plus,
        env_coupling1=args.g_env,
        env_coupling2=args.g_env,
        coupling_type="xx",
        num_clock_states=args.clock_states,
    )

    print("  Running simulation with |0>⊗|0> initial state...")
    results_zero_zero = run_simulation_with_custom_state(
        psi0_zero_zero,
        env_coupling1=args.g_env,
        env_coupling2=args.g_env,
        coupling_type="xx",
        num_clock_states=args.clock_states,
    )

    print("  Running simulation with |+>⊗|0> initial state...")
    results_plus_zero = run_simulation_with_custom_state(
        psi0_plus_zero,
        env_coupling1=args.g_env,
        env_coupling2=args.g_env,
        coupling_type="xx",
        num_clock_states=args.clock_states,
    )

    print("  Running simulation with Bell state...")
    results_bell = run_simulation_with_custom_state(
        psi0_bell,
        env_coupling1=args.g_env,
        env_coupling2=args.g_env,
        coupling_type="xx",
        num_clock_states=args.clock_states,
    )

    # Plot comparative results
    fig, axes = plt.subplots(3, 1, figsize=(12, 15))

    # Time values
    T_values = np.arange(args.clock_states)

    # Plot coherence
    ax1 = axes[0]
    ax1.plot(
        T_values, results_plus_plus["coherence1"], "b-", label="|+>⊗|+>", linewidth=1.5
    )
    ax1.plot(
        T_values, results_zero_zero["coherence1"], "r-", label="|0>⊗|0>", linewidth=1.5
    )
    ax1.plot(
        T_values, results_plus_zero["coherence1"], "g-", label="|+>⊗|0>", linewidth=1.5
    )
    ax1.plot(
        T_values, results_bell["coherence1"], "k-", label="Bell state", linewidth=1.5
    )
    ax1.set_xlabel("Clock Time T")
    ax1.set_ylabel("Coherence |ρ01|")
    ax1.set_title(
        f"Qubit 1 Coherence vs Clock Time for Different Initial States (g_env={args.g_env:.2f})"
    )
    ax1.grid(True)
    ax1.legend()

    # Plot purity
    ax2 = axes[1]
    ax2.plot(
        T_values, results_plus_plus["purity1"], "b-", label="|+>⊗|+>", linewidth=1.5
    )
    ax2.plot(
        T_values, results_zero_zero["purity1"], "r-", label="|0>⊗|0>", linewidth=1.5
    )
    ax2.plot(
        T_values, results_plus_zero["purity1"], "g-", label="|+>⊗|0>", linewidth=1.5
    )
    ax2.plot(T_values, results_bell["purity1"], "k-", label="Bell state", linewidth=1.5)
    ax2.set_xlabel("Clock Time T")
    ax2.set_ylabel("Purity Tr(ρ²)")
    ax2.set_title(
        f"Qubit 1 Purity vs Clock Time for Different Initial States (g_env={args.g_env:.2f})"
    )
    ax2.grid(True)
    ax2.legend()

    # Plot entropy
    ax3 = axes[2]
    ax3.plot(
        T_values, results_plus_plus["entropy"], "b-", label="|+>⊗|+>", linewidth=1.5
    )
    ax3.plot(
        T_values, results_zero_zero["entropy"], "r-", label="|0>⊗|0>", linewidth=1.5
    )
    ax3.plot(
        T_values, results_plus_zero["entropy"], "g-", label="|+>⊗|0>", linewidth=1.5
    )
    ax3.plot(T_values, results_bell["entropy"], "k-", label="Bell state", linewidth=1.5)
    ax3.set_xlabel("Clock Time T")
    ax3.set_ylabel("Entanglement Entropy")
    ax3.set_title(
        f"System-Environment Entanglement Entropy for Different Initial States (g_env={args.g_env:.2f})"
    )
    ax3.grid(True)
    ax3.legend()

    plt.tight_layout()
    plt.savefig(setup_image_path("initial_state_comparison.png"), dpi=300)

    if not args.no_plots:
        plt.show()
    else:
        plt.close(fig)

    print(
        f"Initial state analysis complete. Results saved to '{setup_image_path('initial_state_comparison.png')}'"
    )

    return {
        "results_plus_plus": results_plus_plus,
        "results_zero_zero": results_zero_zero,
        "results_plus_zero": results_plus_zero,
        "results_bell": results_bell,
    }


def run_environment_scaling(args):
    """Analyze how environment size affects decoherence and recurrence.

    Args:
        args: Command-line arguments
    """
    print(f"Running environment size scaling analysis (max size: {args.max_size})...")
    print("WARNING: This analysis requires significant computational resources.")
    print("The simulation time increases exponentially with environment size.")

    # Environment sizes to test
    env_sizes = list(range(1, args.max_size + 1))  # 1 to max_size environment qubits

    results_dict = {}

    # Run simulations for each environment size
    for n_env in env_sizes:
        print(f"  Running simulation with {n_env} environment qubits...")
        results = run_simulation_with_env_size(
            n_env,
            coupling_type=args.coupling_type,
            g_env=0.05,
            num_clock_states=min(
                args.clock_states, 50
            ),  # Limit for computational efficiency
        )
        results_dict[n_env] = results

    # Create plots
    fig, axes = plt.subplots(2, 1, figsize=(12, 10))

    # Time values
    T_values = np.arange(min(args.clock_states, 50))

    # Plot coherence vs environment size
    ax1 = axes[0]
    colors = ["b", "g", "r", "k", "m", "c", "y"]

    for i, n_env in enumerate(env_sizes):
        color_idx = i % len(colors)
        ax1.plot(
            T_values,
            results_dict[n_env]["coherence1"],
            color=colors[color_idx],
            linestyle="-",
            label=f"{n_env} Env. Qubits",
            linewidth=1.5,
        )

    ax1.set_xlabel("Clock Time T")
    ax1.set_ylabel("Coherence |ρ01|")
    ax1.set_title(
        f"Qubit 1 Coherence vs Environment Size ({args.coupling_type.upper()} Coupling)"
    )
    ax1.grid(True)
    ax1.legend()

    # Plot system purity vs environment size
    ax2 = axes[1]

    for i, n_env in enumerate(env_sizes):
        color_idx = i % len(colors)
        ax2.plot(
            T_values,
            results_dict[n_env]["system_purity"],
            color=colors[color_idx],
            linestyle="-",
            label=f"{n_env} Env. Qubits",
            linewidth=1.5,
        )

    ax2.set_xlabel("Clock Time T")
    ax2.set_ylabel("System Purity Tr(ρ²)")
    ax2.set_title(
        f"System Purity vs Environment Size ({args.coupling_type.upper()} Coupling)"
    )
    ax2.grid(True)
    ax2.legend()

    plt.tight_layout()
    plt.savefig(
        setup_image_path(f"environment_scaling_{args.coupling_type}.png"), dpi=300
    )

    if not args.no_plots:
        plt.show()
    else:
        plt.close(fig)

    print(
        f"Environment scaling analysis complete. Results saved to '{setup_image_path(f'environment_scaling_{args.coupling_type}.png')}'"
    )

    return results_dict


def main():
    """Main function."""
    parser = create_parser()
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    if args.analysis_type == "sweep":
        coupling_strength_sweep(args)
    elif args.analysis_type == "compare":
        compare_coupling_types(args)
    elif args.analysis_type == "advanced":
        run_advanced_analysis(args)
    elif args.analysis_type == "initial":
        run_initial_state_analysis(args)
    elif args.analysis_type == "env-scaling":
        run_environment_scaling(args)
    else:
        # No analysis type specified, print help
        parser.print_help()


if __name__ == "__main__":
    main()
