import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from typing import Dict, List, Tuple, Optional, Union, Any

"""
Visualization module for quantum simulation results.
"""


def setup_image_path(filename: str) -> str:
    """Setup the path for saving images.

    Args:
        filename: Image filename

    Returns:
        Complete image path
    """
    import os

    # Create images directory if it doesn't exist
    os.makedirs("images", exist_ok=True)

    return os.path.join("images", filename)


###############################################################################
# Basic Plotting Functions
###############################################################################
def plot_time_series(
    time_values: np.ndarray,
    data_series: Dict[str, np.ndarray],
    title: str,
    xlabel: str = "Time",
    ylabel: str = "Value",
    figsize: Tuple[int, int] = (10, 6),
    save_path: Optional[str] = None,
    show: bool = True,
    line_styles: Optional[Dict[str, str]] = None,
):
    """Plot time series data.

    Args:
        time_values: Array of time values
        data_series: Dictionary mapping series names to data values
        title: Plot title
        xlabel: x-axis label
        ylabel: y-axis label
        figsize: Figure size
        save_path: Path to save the figure
        show: Whether to show the figure
        line_styles: Dictionary mapping series names to line styles
    """
    fig, ax = plt.subplots(figsize=figsize)

    if line_styles is None:
        line_styles = {}

    for name, data in data_series.items():
        style = line_styles.get(name, "-")
        ax.plot(time_values, data, style, label=name)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True)
    ax.legend()

    if save_path:
        plt.savefig(setup_image_path(save_path), dpi=300)

    if show:
        plt.show()
    else:
        plt.close(fig)


def plot_heatmap(
    data: np.ndarray,
    x_values: np.ndarray,
    y_values: np.ndarray,
    title: str,
    xlabel: str,
    ylabel: str,
    colorbar_label: str = "Value",
    figsize: Tuple[int, int] = (8, 6),
    save_path: Optional[str] = None,
    show: bool = True,
    cmap: str = "viridis",
):
    """Plot heatmap.

    Args:
        data: 2D array of data
        x_values: Array of x-axis values
        y_values: Array of y-axis values
        title: Plot title
        xlabel: x-axis label
        ylabel: y-axis label
        colorbar_label: Label for colorbar
        figsize: Figure size
        save_path: Path to save the figure
        show: Whether to show the figure
        cmap: Colormap name
    """
    fig, ax = plt.subplots(figsize=figsize)

    im = ax.imshow(
        data,
        origin="lower",
        extent=[x_values[0], x_values[-1], y_values[0], y_values[-1]],
        aspect="auto",
        cmap=cmap,
    )

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(colorbar_label)

    if save_path:
        plt.savefig(setup_image_path(save_path), dpi=300)

    if show:
        plt.show()
    else:
        plt.close(fig)


def plot_surface_3d(
    data: np.ndarray,
    x_values: np.ndarray,
    y_values: np.ndarray,
    title: str,
    xlabel: str,
    ylabel: str,
    zlabel: str,
    figsize: Tuple[int, int] = (10, 8),
    save_path: Optional[str] = None,
    show: bool = True,
    cmap: str = "viridis",
    azimuth: float = -60,
    elevation: float = 30,
):
    """Plot 3D surface.

    Args:
        data: 2D array of data
        x_values: Array of x-axis values
        y_values: Array of y-axis values
        title: Plot title
        xlabel: x-axis label
        ylabel: y-axis label
        zlabel: z-axis label
        figsize: Figure size
        save_path: Path to save the figure
        show: Whether to show the figure
        cmap: Colormap name
        azimuth: Azimuth angle for view
        elevation: Elevation angle for view
    """
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection="3d")

    X, Y = np.meshgrid(x_values, y_values)
    surf = ax.plot_surface(X, Y, data.T, cmap=cm.get_cmap(cmap), antialiased=True)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    ax.set_title(title)

    ax.view_init(elevation, azimuth)

    fig.colorbar(surf, ax=ax, shrink=0.5)

    if save_path:
        plt.savefig(setup_image_path(save_path), dpi=300)

    if show:
        plt.show()
    else:
        plt.close(fig)


###############################################################################
# Specific Visualization Functions
###############################################################################
def plot_simulation_results(
    results: Dict[str, np.ndarray],
    time_values: Optional[np.ndarray] = None,
    save_path: Optional[str] = None,
    show: bool = True,
):
    """Plot simulation results.

    Args:
        results: Dictionary of simulation results
        time_values: Array of time values (default: range of clock states)
        save_path: Base path for saving figures
        show: Whether to show the figures
    """
    if time_values is None:
        time_values = np.arange(len(results["clock_probs"]))

    # Plot clock state probabilities
    plot_time_series(
        time_values,
        {"Clock Probability": results["clock_probs"]},
        "Clock State Probabilities",
        "Clock Time T",
        "Probability",
        save_path=f"{save_path}_clock_probs.png" if save_path else None,
        show=show,
    )

    # Plot observables
    plot_time_series(
        time_values,
        {
            "Qubit 1 ⟨σz⟩": results["exp_values_Z1"],
            "Qubit 2 ⟨σz⟩": results["exp_values_Z2"],
        },
        "System Observables",
        "Clock Time T",
        "Expectation Value",
        save_path=f"{save_path}_observables.png" if save_path else None,
        show=show,
    )

    # Plot coherence
    plot_time_series(
        time_values,
        {"Qubit 1": results["coherence1"], "Qubit 2": results["coherence2"]},
        "Coherence |ρ01|",
        "Clock Time T",
        "Coherence",
        save_path=f"{save_path}_coherence.png" if save_path else None,
        show=show,
    )

    # Plot purity
    plot_time_series(
        time_values,
        {
            "Qubit 1": results["purity1"],
            "Qubit 2": results["purity2"],
            "System": results.get("system_purity", np.ones_like(results["purity1"])),
        },
        "Purity Tr(ρ²)",
        "Clock Time T",
        "Purity",
        save_path=f"{save_path}_purity.png" if save_path else None,
        show=show,
    )

    # Plot entropy
    plot_time_series(
        time_values,
        {"Entropy": results["entropy"]},
        "System-Environment Entanglement Entropy",
        "Clock Time T",
        "Entropy",
        save_path=f"{save_path}_entropy.png" if save_path else None,
        show=show,
    )


def plot_coupling_strength_sweep(
    sweep_results: Dict[str, Any], save_path: Optional[str] = None, show: bool = True
):
    """Plot results from coupling strength sweep.

    Args:
        sweep_results: Dictionary of sweep results
        save_path: Path for saving figure
        show: Whether to show the figure
    """
    g_env1_values = sweep_results["g_env1_values"]
    g_env2_values = sweep_results["g_env2_values"]
    avg_coherence1 = sweep_results["avg_coherence1"]
    avg_purity1 = sweep_results["avg_purity1"]
    avg_entropy = sweep_results["avg_entropy"]

    fig = plt.figure(figsize=(18, 12))

    # 3D surface plot for average coherence
    ax1 = fig.add_subplot(231, projection="3d")
    X, Y = np.meshgrid(g_env1_values, g_env2_values)
    surf = ax1.plot_surface(X, Y, avg_coherence1.T, cmap=cm.viridis)
    ax1.set_xlabel("g_env1")
    ax1.set_ylabel("g_env2")
    ax1.set_zlabel("Avg Coherence (Qubit 1)")
    ax1.set_title("Average Coherence vs Coupling Strengths")
    fig.colorbar(surf, ax=ax1, shrink=0.5)

    # 3D surface plot for average purity
    ax2 = fig.add_subplot(232, projection="3d")
    surf = ax2.plot_surface(X, Y, avg_purity1.T, cmap=cm.viridis)
    ax2.set_xlabel("g_env1")
    ax2.set_ylabel("g_env2")
    ax2.set_zlabel("Avg Purity (Qubit 1)")
    ax2.set_title("Average Purity vs Coupling Strengths")
    fig.colorbar(surf, ax=ax2, shrink=0.5)

    # 3D surface plot for average entropy
    ax3 = fig.add_subplot(233, projection="3d")
    surf = ax3.plot_surface(X, Y, avg_entropy.T, cmap=cm.viridis)
    ax3.set_xlabel("g_env1")
    ax3.set_ylabel("g_env2")
    ax3.set_zlabel("Avg Entropy")
    ax3.set_title("Average Entanglement Entropy vs Coupling Strengths")
    fig.colorbar(surf, ax=ax3, shrink=0.5)

    # 2D heatmap for average coherence
    ax4 = fig.add_subplot(234)
    im = ax4.imshow(
        avg_coherence1,
        origin="lower",
        extent=[
            g_env1_values[0],
            g_env1_values[-1],
            g_env2_values[0],
            g_env2_values[-1],
        ],
    )
    ax4.set_xlabel("g_env1")
    ax4.set_ylabel("g_env2")
    ax4.set_title("Average Coherence Heatmap (Qubit 1)")
    fig.colorbar(im, ax=ax4)

    # 2D heatmap for average purity
    ax5 = fig.add_subplot(235)
    im = ax5.imshow(
        avg_purity1,
        origin="lower",
        extent=[
            g_env1_values[0],
            g_env1_values[-1],
            g_env2_values[0],
            g_env2_values[-1],
        ],
    )
    ax5.set_xlabel("g_env1")
    ax5.set_ylabel("g_env2")
    ax5.set_title("Average Purity Heatmap (Qubit 1)")
    fig.colorbar(im, ax=ax5)

    # 2D heatmap for average entropy
    ax6 = fig.add_subplot(236)
    im = ax6.imshow(
        avg_entropy,
        origin="lower",
        extent=[
            g_env1_values[0],
            g_env1_values[-1],
            g_env2_values[0],
            g_env2_values[-1],
        ],
    )
    ax6.set_xlabel("g_env1")
    ax6.set_ylabel("g_env2")
    ax6.set_title("Average Entropy Heatmap")
    fig.colorbar(im, ax=ax6)

    plt.tight_layout()

    if save_path:
        plt.savefig(setup_image_path(save_path), dpi=300)

    if show:
        plt.show()
    else:
        plt.close(fig)


def plot_coupling_type_comparison(
    results_xx: Dict[str, np.ndarray],
    results_zz: Dict[str, np.ndarray],
    results_mixed: Dict[str, np.ndarray],
    g_env: float = 0.05,
    save_path: Optional[str] = None,
    show: bool = True,
):
    """Plot comparison of different coupling types.

    Args:
        results_xx: Results from xx coupling simulation
        results_zz: Results from zz coupling simulation
        results_mixed: Results from mixed coupling simulation
        g_env: Environment coupling strength
        save_path: Path for saving figure
        show: Whether to show the figure
    """
    fig, axes = plt.subplots(3, 1, figsize=(10, 12))

    # Time values
    T_values = np.arange(len(results_xx["coherence1"]))

    # Plot coherence
    ax1 = axes[0]
    ax1.plot(
        T_values, results_xx["coherence1"], "b-o", label="XX Coupling", markersize=4
    )
    ax1.plot(
        T_values, results_zz["coherence1"], "r-o", label="ZZ Coupling", markersize=4
    )
    ax1.plot(
        T_values,
        results_mixed["coherence1"],
        "g-o",
        label="Mixed Coupling",
        markersize=4,
    )
    ax1.set_xlabel("Clock Time T")
    ax1.set_ylabel("Coherence |ρ01|")
    ax1.set_title(f"Qubit 1 Coherence vs Clock Time (g_env={g_env:.2f})")
    ax1.grid(True)
    ax1.legend()

    # Plot purity
    ax2 = axes[1]
    ax2.plot(T_values, results_xx["purity1"], "b-o", label="XX Coupling", markersize=4)
    ax2.plot(T_values, results_zz["purity1"], "r-o", label="ZZ Coupling", markersize=4)
    ax2.plot(
        T_values, results_mixed["purity1"], "g-o", label="Mixed Coupling", markersize=4
    )
    ax2.set_xlabel("Clock Time T")
    ax2.set_ylabel("Purity Tr(ρ²)")
    ax2.set_title(f"Qubit 1 Purity vs Clock Time (g_env={g_env:.2f})")
    ax2.grid(True)
    ax2.legend()

    # Plot entropy
    ax3 = axes[2]
    ax3.plot(T_values, results_xx["entropy"], "b-o", label="XX Coupling", markersize=4)
    ax3.plot(T_values, results_zz["entropy"], "r-o", label="ZZ Coupling", markersize=4)
    ax3.plot(
        T_values, results_mixed["entropy"], "g-o", label="Mixed Coupling", markersize=4
    )
    ax3.set_xlabel("Clock Time T")
    ax3.set_ylabel("Entanglement Entropy")
    ax3.set_title(f"System-Environment Entanglement Entropy (g_env={g_env:.2f})")
    ax3.grid(True)
    ax3.legend()

    plt.tight_layout()

    if save_path:
        plt.savefig(setup_image_path(save_path), dpi=300)

    if show:
        plt.show()
    else:
        plt.close(fig)


def plot_frequency_analysis(
    freq_results: Dict[str, Dict[str, Tuple[np.ndarray, np.ndarray]]],
    save_path: Optional[str] = None,
    show: bool = True,
):
    """Plot frequency analysis results.

    Args:
        freq_results: Dictionary of frequency analysis results
        save_path: Path for saving figure
        show: Whether to show the figure
    """
    fig, axes = plt.subplots(2, 1, figsize=(12, 10))

    # Plot coherence frequency components
    ax1 = axes[0]

    for coupling_type, (freq, amp) in freq_results["coherence"].items():
        color = "b" if coupling_type == "xx" else "r" if coupling_type == "zz" else "g"
        ax1.plot(freq, amp, f"{color}-", label=f"{coupling_type.upper()} Coupling")

    ax1.set_xlabel("Frequency (1/time)")
    ax1.set_ylabel("Amplitude")
    ax1.set_title("Frequency Components of Coherence Signal")
    ax1.grid(True)
    ax1.legend()

    # Plot purity frequency components
    ax2 = axes[1]

    for coupling_type, (freq, amp) in freq_results["purity"].items():
        color = "b" if coupling_type == "xx" else "r" if coupling_type == "zz" else "g"
        ax2.plot(freq, amp, f"{color}-", label=f"{coupling_type.upper()} Coupling")

    ax2.set_xlabel("Frequency (1/time)")
    ax2.set_ylabel("Amplitude")
    ax2.set_title("Frequency Components of Purity Signal")
    ax2.grid(True)
    ax2.legend()

    plt.tight_layout()

    if save_path:
        plt.savefig(setup_image_path(save_path), dpi=300)

    if show:
        plt.show()
    else:
        plt.close(fig)
