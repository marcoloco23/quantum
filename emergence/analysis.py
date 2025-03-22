import numpy as np
from scipy.fft import fft, fftfreq
from scipy.signal import correlate
from typing import Dict, List, Tuple, Optional, Union, Any, Callable

"""
Analysis module for quantum simulations.
"""


###############################################################################
# Signal processing utilities
###############################################################################
def normalized_autocorr(x: np.ndarray) -> np.ndarray:
    """Compute normalized autocorrelation of a signal.

    Args:
        x: Input signal

    Returns:
        Normalized autocorrelation function
    """
    # Zero-mean the signal
    x_centered = x - np.mean(x)

    # Compute autocorrelation
    result = correlate(x_centered, x_centered, mode="full")

    # Normalize to have autocorrelation=1 at lag=0 and take positive lags
    return result[len(result) // 2 :] / result[len(result) // 2]


def compute_fft(signal: np.ndarray, dt: float = 0.2) -> Tuple[np.ndarray, np.ndarray]:
    """Compute FFT of a signal.

    Args:
        signal: Input time-domain signal
        dt: Time step

    Returns:
        Tuple of (frequencies, amplitudes)
    """
    # Length of signal
    n = len(signal)

    # Compute FFT
    fft_vals = fft(signal)

    # Compute frequency axis
    freq = fftfreq(n, d=dt)

    # Only positive frequencies up to Nyquist frequency
    positive_freq_idx = np.arange(1, n // 2)

    # Extract frequencies and amplitudes
    frequencies = freq[positive_freq_idx]
    amplitudes = 2.0 / n * np.abs(fft_vals[positive_freq_idx])

    return frequencies, amplitudes


def compute_distance_matrix(signal: np.ndarray) -> np.ndarray:
    """Compute distance matrix for recurrence analysis.

    Args:
        signal: Input time-domain signal

    Returns:
        Distance matrix
    """
    n = len(signal)
    dist_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(n):
            dist_matrix[i, j] = abs(signal[i] - signal[j])

    return dist_matrix


def find_recurrence_times(ac: np.ndarray, threshold: float = 0.5) -> List[int]:
    """Find recurrence times from autocorrelation function.

    Args:
        ac: Autocorrelation function
        threshold: Threshold for peak detection

    Returns:
        List of indices corresponding to recurrence times
    """
    peaks = []
    for i in range(2, len(ac) - 1):
        if ac[i] > threshold and ac[i] > ac[i - 1] and ac[i] > ac[i + 1]:
            peaks.append(i)
    return peaks


###############################################################################
# Statistical analysis
###############################################################################
def compute_average_metrics(results_list: List[Dict], metric_key: str) -> np.ndarray:
    """Compute average of a metric across multiple simulations.

    Args:
        results_list: List of simulation result dictionaries
        metric_key: Key of the metric to average

    Returns:
        Average values
    """
    # Extract values for the specified metric
    metric_values = [results[metric_key] for results in results_list]

    # Convert to numpy array and average
    return np.mean(np.array(metric_values), axis=0)


def compute_metric_statistics(
    results: Dict[str, np.ndarray], metric_key: str
) -> Dict[str, float]:
    """Compute statistics for a metric.

    Args:
        results: Simulation results dictionary
        metric_key: Key of the metric to analyze

    Returns:
        Dictionary of statistics
    """
    metric = results[metric_key]

    stats = {
        "mean": np.mean(metric),
        "std": np.std(metric),
        "min": np.min(metric),
        "max": np.max(metric),
        "initial": metric[0] if len(metric) > 0 else None,
        "final": metric[-1] if len(metric) > 0 else None,
    }

    return stats


def compare_metrics(
    results_dict: Dict[str, Dict[str, np.ndarray]], metric_key: str
) -> Dict[str, Dict[str, float]]:
    """Compare a metric across multiple simulations.

    Args:
        results_dict: Dictionary mapping simulation names to result dictionaries
        metric_key: Key of the metric to compare

    Returns:
        Dictionary mapping simulation names to metric statistics
    """
    comparison = {}

    for sim_name, results in results_dict.items():
        comparison[sim_name] = compute_metric_statistics(results, metric_key)

    return comparison


###############################################################################
# Frequency analysis
###############################################################################
def analyze_frequency_components(
    results_xx: Dict[str, np.ndarray],
    results_zz: Dict[str, np.ndarray],
    results_mixed: Dict[str, np.ndarray],
    dt: float = 0.2,
) -> Dict[str, Dict[str, Tuple[np.ndarray, np.ndarray]]]:
    """Analyze frequency components of coherence and purity signals.

    Args:
        results_xx: Results from xx coupling simulation
        results_zz: Results from zz coupling simulation
        results_mixed: Results from mixed coupling simulation
        dt: Time step

    Returns:
        Dictionary of frequency analysis results
    """
    # Time domain signals
    coherence_xx = results_xx["coherence1"]
    coherence_zz = results_zz["coherence1"]
    coherence_mixed = results_mixed["coherence1"]

    purity_xx = results_xx["purity1"]
    purity_zz = results_zz["purity1"]
    purity_mixed = results_mixed["purity1"]

    # Compute FFTs
    freq_coherence_xx, fft_coherence_xx = compute_fft(coherence_xx, dt)
    freq_coherence_zz, fft_coherence_zz = compute_fft(coherence_zz, dt)
    freq_coherence_mixed, fft_coherence_mixed = compute_fft(coherence_mixed, dt)

    freq_purity_xx, fft_purity_xx = compute_fft(purity_xx, dt)
    freq_purity_zz, fft_purity_zz = compute_fft(purity_zz, dt)
    freq_purity_mixed, fft_purity_mixed = compute_fft(purity_mixed, dt)

    # Return all results
    return {
        "coherence": {
            "xx": (freq_coherence_xx, fft_coherence_xx),
            "zz": (freq_coherence_zz, fft_coherence_zz),
            "mixed": (freq_coherence_mixed, fft_coherence_mixed),
        },
        "purity": {
            "xx": (freq_purity_xx, fft_purity_xx),
            "zz": (freq_purity_zz, fft_purity_zz),
            "mixed": (freq_purity_mixed, fft_purity_mixed),
        },
    }


###############################################################################
# Recurrence analysis
###############################################################################
def analyze_recurrence_dynamics(
    results_xx: Dict[str, np.ndarray],
    results_zz: Dict[str, np.ndarray],
    results_mixed: Dict[str, np.ndarray],
) -> Dict[str, Any]:
    """Analyze recurrence dynamics and temporal correlations.

    Args:
        results_xx: Results from xx coupling simulation
        results_zz: Results from zz coupling simulation
        results_mixed: Results from mixed coupling simulation

    Returns:
        Dictionary of recurrence analysis results
    """
    # Extract time series data
    coherence_xx = results_xx["coherence1"]
    coherence_zz = results_zz["coherence1"]
    coherence_mixed = results_mixed["coherence1"]

    purity_xx = results_xx["purity1"]
    purity_zz = results_zz["purity1"]
    purity_mixed = results_mixed["purity1"]

    entropy_xx = results_xx["entropy"]
    entropy_zz = results_zz["entropy"]
    entropy_mixed = results_mixed["entropy"]

    # Compute auto-correlation functions
    ac_coherence_xx = normalized_autocorr(coherence_xx)
    ac_coherence_zz = normalized_autocorr(coherence_zz)
    ac_coherence_mixed = normalized_autocorr(coherence_mixed)

    ac_purity_xx = normalized_autocorr(purity_xx)
    ac_purity_zz = normalized_autocorr(purity_zz)
    ac_purity_mixed = normalized_autocorr(purity_mixed)

    ac_entropy_xx = normalized_autocorr(entropy_xx)
    ac_entropy_zz = normalized_autocorr(entropy_zz)
    ac_entropy_mixed = normalized_autocorr(entropy_mixed)

    # Compute distance matrices for recurrence plots
    dist_coherence_xx = compute_distance_matrix(coherence_xx)
    dist_purity_xx = compute_distance_matrix(purity_xx)
    dist_entropy_xx = compute_distance_matrix(entropy_xx)

    # Find recurrence times
    rec_times_coherence_xx = find_recurrence_times(ac_coherence_xx)
    rec_times_purity_xx = find_recurrence_times(ac_purity_xx)
    rec_times_entropy_xx = find_recurrence_times(ac_entropy_xx)

    # Return all results
    return {
        "autocorrelation": {
            "coherence": {
                "xx": ac_coherence_xx,
                "zz": ac_coherence_zz,
                "mixed": ac_coherence_mixed,
            },
            "purity": {
                "xx": ac_purity_xx,
                "zz": ac_purity_zz,
                "mixed": ac_purity_mixed,
            },
            "entropy": {
                "xx": ac_entropy_xx,
                "zz": ac_entropy_zz,
                "mixed": ac_entropy_mixed,
            },
        },
        "distance_matrices": {
            "coherence_xx": dist_coherence_xx,
            "purity_xx": dist_purity_xx,
            "entropy_xx": dist_entropy_xx,
        },
        "recurrence_times": {
            "coherence_xx": rec_times_coherence_xx,
            "purity_xx": rec_times_purity_xx,
            "entropy_xx": rec_times_entropy_xx,
        },
    }


###############################################################################
# Parameter sweep analysis
###############################################################################
def analyze_parameter_sweep(
    results_grid: np.ndarray,
    param1_values: np.ndarray,
    param2_values: np.ndarray,
    metric_key: str,
) -> Tuple[np.ndarray, Tuple[int, int]]:
    """Analyze results from a parameter sweep.

    Args:
        results_grid: 2D grid of simulation results
        param1_values: Values for first parameter
        param2_values: Values for second parameter
        metric_key: Key of the metric to analyze

    Returns:
        Tuple of (average metric values, optimal parameter indices)
    """
    # Initialize result array
    avg_values = np.zeros((len(param1_values), len(param2_values)))

    # Compute average values for each parameter combination
    for i, _ in enumerate(param1_values):
        for j, _ in enumerate(param2_values):
            result = results_grid[i, j]
            avg_values[i, j] = np.mean(result[metric_key])

    # Find optimal parameters (maximum average value)
    opt_i, opt_j = np.unravel_index(np.argmax(avg_values), avg_values.shape)

    return avg_values, (opt_i, opt_j)
