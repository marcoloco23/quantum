import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple, Any
from tqdm import tqdm


# ---------- Helpers -------------------------------------------------
def kron(*ops: np.ndarray) -> np.ndarray:
    """Kronecker product of many small matrices."""
    out: np.ndarray = np.array([[1]], dtype=complex)
    for op in ops:
        out = np.kron(out, op)
    return out


def bell_pair() -> np.ndarray:
    """Return |Φ+> Bell pair state vector for 2 qubits."""
    v: np.ndarray = np.zeros(4, dtype=complex)
    v[0] = v[3] = 1 / np.sqrt(2)
    return v


def random_state(n_qubits: int) -> np.ndarray:
    """Haar‑random pure state on n qubits (dimension 2^n)."""
    dim: int = 2**n_qubits
    psi: np.ndarray = np.random.randn(dim) + 1j * np.random.randn(dim)
    psi /= np.linalg.norm(psi)
    return psi


def density_matrix(psi: np.ndarray) -> np.ndarray:
    """Pure state density matrix |psi><psi|."""
    return np.outer(psi, np.conjugate(psi))


def partial_trace(rho: np.ndarray, keep: List[int], n: int, dim: int = 2) -> np.ndarray:
    """
    Trace out all qubits not in `keep` (list of indices) from an n‑qubit density matrix rho.
    Returns the reduced density matrix on the subsystem indexed by `keep`.
    """
    keep = sorted(keep)
    new_shape: List[int] = [dim] * (2 * n)
    rho_reshaped: np.ndarray = rho.reshape(new_shape)
    axes_keep: List[int] = keep + [i + n for i in keep]
    axes_trace: List[int] = [i for i in range(2 * n) if i not in axes_keep]
    perm: List[int] = axes_keep + axes_trace
    rho_perm: np.ndarray = np.transpose(rho_reshaped, perm)
    k: int = len(keep)
    dimA: int = dim**k
    dim_rest: List[int] = [dim] * len(axes_trace)
    rho_perm = rho_perm.reshape(dimA, dimA, *dim_rest)
    while rho_perm.ndim > 2:
        rho_perm = np.trace(rho_perm, axis1=-2, axis2=-1)
    return rho_perm


def entanglement_entropy(rho_sub: np.ndarray) -> float:
    """Von Neumann entropy (base‑2) of a density matrix."""
    eigs: np.ndarray = np.linalg.eigvalsh(rho_sub)
    eigs = eigs[eigs > 1e-12]
    return float(-np.sum(eigs * np.log2(eigs)))


def bell_pairs_state(N: int) -> np.ndarray:
    """
    Build an N‑qubit state that is a product of Bell pairs on (0,1),(2,3),...
    N must be even.
    """
    assert N % 2 == 0, "N must be even for Bell-pair product"
    n_pairs: int = N // 2
    dim: int = 1 << N
    psi: np.ndarray = np.zeros(dim, dtype=complex)
    for pair_bits in range(1 << n_pairs):
        idx: int = 0
        for p in range(n_pairs):
            bit: int = (pair_bits >> p) & 1
            idx |= bit << (2 * p)
            idx |= bit << (2 * p + 1)
        psi[idx] = 1
    psi /= np.linalg.norm(psi)
    return psi


def project_qubits_zero(psi: np.ndarray, k: int, N: int) -> np.ndarray:
    """
    Project the first k qubits onto |0> by zeroing out amplitudes where
    any of those qubits is 1, then renormalize.
    """
    psi2: np.ndarray = psi.copy()
    dim: int = 1 << N
    mask: int = (1 << k) - 1  # bits 0..k-1
    for idx in range(dim):
        if idx & mask:
            psi2[idx] = 0.0
    norm: float = np.linalg.norm(psi2)
    if norm < 1e-16:
        return psi2
    return psi2 / norm


def reduced_density(psi: np.ndarray, L: int, N: int) -> np.ndarray:
    """Return the reduced density matrix of the first L qubits."""
    dimA: int = 1 << L
    dimB: int = 1 << (N - L)
    psi_resh: np.ndarray = psi.reshape(dimA, dimB)
    return psi_resh @ psi_resh.conj().T


# ---------- Simulation Parameters ----------------------------------
N: int = 12  # total qubits
Ls: np.ndarray = np.arange(1, N)  # block sizes from 1 to N-1

# Build the two "vacuum" states
psi_area: np.ndarray = bell_pairs_state(N)  # area‑law state
psi_vol: np.ndarray = random_state(N)  # volume‑law (Haar random)


def entropy_curve(psi: np.ndarray) -> List[float]:
    """Compute the entanglement entropy curve S(L) for a given state psi, with progress bar."""
    return [
        entanglement_entropy(reduced_density(psi, int(L), N))
        for L in tqdm(Ls, desc="Computing S(L)")
    ]


S_area: List[float] = entropy_curve(psi_area)
S_vol: List[float] = entropy_curve(psi_vol)

results: List[Tuple[int, List[float], np.ndarray]] = []
for k in tqdm([1, 2, 3], desc="Projecting qubits"):
    psi_mass: np.ndarray = project_qubits_zero(psi_vol, k, N)
    S_mass: List[float] = entropy_curve(psi_mass)
    deltaS: np.ndarray = np.array(S_vol) - np.array(S_mass)
    results.append((k, S_mass, deltaS))

# ---------- Plot: Entropy Scaling ----------------------------------
plt.figure(figsize=(8, 5))
plt.plot(Ls, S_area, "o-", label="Area‑law (Bell pairs)")
plt.plot(Ls, S_vol, "o-", label="Volume‑law (random state)")
for k, S_mass, _ in results:
    plt.plot(Ls, S_mass, "o--", label=f"After projecting {k} qubit(s)")
plt.xlabel("Block size $L$")
plt.ylabel("Entanglement entropy $S(L)$ (bits)")
plt.title(f"Entanglement entropy vs block size (N={N} qubits)")
plt.legend(loc="upper left")
plt.grid(True)
plt.tight_layout()

# ---------- Plot: Displaced Entropy -------------------------------
plt.figure(figsize=(8, 5))
for k, _, deltaS in results:
    plt.plot(Ls, deltaS, "o-", label=f"ΔS after projecting {k} qubit(s)")
plt.xlabel("Block size $L$")
plt.ylabel("Entropy displacement ΔS(L) (bits)")
plt.title('Entropy displaced by "mass" insertion')
plt.legend(loc="upper left")
plt.grid(True)
plt.tight_layout()
plt.show()
