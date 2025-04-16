import numpy as np
import matplotlib.pyplot as plt

# ---------- Helpers -------------------------------------------------
def kron(*ops):
    """Kronecker product of many small matrices."""
    out = np.array([[1]], dtype=complex)
    for op in ops:
        out = np.kron(out, op)
    return out

def bell_pair():
    """Return |Φ+> Bell pair state vector for 2 qubits."""
    v = np.zeros(4, dtype=complex)
    v[0] = v[3] = 1/np.sqrt(2)
    return v

def random_state(n_qubits):
    """Haar‑random pure state on n qubits (dimension 2^n)."""
    dim = 2**n_qubits
    psi = np.random.randn(dim) + 1j*np.random.randn(dim)
    psi /= np.linalg.norm(psi)
    return psi

def density_matrix(psi):
    """Pure state density matrix |psi><psi|."""
    return np.outer(psi, np.conjugate(psi))

def partial_trace(rho, keep, n, dim=2):
    """
    Trace out all qubits not in `keep` (list of indices) from an n‑qubit density matrix rho.
    Returns the reduced density matrix on the subsystem indexed by `keep`.
    """
    keep = sorted(keep)
    # Reshape rho into shape [dim]*2n
    new_shape = [dim] * (2*n)
    rho_reshaped = rho.reshape(new_shape)
    # Axes to keep: for each qubit i in "keep", axes i and i+n
    axes_keep = keep + [i + n for i in keep]
    axes_trace = [i for i in range(2*n) if i not in axes_keep]
    # Permute so that kept axes come first
    perm = axes_keep + axes_trace
    rho_perm = np.transpose(rho_reshaped, perm)
    k = len(keep)
    dimA = dim**k
    dim_rest = [dim] * len(axes_trace)
    # Collapse kept vs traced axes
    rho_perm = rho_perm.reshape(dimA, dimA, *dim_rest)
    # Now trace out the rest axes pairwise
    while rho_perm.ndim > 2:
        rho_perm = np.trace(rho_perm, axis1=-2, axis2=-1)
    return rho_perm

def entanglement_entropy(rho_sub):
    """Von Neumann entropy (base‑2) of a density matrix."""
    eigs = np.linalg.eigvalsh(rho_sub)
    eigs = eigs[eigs > 1e-12]
    return -np.sum(eigs * np.log2(eigs))

def bell_pairs_state(N):
    """
    Build an N‑qubit state that is a product of Bell pairs on (0,1),(2,3),...
    N must be even.
    """
    assert N % 2 == 0, "N must be even for Bell-pair product"
    n_pairs = N // 2
    dim = 1 << N
    psi = np.zeros(dim, dtype=complex)
    # For each assignment of each pair (00 or 11)
    for pair_bits in range(1 << n_pairs):
        idx = 0
        for p in range(n_pairs):
            bit = (pair_bits >> p) & 1
            # set bits 2p and 2p+1 to `bit`
            idx |= bit << (2*p)
            idx |= bit << (2*p + 1)
        psi[idx] = 1
    psi /= np.linalg.norm(psi)
    return psi

def project_qubits_zero(psi, k, N):
    """
    Project the first k qubits onto |0> by zeroing out amplitudes where
    any of those qubits is 1, then renormalize.
    """
    psi2 = psi.copy()
    dim = 1 << N
    mask = (1 << k) - 1  # bits 0..k-1
    for idx in range(dim):
        if idx & mask:  # if any of the first k bits is 1
            psi2[idx] = 0.0
    norm = np.linalg.norm(psi2)
    if norm < 1e-16:
        return psi2
    return psi2 / norm

def reduced_density(psi, L, N):
    """Return the reduced density matrix of the first L qubits."""
    dimA = 1 << L
    dimB = 1 << (N - L)
    psi_resh = psi.reshape(dimA, dimB)
    return psi_resh @ psi_resh.conj().T

# ---------- Simulation Parameters ----------------------------------
N = 12                     # total qubits
Ls = np.arange(1, N)       # block sizes from 1 to N-1

# Build the two "vacuum" states
psi_area = bell_pairs_state(N)  # area‑law state
psi_vol  = random_state(N)      # volume‑law (Haar random)

# Function to compute entropy curve S(L)
def entropy_curve(psi):
    return [entanglement_entropy(reduced_density(psi, L, N)) for L in Ls]

# Compute entanglement entropy for both vacua
S_area = entropy_curve(psi_area)
S_vol  = entropy_curve(psi_vol)

# Simulate "mass insertion" by projecting first k qubits, for k=1,2,3
results = []
for k in [1, 2, 3]:
    psi_mass = project_qubits_zero(psi_vol, k, N)
    S_mass   = entropy_curve(psi_mass)
    deltaS   = np.array(S_vol) - np.array(S_mass)
    results.append((k, S_mass, deltaS))

# ---------- Plot: Entropy Scaling ----------------------------------
plt.figure(figsize=(8,5))
plt.plot(Ls, S_area, 'o-', label='Area‑law (Bell pairs)')
plt.plot(Ls, S_vol,  'o-', label='Volume‑law (random state)')
for k, S_mass, _ in results:
    plt.plot(Ls, S_mass, 'o--', label=f'After projecting {k} qubit(s)')
plt.xlabel('Block size $L$')
plt.ylabel('Entanglement entropy $S(L)$ (bits)')
plt.title(f'Entanglement entropy vs block size (N={N} qubits)')
plt.legend(loc='upper left')
plt.grid(True)
plt.tight_layout()

# ---------- Plot: Displaced Entropy -------------------------------
plt.figure(figsize=(8,5))
for k, _, deltaS in results:
    plt.plot(Ls, deltaS, 'o-', label=f'ΔS after projecting {k} qubit(s)')
plt.xlabel('Block size $L$')
plt.ylabel('Entropy displacement ΔS(L) (bits)')
plt.title('Entropy displaced by "mass" insertion')
plt.legend(loc='upper left')
plt.grid(True)
plt.tight_layout()
plt.show()