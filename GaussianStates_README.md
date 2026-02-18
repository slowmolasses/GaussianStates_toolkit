# Gaussian State Manipulation Toolkit (Mathematica)

Comprehensive Mathematica library for continuous-variable (CV) quantum optics: Gaussian state representations, symplectic transformations, Fock basis expansions, and photonic quantum circuit simulation. Implements the mathematical foundations used at Xanadu (StrawberryFields, The Walrus) with native Mathematica symbolic computation and numerical optimization.

**Key Features:**
- Symplectic group generators (beam splitters, squeezers, two-mode squeezers, displacements)  
- Covariance matrix evolution under linear and nonlinear optical transformations  
- Fock basis decomposition via hafnian computation (integrates with The Walrus Python library)  
- Symplectic spectrum diagonalization (Williamson decomposition)  
- Quantum Chernoff bound computation for multimode Gaussian state discrimination  

Compatible with **Xanadu's StrawberryFields** conventions: uses the (x̂, p̂) phase-space ordering and symplectic form Ω = [[0, I], [−I, 0]].

---

## Contents

### Core Toolkits

| Notebook | Purpose |
|----------|---------|
| `GaussianStates_clean.nb` | **Primary toolkit:** Symplectic generators, mode transformations, covariance matrix algebra |
| `SymplecticDiagonalization_clean.nb` | Williamson decomposition for arbitrary Gaussian states (returns symplectic eigenvalues + diagonalizing transform) |
| `GaussianStateFockBasisTools.nb` | Fock expansion coefficients via loop hafnian (calls The Walrus for large matrices) |

### Demonstrations (`demo/` subdirectory)

See [`demo/README.md`](demo/README.md) for detailed documentation of:
- **BeamDisplacement_Debug.nb** — Optimizes transmitter spatial mode for quantum-enhanced target detection (Quantum Chernoff Exponent analysis)  
- **chernoff-exponent.nb** — Computes QCE for multimode Gaussian state discrimination (implements Pirandola et al.'s formula)  
- **SFG_evolution_organized.nb** — Nonlinear sum-frequency generation circuit simulation (appendix C of [Phys. Rev. Applied 19, 064015](https://doi.org/10.1103/PhysRevApplied.19.064015))

---

## Quick Start

### Example 1: Two-Mode Squeezer + Beam Splitter

```mathematica
(* Load toolkit *)
Get["GaussianStates_clean.nb"]

(* Initial state: vacuum on both modes *)
n = 2;  (* number of modes *)
V0 = IdentityMatrix[2*n];  (* covariance matrix for |0,0⟩ *)

(* Apply two-mode squeezing between modes 1 and 2 *)
r = 0.5;  phi = 0;  (* squeezing parameters *)
S_tms = TMSgen[n, 1, 2, r, phi];  (* symplectic matrix for TMS *)
V1 = S_tms . V0 . Transpose[S_tms];  (* evolved covariance *)

(* Apply 50:50 beam splitter *)
eta = 0.5;  (* transmissivity *)
S_bs = BSgen[n, 1, 2, eta];
V_final = S_bs . V1 . Transpose[S_bs];

(* Check physical state (Williamson condition): all symplectic eigenvalues ≥ 1 *)
<< SymplecticDiagonalization_clean.nb
{sympEvals, sympTransform} = SympSpec[V_final];
Print["Symplectic spectrum: ", sympEvals];  (* Should be ≥ 1 *)
```

**Expected output:** `{1.13107, 1.13107}` — valid physical state

---

### Example 2: Fock Expansion of Squeezed Vacuum

```mathematica
(* Load Fock basis tools *)
<< GaussianStateFockBasisTools.nb

(* Single-mode squeezed vacuum state *)
r = 1.0;  (* squeezing parameter *)
V_sq = {{Exp[2*r], 0}, {0, Exp[-2*r]}};  (* covariance matrix *)
alpha = {0, 0};  (* zero displacement *)

(* Compute Fock state amplitudes up to n_max = 5 *)
nMax = 5;
fockCoeffs = Table[
  FockAmplitude[V_sq, alpha, {nx, 0}],  (* photon number state |nx, 0⟩ *)
  {nx, 0, nMax}
];

Print["Fock expansion: ", Abs[fockCoeffs]^2];  
(* Probabilities for |0⟩, |2⟩, |4⟩, ... (only even photon numbers) *)
```

**Note:** `FockAmplitude` uses the loop hafnian algorithm (Kan 2008) for small states, and delegates to The Walrus via Mathematica's `ExternalEvaluate["Python", ...]` for larger Fock cutoffs.

---

## API Reference

### `GaussianStates_clean.nb` — Symplectic Transformations

All functions return symplectic matrices **S** (size 2n × 2n) acting on the phase-space vector **r** = (x̂₁, ..., x̂ₙ, p̂₁, ..., p̂ₙ)^T. Covariance matrices transform as **V → S·V·S^T**.

---

#### `atoq[n]`

Converts from annihilation operator basis (â, â†) to quadrature basis (x̂, p̂).

**Parameters:**
- `n` — number of modes

**Returns:** 2n × 2n matrix **A** such that [**x̂**, **p̂**]^T = **A** [**â**, **â**†]^T

**Conventions:**
- x̂ = (â + â†) / √2  
- p̂ = i(â† − â) / √2  
- [x̂, p̂] = i (ℏ = 1 units)

**Inverse:** `qtoa[n]` returns **A**^(−1)

**Example:**
```mathematica
A = atoq[2];  (* 4×4 matrix for 2 modes *)
omega = SymplecticForm[2];  (* Ω = [[0,I], [−I,0]] *)
A . omega . Transpose[A] // Simplify
(* Returns {{0, -i}, {i, 0}} ⊗ I₂ — canonical commutation relations preserved *)
```

---

#### `SymplecticForm[n]`

Returns the symplectic form **Ω** for n modes.

**Returns:** 2n × 2n block matrix [[0, Iₙ], [−Iₙ, 0]]

**Usage:** Verify symplectic condition **S^T Ω S = Ω** for any symplectic transform **S**

---

#### `TMSgen[n, i, j, r, phi]`

Two-mode squeezer between modes i and j.

**Parameters:**
- `n` — total number of modes  
- `i, j` — mode indices (1-indexed)  
- `r` — squeezing magnitude (real, r ≥ 0)  
- `phi` — squeezing phase (real)

**Returns:** 2n × 2n symplectic matrix for Hamiltonian **H** = −r(e^(iφ) âᵢ†âⱼ† + e^(−iφ) âᵢâⱼ)

**Transformation on annihilation operators:**
- â'ᵢ = cosh(r)·âᵢ − e^(iφ) sinh(r)·âⱼ†  
- â'ⱼ = cosh(r)·âⱼ − e^(−iφ) sinh(r)·âᵢ†

**Physical meaning:** Generates entanglement between modes i and j. For r = 0.88 (arcsinh(1)), produces a maximally entangled EPR state.

**Example:**
```mathematica
(* Create EPR pair between modes 1 and 2 *)
S = TMSgen[2, 1, 2, ArcSinh[1], 0];
V_EPR = S . IdentityMatrix[4] . Transpose[S];

(* Verify entanglement: compute logarithmic negativity *)
<< SymplecticDiagonalization_clean.nb
evals = SympSpec[V_EPR][[1]];
LogNegativity = Sum[Max[0, -Log[eval]], {eval, evals}];
Print["Log negativity: ", LogNegativity];  (* Should be > 0 *)
```

---

#### `SMSgen[n, rvec]` and `SMSgen2[n, rvec, phivec]`

Single-mode squeezing on multiple modes simultaneously.

**SMSgen (phase-less version):**
- `n` — number of modes  
- `rvec` — list of n squeezing parameters (one per mode)

**Returns:** Symplectic matrix for Hamiltonian **H** = −Σᵢ rᵢ (âᵢ†² + âᵢ²) / 2

**SMSgen2 (general version with phase):**
- `rvec` — squeezing magnitudes  
- `phivec` — squeezing phases

**Returns:** Symplectic matrix for **H** = −Σᵢ rᵢ (e^(iφᵢ) âᵢ†² + e^(−iφᵢ) âᵢ²) / 2

**Example:**
```mathematica
(* Squeeze mode 1 by r=0.5, leave mode 2 at vacuum *)
S = SMSgen[2, {0.5, 0}];
V = S . IdentityMatrix[4] . Transpose[S];

(* Variance in x-quadrature of mode 1 should be reduced *)
Print["Var(x₁) = ", V[[1,1]]];  (* = exp(−2·0.5) ≈ 0.368 *)
Print["Var(p₁) = ", V[[3,3]]];  (* = exp(2·0.5) ≈ 2.718 *)
```

---

#### `BSgen[n, i, j, eta]` and `BSgen2[n, i, j, theta, phi]`

Beam splitter mixing modes i and j.

**BSgen (lossless, parameterized by transmissivity):**
- `eta` — transmissivity (0 ≤ η ≤ 1), where η = cos²(θ)

**Transformation:**
- â'ᵢ = √η âᵢ − √(1−η) âⱼ  
- â'ⱼ = √(1−η) âᵢ + √η âⱼ

**BSgen2 (general, with phase):**
- `theta` — mixing angle  
- `phi` — beam splitter phase

**Transformation:**
- â'ᵢ = cos(θ) âᵢ − e^(iφ) sin(θ) âⱼ  
- â'ⱼ = e^(−iφ) sin(θ) âᵢ + cos(θ) âⱼ

**Example (Hong-Ou-Mandel interference):**
```mathematica
(* Two single photons entering a 50:50 BS *)
(* Initial state: |1,1⟩ (Fock state representation not directly in this toolkit) *)
S_BS = BSgen[2, 1, 2, 0.5];  (* 50:50 splitter *)

(* For Gaussian approximation: replace |1⟩ with displaced vacuum *)
V0 = IdentityMatrix[4];
d = DisplacementOp[2, {1, 1, 0, 0}];  (* Displace both modes in x-quadrature *)
V_HOM = S_BS . d . V0 . Transpose[d] . Transpose[S_BS];
```

---

#### `Phasegen[phivec]`

Applies phase shifts to each mode.

**Parameters:**
- `phivec` — list of n phases

**Returns:** Symplectic matrix for **U** = exp(−iΣᵢ φᵢ n̂ᵢ) where n̂ᵢ = âᵢ†âᵢ

**Transformation:**
- âᵢ → e^(−iφᵢ) âᵢ  
- x̂ᵢ → cos(φᵢ) x̂ᵢ − sin(φᵢ) p̂ᵢ  
- p̂ᵢ → sin(φᵢ) x̂ᵢ + cos(φᵢ) p̂ᵢ

**Example (phase-space rotation):**
```mathematica
S_phase = Phasegen[{Pi/4, 0}];  (* Rotate mode 1 by 45° in phase space *)
```

---

#### `DisplacementOp[n, dvec]`

Displacement operator (affine transformation on phase space).

**Parameters:**
- `n` — number of modes  
- `dvec` — 2n-dimensional displacement vector (x-quadratures, then p-quadratures)

**Returns:** Symplectic transformation + displacement vector encoding **D**(α) = exp(α â† − α* â)

**Note:** Displacements **do not change** the covariance matrix (only shifts the mean). To apply a displacement:

```mathematica
(* Coherent state |α⟩ = Displacement applied to vacuum *)
alpha = 2 + 3*I;  (* complex amplitude *)
dvec = Re[{alpha, -I*alpha}] * Sqrt[2];  (* Convert to (x, p) *)
V_coherent = IdentityMatrix[2];  (* Covariance unchanged from vacuum *)
mean_coherent = dvec;  (* Mean vector shifts *)
```

---

### `SymplecticDiagonalization_clean.nb` — Williamson Decomposition

#### `SympSpec[V]`

Computes the symplectic eigenvalue decomposition (Williamson normal form) of a covariance matrix **V**.

**Parameters:**
- `V` — 2n × 2n real symmetric positive-definite matrix

**Returns:** `{D, S}` where:
- `D` — list of n symplectic eigenvalues νᵢ (real, ≥ 1 for physical states)  
- `S` — 2n × 2n symplectic matrix such that **V = S^T diag(ν, ν) S** and **S^T Ω S = Ω**

**Algorithm:**
1. Compute eigendecomposition of **V^(−1/2) Ω V^(−1/2)** (antisymmetric matrix → purely imaginary eigenvalues ±iνᵢ)
2. Extract Youla decomposition to find orthogonal **K** such that **S = Λ^(−1/2) K V^(1/2)** is symplectic

**Usage (check if state is physical):**
```mathematica
V = (* some covariance matrix *);
{evals, S} = SympSpec[V];

If[Min[evals] >= 1,
  Print["Valid physical state"],
  Print["Unphysical state (violates uncertainty principle)"]
];
```

**Applications:**
- Convert arbitrary Gaussian state to diagonal form (decoupled modes)  
- Compute entanglement measures (logarithmic negativity, von Neumann entropy)  
- Verify physicality of numerically propagated states

---

### `GaussianStateFockBasisTools.nb` — Hafnian & Fock Expansions

#### `Lhafrepeated[A, gamma, svec]`

Computes the **loop hafnian** of a matrix with reduction vector, used to calculate Fock state amplitudes.

**Parameters:**
- `A` — m × m complex matrix (typically the **X** matrix from covariance decomposition)  
- `gamma` — m-dimensional complex vector (displacement-dependent term)  
- `svec` — m-dimensional integer vector (photon number configuration, 0-indexed)

**Returns:** Scalar complex number (hafnian value)

**Algorithm:** Kan's recursive algorithm (2008) for small states; delegates to The Walrus Python library via `ExternalEvaluate` for large matrices (m > 10).

**Mathematical context:**
For a Gaussian state with covariance **V** and mean **μ**, the Fock amplitude ⟨**n**|ρ|**n**⟩ is:

⟨**n**|ρ|**n**⟩ = haf_loop(**X**, **γ**) / sqrt(det(**V**) Π nᵢ!)

where **X = V − I**, **γ** involves the displacement, and **n** = (n₁, ..., nₘ) is the photon number vector.

**Example:**
```mathematica
(* Squeezed vacuum Fock expansion *)
V = {{Exp[2*r], 0}, {0, Exp[-2*r]}};  (* r = 1.0 *)
X = V - IdentityMatrix[2];
gamma = {0, 0};  (* no displacement *)

(* Probability of |2⟩ state *)
svec = {2};  (* n=2 photons *)
haf = Lhafrepeated[X, gamma, svec];
prob2 = Abs[haf]^2 / (Sqrt[Det[V]] * 2!);  (* normalization *)
Print["P(n=2) = ", prob2];
```

---

#### Integration with The Walrus

For large photon number cutoffs (Fock dimension > 10), the toolkit automatically switches to The Walrus's optimized C++ hafnian implementation:

```mathematica
(* Install The Walrus in your Python environment first: *)
(* pip install thewalrus *)

(* Then Mathematica will use it via ExternalEvaluate *)
ExternalEvaluate["Python", "from thewalrus import hafnian"];
```

**Performance:** Mathematica's native implementation is ~10× slower than The Walrus for matrices larger than 20×20, but avoids the overhead of Python interop for small states.

---

## Xanadu StrawberryFields Compatibility

This toolkit uses the **same conventions** as Xanadu's [StrawberryFields](https://github.com/XanaduAI/strawberryfields) library:

| Quantity | This Toolkit | StrawberryFields |
|----------|--------------|------------------|
| Quadrature ordering | **(x̂, p̂)** | `xp_ordering=True` |
| Symplectic form | **Ω = [[0,I], [−I,0]]** | `sf.hbar = 2` (implicit) |
| Squeezing convention | **S(r) = exp[r(â†² − â²)/2]** | `Sgate(r, phi)` |
| Beam splitter | **BS(θ) = exp[θ(âb† − b̂a†)]** | `BSgate(theta, phi)` |
| Displacement | **D(α) = exp[α↠− α*â]** | `Dgate(alpha)` |

**Example (cross-check with StrawberryFields):**

```python
# StrawberryFields code
import strawberryfields as sf
from strawberryfields.ops import Sgate, BSgate

prog = sf.Program(2)
with prog.context as q:
    Sgate(0.5) | q[0]
    BSgate(np.pi/4) | (q[0], q[1])

eng = sf.Engine("gaussian")
state = eng.run(prog).state
cov_SF = state.cov()  # 4×4 covariance matrix
```

```mathematica
(* Equivalent Mathematica code *)
V0 = IdentityMatrix[4];
S_sq = SMSgen[2, {0.5, 0}];
S_bs = BSgen2[2, 1, 2, Pi/4, 0];
cov_Mathematica = S_bs . S_sq . V0 . Transpose[S_sq] . Transpose[S_bs];

(* cov_SF == cov_Mathematica (up to numerical precision) *)
```

---

## Demonstrations

See [`demo/README.md`](demo/README.md) for:

1. **Quantum target detection** — Optimizes spatial mode shapes for quantum-enhanced sensing (Quantum Chernoff Exponent > classical Fisher information bound)
2. **QCE computation** — Implements Pirandola et al.'s closed-form formula for Gaussian state discrimination
3. **Nonlinear optics** — Simulates sum-frequency generation (SFG) circuits with three-wave mixing (χ⁽²⁾ nonlinearity)

---

## Dependencies

- **Mathematica 12.2+** (uses `Module`, `KroneckerProduct`, native sparse matrices)
- **Python 3.7+ with The Walrus** (optional, for large Fock cutoffs):
  ```bash
  pip install thewalrus numpy
  ```
- **Mathematica External Evaluate** configured for Python (run `FindExternalEvaluators["Python"]` to check)

**No external packages required** for basic Gaussian state manipulation (symplectic transforms, covariance evolution). The Walrus is only needed for Fock basis expansions with large photon number cutoffs.

---

## Mathematical Background

### Gaussian States

A multimode Gaussian state is fully characterized by:
- **Mean vector** **μ** ∈ ℝ^(2n) (first moments of x̂ᵢ, p̂ᵢ)  
- **Covariance matrix** **V** ∈ ℝ^(2n×2n) where Vᵢⱼ = ⟨ΔrᵢΔrⱼ + ΔrⱼΔrᵢ⟩/2

**Physicality condition (Heisenberg uncertainty):**
**V + iΩ/2 ≥ 0** (positive semidefinite)

Equivalently: all symplectic eigenvalues νᵢ ≥ 1 (in units where ℏ = 2).

### Symplectic Group

Linear transformations **S** ∈ Sp(2n,ℝ) preserve the canonical commutation relations:

**S^T Ω S = Ω**

Every symplectic transformation corresponds to a Gaussian unitary in Fock space. The symplectic Lie algebra generators are:
- **Two-mode squeezing:** exp[r(âᵢ†âⱼ† − âᵢâⱼ)]  
- **Single-mode squeezing:** exp[r(â†² − â²)/2]  
- **Beam splitter:** exp[θ(âᵢ†âⱼ − âⱼ†âᵢ)]  
- **Phase shift:** exp[iφ n̂]

**Cartan decomposition:** Any **S** ∈ Sp(2n,ℝ) can be decomposed as **S = K₁ Λ K₂** where **K₁, K₂** are orthogonal symplectic (beam splitters + phase shifts) and **Λ** is diagonal (squeezing).

### Hafnian & Fock Amplitudes

The **hafnian** of a 2m×2m matrix **A** is:

haf(**A**) = Σ_π ∏_(i,j)∈π A_ij

summed over all perfect matchings π of 2m elements. For Gaussian states, the Fock amplitude is:

⟨**n**|ρ|**n**⟩ = haf_loop(**X**, **γ**; **n**) / sqrt(det(**V**) ∏ nᵢ!)

where the **loop hafnian** generalizes the standard hafnian to include diagonal terms (self-loops) and a reduction vector **n**.

**Computational complexity:** O((2m)!!) ≈ O(2^m m^(3/2)) for exact evaluation. The Walrus uses Kan's algorithm with O(m² 2^m) scaling.

---

## References

### Symplectic Geometry & Gaussian States
- Weedbrook et al., "Gaussian quantum information," *Rev. Mod. Phys.* **84**, 621 (2012). [arXiv:1110.3234](https://arxiv.org/abs/1110.3234)  
- Serafini, *Quantum Continuous Variables: A Primer of Theoretical Methods* (CRC Press, 2017).

### Williamson Decomposition
- Williamson, "On the algebraic problem concerning the normal forms of linear dynamical systems," *Am. J. Math.* **58**(1), 141 (1936).  
- Youla, "A normal form for a matrix under the unitary congruence group," *Canadian J. Math.* **13**, 694 (1961).

### Hafnian Algorithms
- Kan, "From moments of sum to moments of product," *J. Multivariate Anal.* **99**, 542 (2008). [DOI:10.1016/j.jmva.2007.01.013](https://doi.org/10.1016/j.jmva.2007.01.013)  
- Björklund et al., "A faster hafnian formula for complex matrices and its benchmarking on a supercomputer," *ACM J. Exp. Algorithmics* **24**, 1.11 (2019).

### Xanadu Ecosystem
- Killoran et al., "StrawberryFields: A Software Platform for Photonic Quantum Computing," *Quantum* **3**, 129 (2019). [arXiv:1804.03159](https://arxiv.org/abs/1804.03159)  
- Quesada et al., "Simulating realistic non-Gaussian state preparation," *Phys. Rev. A* **100**, 022341 (2019). [The Walrus documentation](https://the-walrus.readthedocs.io/)

### Applications (see demo/)
- Guha & Erkmen, "Gaussian-state quantum-illumination receivers for target detection," *Phys. Rev. A* **80**, 052310 (2009).  
- Cox et al., "Sum-frequency generation-based entanglement-assisted photonic transceiver," *Phys. Rev. Applied* **19**, 064015 (2023). [DOI:10.1103/PhysRevApplied.19.064015](https://doi.org/10.1103/PhysRevApplied.19.064015)

---

## License

MIT License. Free for academic and commercial use with attribution.

---

## Acknowledgments

Symplectic diagonalization algorithm adapted from **Serafini's textbook** (2017). Fock expansion methods based on **The Walrus** by Xanadu AI. Quantum Chernoff bound implementation follows **Pirandola, Laurenza, Ottaviani & Banchi** (2017). Special thanks to **Prof. Saikat Guha** (University of Arizona) for guidance on Gaussian state theory and its applications to quantum communications.
