# Processing bounds from the noisy envelope distribution

[Model](model.md) · [Mean-rank ledger](control_cost.md) · [BB84 interface](guessing_processing.md)

The actual product-noise distribution supplies a stronger processing bound than reducing every cut to its mean envelope rank. The quantity used here is a logarithmic exponential moment. It contracts under the same noise, grows by at most one across a T gate, and admits an all-cut floor from either delayed-basis BB84 guessing or Bell recovery fidelity. The resulting gate and interval laws retain the original device model.

In the growing-block, long-duration regime of Corollary M3, a fixed positive whole-string score and fixed noise strength give a leading necessary charge approaching $`2k(1-q)D/q`$. The convenient mean-rank score floor instead gives the coefficient $`2pk(1-q)D/q`$. Both have order $`k\gamma D`$ at weak noise. This result improves the dependence on a nonvanishing score; it does not determine optimal physical gate cost.

## 1. Model, chronology and one full-output envelope

The device receives only the uniformly random $`k`$-qubit BB84 challenge. The verifier keeps $`X`$ secret and discloses only $`\Theta`$ after all charged storage ends. Let $`p=\Pr[\widehat X=X]`$ be the unconditional whole-string score, with aborts counted as failures.

Use the [fixed-circuit model](model.md): all $`n`$ persistent data, controller, syndrome, resource and workspace qubits undergo $`D`$ independent product-depolarizing rounds. Write $`q=1-\gamma\in(0,1)`$. Fixed Clifford operations and input-independent replacement resets are free; actual interior T/T-dagger gates are charged. The budget counts the prescribed circuit slots; it is not an average charge obtained by randomly choosing between different storage circuits. Encoding and basis-dependent final decoding are arbitrary CPTP endpoints. No protected persistent record or uncharged nonlinear controller is supplied.

Initialize a common backward envelope at the **complete physical output before basis-dependent readout**, with rank $`2n`$. At primitive cut $`t`$, sample the independent future identity/replacement choices $`b`$ and write

```math
r_{b,t}=\dim V_{b,t},\qquad
M_{z,t}=\mathbb E_b z^{r_{b,t}},\qquad
x_{z,t}=\log_z M_{z,t},\quad z\gt 1.
\qquad\text{(M1)}
```

The flags describe channels in a proof decomposition. They are not given to the device. The [branchwise CPTP factorization and common-decoder averaging](guessing_processing.md#2-one-envelope-and-a-quantum-plus-classical-cptp-bottleneck) remain the interface to the operational score.

## 2. Product noise contracts the logarithmic moment

### Lemma M1 — Conditional exponential-moment bound

Let $`V\leq\mathbb F_2^{2n}`$ have rank $`r`$. A noise layer independently retains each site with probability $`q`$ and produces the backward envelope $`V'=V\cap W_S`$, where $`S`$ is the retained-site set. For every real $`z\geq1`$,

```math
\mathbb E\!\left[z^{\dim V'}\mid V\right]
\leq (1-q+qz^2)^{r/2}.
\qquad\text{(M2)}
```

**Proof.** Choose $`r`$ coordinate bits on which restriction from $`V`$ is injective. Let $`d_j\in\{0,1,2\}`$ count the selected bits on physical site $`j`$, so $`\sum_jd_j=r`$. Every vector of $`V\cap W_S`$ has zero selected bits outside $`S`$. The same restriction remains injective, giving

```math
\dim(V\cap W_S)\leq\sum_{j\in S}d_j.
```

Independence bounds the conditional moment by $`\prod_j(1-q+qz^{d_j})`$. Put $`A=1-q+qz^2`$ and $`B=1-q+qz`$. The identity

```math
A-B^2=q(1-q)(z-1)^2\geq0
```

shows that each singleton pivot contributes at most $`A^{1/2}`$ and each pair contributes $`A`$. The product is at most $`A^{r/2}`$. $`\square`$

Define

```math
\alpha_z=\frac{\ln(1-q+qz^2)}{2\ln z},
\qquad q\lt \alpha_z\lt 1.
\qquad\text{(M3)}
```

The next backward noise layer is independent of the incoming suffix envelope. Lemma M1 and concavity of $`u^{\alpha_z}`$ therefore give

```math
M_z'\leq\mathbb E[(z^r)^{\alpha_z}]
\leq M_z^{\alpha_z},\qquad
x_z'\leq\alpha_z x_z.
\qquad\text{(M4)}
```

A free Clifford preserves rank and a reset projection cannot increase it. A primitive with proved rank-growth charge $`b_g`$ satisfies $`M_z'\leq z^{b_g}M_z`$, hence $`x_z'\leq x_z+b_g`$. In particular, a T/T-dagger gate adds at most one. The terminal value is $`x_z=2n`$, and $`0\leq x_z\leq2n`$ at every cut.

These statements use the actual independent noise probabilities. They require no concentration hypothesis. Surviving rank is not generally stochastically dominated by the paired-coordinate example; Lemma M1 is an exponential-moment inequality.

## 3. The direct BB84 profile supplies an all-cut floor

The [quantum-plus-classical bottleneck](guessing_processing.md#3-from-branch-scores-to-an-all-cut-rank-floor), using the companion's [finite BB84 theorem](https://github.com/GoGoKo699/Lattice-Memory-Certification/blob/84d33890bc2bd8312325fc616524963e9c025db6/science/guessing-bound.md), gives

```math
p_{b,t}\leq h_k(r_{b,t}),\qquad
h_k(r)=
\begin{cases}
2^{r-k},&0\leq r\leq5k/6,\\
2^{(r-2k)/7},&5k/6\leq r\leq2k,\\
1,&r\geq2k.
\end{cases}
\qquad\text{(M5)}
```

Here $`p_{b,t}`$ is the actual common decoder's score; $`p=\mathbb E_b p_{b,t}`$ at every primitive cut. The upper profile need not be attainable.

Fix $`1\lt z\leq2^{1/7}`$ and set

```math
b_0=2^{-k},\quad H=2^{-k/6},\quad
Y_1=z^{5k/6},\quad Y_2=z^{2k}.
```

As functions of $`Y=z^r`$, both exponential pieces of (M5) are convex. Let $`P_z`$ be the least concave polygon through or above

```math
(1,b_0),\qquad (Y_1,H),\qquad(Y_2,1),
\qquad\text{(M6)}
```

continued constantly at one. Keep the middle vertex exactly when

```math
H\geq b_0+(1-b_0)\frac{Y_1-1}{Y_2-1};
```

otherwise use the endpoint chord. Each exponential piece lies below its own endpoint secant, so $`P_z(z^r)\geq h_k(r)`$. Concavity and actual-score averaging imply

```math
p\leq\mathbb E_bP_z(z^{r_{b,t}})
\leq P_z(M_{z,t}).
\qquad\text{(M7)}
```

For $`p\gt b_0`$, invert this polygon and put $`a_{z,\mathrm{dir}}=\log_z L_z(p)`$, where

```math
L_z(p)=1+(p-b_0)\frac{Y_2-1}{1-b_0}
```

if the middle vertex is dropped. If it is retained, use

```math
L_z(p)=
\begin{cases}
1+(p-b_0)\dfrac{Y_1-1}{H-b_0},&b_0\lt p\leq H,\\
Y_1+(p-H)\dfrac{Y_2-Y_1}{1-H},&H\leq p\leq1.
\end{cases}
\qquad\text{(M8)}
```

Set $`a_{z,\mathrm{dir}}=0`$ at or below chance. Equations (M7)–(M8) prove $`x_{z,t}\geq a_{z,\mathrm{dir}}`$ at **every** primitive cut. The floor is positive at every above-chance score, including scores below one half. As $`z\downarrow1`$, it recovers the earlier rank polygon or chord, while $`\alpha_z\to q`$.

For a simple asymptotic bound, (M5) also obeys

```math
h_k(r)\leq\min\{1,2^{(r-2k)/7}\}.
\qquad\text{(M9)}
```

In the first region this is equivalent to $`r\leq5k/6`$; in the second it is equality. If $`1\lt z\leq2^{1/7}`$, then $`p\leq z^{-2k}M_z`$, giving the simpler floor $`x_z\geq2k+\log_zp`$ when positive. At $`z=2^{1/7}`$ and $`p\geq H`$, the full polygon gives exactly $`2k+7\log_2p`$.

Integer algebra sizes, the width constraint and other proved branch caps can improve (M6). Any concave majorant of those branch caps in $`Y=z^r`$ may replace $`P_z`$; Section 6 gives an exact finite certificate.

## 4. A complementary recovery floor

Let $`\Phi`$ be the averaged challenge-to-full-output channel. The [normalized DFW/Renes reduction](guessing_processing.md#5-comparison-with-existing-fidelity-and-recovery-reductions) supplies a lower bound $`F_0\geq4^{-k}`$ on its optimum squared Bell recovery fidelity. With $`d=2^k`$ and $`p\geq1/d`$, one may take

```math
F_0=\max\left\{4^{-k},\frac p d,
\max_{0\leq\ell\leq k}
\frac{[g_d(p)-\tau_\ell]_+}{A_\ell}\right\},
\qquad\text{(M10)}
```

where $`g_d(p)=p^2+(1-p)^2/(d-1)`$, and $`\tau_\ell,A_\ell`$ are the coefficients defined in that reduction.

Choose one optimal CPTP recovery of the averaged channel. At every cut append that **same** recovery to each sampled tail. It remains independent of the mathematical flags. The full-output envelope covers it, and the superoperator-rank lemma bounds the squared branch fidelity by $`\min\{1,2^{r-2k}\}`$. Fidelity is linear in the recovered channel, so

```math
F_0\leq\mathbb E_b\min\{1,2^{r_{b,t}-2k}\}.
```

For $`1\lt z\leq2`$, the exponential is convex in $`Y=z^r`$. Its endpoint secant yields the all-cut floor

```math
a_{z,\mathrm{rec}}=
\log_z\!\left[1+
\frac{F_0-4^{-k}}{1-4^{-k}}(z^{2k}-1)\right].
\qquad\text{(M11)}
```

This is positive for every above-chance score because measure-and-prepare already gives $`F_0\geq p/d\gt 4^{-k}`$. It may be optimized over a larger parameter range than the direct polygon. One can also apply the recovery bounds branch by branch before constructing the moment polygon; averaging the recovery premise first can lose strength. No claim of universal dominance between these alternatives is needed.

## 5. Processing and contiguous intervals

### Theorem M2 — Moment processing ledger

Fix $`z\gt 1`$ and any positive proved all-cut floor $`a_z`$, from Sections 3–4 or a certified branch majorant. Put $`M=2n`$ and $`\kappa_z=-\ln\alpha_z`$. The [existing scalar processing proofs](control_cost.md#4-logarithmic-and-additive-resource-laws) apply to $`x_z`$ with retention $`\alpha_z`$, giving

```math
B\geq a_z[D\kappa_z-\ln(M/a_z)]_+,
```

and, for actual Clifford+T implementations,

```math
m_T\geq
\left\lceil\frac{[D\kappa_z-\ln(M/a_z)]_+}
{\ln(1+1/a_z)}\right\rceil.
\qquad\text{(M12)}
```

The discrete injection bound is

```math
B\geq\mathcal L(D,M,a_z,\alpha_z),
\qquad\text{(M13)}
```

with the same $`\mathcal L`$ as in [Theorem C6](control_cost.md#theorem-c6--discrete-processing-budget). Round (M13) upward when $`B=m_T`$. A noise layer also requires $`a_z\leq2n\alpha_z`$. These are necessary conditions; scalar optimality does not give a physical construction.

For a contiguous interval, use its own number of noise rounds and an upper bound on its entering **logarithmic moment** $`x_z`$. The universal reserve $`2n`$ is valid. A deterministic rank cap or a direct moment bound may improve it. An upper bound on $`\mathbb E r`$ alone cannot replace this reserve, since $`\log_z\mathbb E z^r\geq\mathbb E r`$. The same global floor $`a_z`$ holds inside every interval, and independence of its noise from the suffix proves the same contraction.

Take the strongest bound over proved floors, admissible $`z`$, the two ledgers and the earlier mean-rank bounds. The limiting direct family includes the earlier score ledger; a particular finite parameter choice need not dominate it.

## 6. A certified improvement at the original parameters

For $`k=12,p=2/5,n=120,q=99/100`$, take $`z=27/25`$. A finite majorant uses the output-dimension and direct quantum-dimension bounds together with the established purely classical BB84 cap. The [verifier](../verification/verify_moment_processing.py) proves, using exact rational arithmetic, that all 169 integer pairs $`(s,c)`$ with $`2s+c\leq24`$ lie below the concave polygon through

```math
(1,1/4096),\qquad(z^{11},1/4),\qquad(z^{24},1).
\qquad\text{(M14)}
```

The classical input to this certificate is $`p\leq\beta^k`$ for a fully classical memory, where $`\beta=(1+1/\sqrt2)/2`$. Here is a direct proof under the present input promise. For a deterministic decoder $`f`$, let $`S_f=\sum_\theta\rho_{f(\theta),\theta}`$. The Gram matrix of the selected BB84 vectors has absolute row sums at most $`(1+1/\sqrt2)^k`$: basis strings at Hamming distance $`h`$ contribute overlap at most $`2^{-h/2}`$. Thus $`\|S_f\|_\infty\leq(1+1/\sqrt2)^k`$. A classical encoder is a POVM $`\{E_c\}`$ with deterministic postprocessing $`f_c`$, and therefore

```math
p=4^{-k}\sum_c\mathrm{Tr}(E_cS_{f_c})
\leq4^{-k}(1+1/\sqrt2)^k\mathrm{Tr}I
=\beta^k.
```

This established cap allows an unrestricted classical alphabet, so it remains valid with counted finite classical memory. The finite certificate needs 84 output-dimension witnesses, 84 direct quantum-dimension witnesses and one classical witness; no recovery inequality is needed for this particular polygon.

Above rank 24, the constant continuation at one suffices. Since the score is $`2/5`$, inversion gives

```math
M_z\geq\frac{4z^{11}+z^{24}}5,
\qquad a_z=14.84082963327\ldots,
\qquad\alpha_z=0.990724980097\ldots.
```

Certified logarithm enclosures then give:

| Charged storage | Earlier mean-rank count | Certified moment count |
|---|---:|---:|
| 2,000 noise rounds | 171 | 244 |
| 700-round interval, entering moment reserve 240 | 39 | 58 |

The new unrounded individual-gate expressions lie in $`(243.117317956157,243.117317956158)`$ and $`(57.347537146205,57.347537146206)`$, respectively. These are valid certified lower counts, not globally optimized parameters or achievable physical costs.

The earlier recovery comparison remains correct for the convenient **mean-rank** floors. The present theorem uses actual product-noise moments instead. The specified normalized recovery family can also be fed into that stronger ledger **branchwise before averaging**. The verifier includes the squared-recovery caps, all normalized spectral thresholds, the Renes conversion, measure-and-prepare and output-dimension caps, and the established classical bound. The additional dimension-lifting cap $`p\leq2^s\beta^k`$ is inactive at the certificate's rank-16 witness and does not remove it. This cap follows by tracing out the quantum register to obtain a classical encoder and bounding quantum guessing in each classical block by its dimension times the largest posterior. A rational mixture witness, 67 exact parameter intervals and an analytic tail certify that, for these same parameters and all $`z\gt 1`$, its individual-gate and injection formulas yield rounded necessary counts at most 223 and 52. The direct finite polygon yields 244 and 58. The numbers 223 and 52 are ceilings on the strength of those converse formulas, not physical upper gate counts. This comparison specifies its bound family and does not assert dominance over every possible recovery theorem. At scores at or below one half, preserving this order of operations is particularly important.

## 7. Long storage at a fixed positive score

### Corollary M3 — Leading charge without a fixed-score prefactor

Let $`k\to\infty`$, $`n\geq k`$, $`q\in(0,1)`$ be fixed, and suppose $`D\to\infty`$ with $`\ln(n/k)=o(D)`$. If $`p\geq p_0\gt 0`$, every feasible family obeys

```math
B\geq(2-o(1))k\frac{1-q}{q}D,
\qquad m_T\geq(2-o(1))k\frac{1-q}{q}D
\quad\text{for Clifford+T circuits}.
\qquad\text{(M15)}
```

**Proof.** Choose $`\lambda=\ln z=1/\sqrt k`$, eventually below $`\ln2/7`$. Equation (M9) gives

```math
a_z\geq2k+\frac{\ln p_0}{\lambda}
=2k-O_{p_0}(\sqrt k).
```

For $`f(t)=\ln(1-q+qe^t)`$ on $`0\leq t\leq2\lambda`$, one has $`f'(0)=q`$ and $`f''(t)\leq q(1-q)e^{2\lambda}`$. Taylor's remainder therefore implies

```math
q\leq\alpha_z\leq q[1+(1-q)e^{2\lambda}\lambda],
\qquad\alpha_z\longrightarrow q.
```

The scalar coasting time $`t_0=\lfloor\ln(2n/a_z)/(-\ln\alpha_z)\rfloor`$ is $`o(D)`$. The exact injection relaxation is at least

```math
(D-t_0-1)a_z\frac{1-\alpha_z}{\alpha_z},
```

which proves (M15). $`\square`$

The earlier displayed mean-rank score floor has leading coefficient $`2p_0k(1-q)/q`$. Equation (M15) removes that fixed-success prefactor. The weak-noise order remains $`\Theta(k\gamma D)`$. The result concerns feasible width/time families and does not guarantee indefinite storage at fixed width.

Equivalently, under these conditions, if $`\limsup m_T/(kD)\lt 2(1-q)/q`$, the unconditional whole-string score must tend to zero. This is a vanishing-score consequence; the argument proves neither an exponential decay rate nor a matching capacity or circuit construction.

## 8. Proof dependencies and limits

The Pauli envelope, fixed primitive charges, CPTP factorization and scalar ledgers are established in this repository. The finite BB84 quantum-dimension theorem belongs to Lattice Memory Certification; DFW/Renes supply the complementary recovery bounds. The additional step here is the conditional moment contraction and its composition with the score profiles.

The earlier optimality of the rank-only chord does not obstruct this theorem: it applied after replacing the noise distribution by arbitrary mixtures constrained only by $`\mathbb E r`$. Here the actual product-noise transition controls $`\mathbb E z^r`$ before its logarithm is taken.

All persistent-controller, noise independence and endpoint assumptions remain premises. A score verifies none of them. Mathematical flags stay unavailable to the device, and no LWE preparation or resource-preserving compilation is imported. The [verification guide](../verification/README.md) distinguishes exact certificates, finite checks and written all-parameter proofs.
