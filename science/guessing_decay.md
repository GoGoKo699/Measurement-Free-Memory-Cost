# An explicit decay guarantee for delayed-basis guessing

[Moment dynamics](moment_processing.md) · [Discrete ledger](control_cost.md) · [BB84 interface](guessing_processing.md)

The existing moment ledger gives an explicit exponential upper bound on whole-string guessing success below a necessary processing rate. This chapter inverts that ledger at finite width and time, then takes a controlled limit. It is a quantitative corollary of the established moment theorem, rather than an achievable processing capacity or a new uncertainty relation.

## 1. Premises and notation

Use the unchanged [memory model](model.md). The device receives only the $`k`$ BB84 challenge qubits; $`X`$ stays secret and only the basis string is disclosed after storage. The score $`p`$ is unconditional exact-string success, and aborts fail. All $`n`$ persistent physical qubits, including controller and workspace registers, undergo $`D`$ independent product-depolarizing rounds. Put $`q=1-\gamma\in(0,1)`$. Encoding and final decoding remain arbitrary free CPTP endpoints. The budget $`m_T`$ counts prescribed interior T/T-dagger slots in the fixed circuit. Mathematical noise flags remain unavailable to the controller.

For $`1\lt z\le2^{1/7}`$, write

```math
\lambda=\ln z,\qquad
\alpha=\alpha_z=\frac{\ln(1-q+qz^2)}{2\ln z},\qquad
M=2n.
\qquad\text{(D1)}
```

The common full-output logarithmic moment starts at $`M`$, contracts by at most $`\alpha`$ per backward noise layer, and increases by at most one per charged gate. Its direct all-cut floor is

```math
a_z=\left[2k+\frac{\ln p}{\lambda}\right]_+.
\qquad\text{(D2)}
```

For $`p=0`$, all upper bounds below hold trivially; use $`-\ln0=+\infty`$ in the asymptotic statement. Equation (D2) follows from the branch profile and score averaging in the moment chapter, before basis-dependent readout. It is not a floor on each individual noise branch.

## 2. Finite width and time

### Corollary D1 — Finite inverse of the injection ledger

Let $`n\ge k\ge1`$, $`D\ge1`$ and $`m_T\ge0`$. Define

```math
\begin{aligned}
U_z(D,n,m_T)
=\min\left\{M\alpha,\;
\min_{1\le\ell\le D}
\frac{\alpha(m_T+M\alpha^{\ell-1})}
{1+(D-\ell)(1-\alpha)}\right\}.
\end{aligned}
\qquad\text{(D3)}
```

Every permitted circuit obeys

```math
p\le\min\{1,z^{U_z-2k}\},
\qquad
-\ln p\ge(\ln z)[2k-U_z]_+.
\qquad\text{(D4)}
```

Take the strongest bound over admissible $`z`$. For $`D=0`$, set $`U_z=M`$; this gives the trivial upper bound one because $`n\ge k`$.

**Proof.** If $`a_z\gt0`$, the discrete injection ledger gives, for every integer $`1\le\ell\le D`$,

```math
m_T\ge
\frac{a_z}{\alpha}
+(D-\ell)\frac{(1-\alpha)a_z}{\alpha}
-M\alpha^{\ell-1}.
\qquad\text{(D5)}
```

Rearranging gives each inner upper bound in (D3). Immediately after traversing any noise layer backward, the moment is at most $`M\alpha`$, since every incoming envelope has rank at most $`M`$. The same all-cut floor therefore also satisfies $`a_z\le M\alpha`$. Consequently $`2k+\ln p/\lambda\le U_z`$. If $`a_z=0`$, that last inequality still holds because its left side is nonpositive and $`U_z\ge0`$. This proves (D4), including scores at or below chance. $`\square`$

The term $`M\alpha^{\ell-1}`$ retains the complete initial-width allowance. No assumption that the initial reserve is already of order $`k`$ is hidden in this finite statement. The minimum in (D3), without the separate width cutoff, is exactly the inverse of the necessary continuous scalar injection constraints. Its exactness is scalar; it supplies no physical circuit. The individual-gate logarithmic ledger and integer branch profiles may further strengthen it.

A convenient unoptimized form, valid for every $`1\le L\le D`$, is

```math
-\ln p\ge(\ln z)
\left[2k-
\frac{\alpha\bigl(m_T+2n\alpha^{L-1}\bigr)}
{1+(D-L)(1-\alpha)}\right]_+.
\qquad\text{(D6)}
```

This makes the tradeoff between coasting away the initial width and retaining a long charged interval explicit. The same proof works on an existing contiguous interval when its reserve bounds the logarithmic moment itself.

## 3. A uniform exponential consequence along growing families

### Corollary D2 — Explicit necessary-rate decay exponent

Fix $`q\in(0,1)`$ and a finite real $`\rho\ge0`$. Consider any sequence of permitted circuits, indexed by $`i`$, with integer parameters satisfying

```math
k_i\longrightarrow\infty,\quad D_i\longrightarrow\infty,\quad
n_i\ge k_i,\quad
\ln(n_i/k_i)=o(D_i),\quad
\limsup_i\frac{m_{T,i}}{k_iD_i}\le\rho.
\qquad\text{(D7)}
```

Then their unconditional scores satisfy

```math
\liminf_i\frac{-\ln p_i}{k_i}\ge E_q(\rho),
\qquad
E_q(\rho)=
\sup_{1\lt z\le2^{1/7}}
(\ln z)
\left[2-\rho\frac{\alpha_z}{1-\alpha_z}\right]_+.
\qquad\text{(D8)}
```

**Proof.** Fix one admissible $`z\gt1`$, independently of $`i`$. Put $`\kappa_z=-\ln\alpha_z\gt0`$ and choose

```math
L_i=1+\left\lceil\frac{\ln(2n_i/k_i)}{\kappa_z}\right\rceil.
\qquad\text{(D9)}
```

Then $`L_i=o(D_i)`$, so eventually $`1\le L_i\le D_i`$, and $`2n_i\alpha_z^{L_i-1}\le k_i`$. Equation (D3) implies

```math
\frac{U_{z,i}}{k_i}
\le
\frac{\alpha_z\bigl(m_{T,i}/(k_iD_i)+1/D_i\bigr)}
{(1-\alpha_z)(1-L_i/D_i)+1/D_i}.
```

The denominator tends to $`1-\alpha_z\gt0`$, giving $`\limsup_i U_{z,i}/k_i\le\rho\alpha_z/(1-\alpha_z)`$. Apply (D4), continuity of the positive part, and then take the supremum over fixed $`z`$. This does not exchange an optimized finite parameter with a limit: the same liminf is at least every one of the displayed fixed-parameter bounds. Subsequences with $`p_i=0`$ cause no problem under the extended logarithm convention. $`\square`$

Thus, whenever $`E_q(\rho)\gt0`$, for every $`\eta\gt0`$ the score eventually obeys $`p_i\le\exp[-k_i(E_q(\rho)-\eta)]`$. Both growth conditions in (D7) matter: the proof must have enough rounds to dissipate the initial width allowance. The noise parameter is fixed in this corollary; simultaneous weak-noise limits require additional uniform estimates.

## 4. Positivity, endpoints and certified optimization

For fixed $`q\in(0,1)`$, the function $`\alpha_z`$ strictly increases from $`q`$ as $`z`$ increases above one. Indeed, with $`t=2\ln z`$ and $`f(t)=\ln(1-q+qe^t)`$, strict convexity, $`f(0)=0`$, and $`f'(0)=q`$ make $`f(t)/t`$ strictly increasing. Consequently

```math
E_q(\rho)\gt0
\quad\Longleftrightarrow\quad
0\le\rho\lt\frac{2(1-q)}q.
\qquad\text{(D10)}
```

At the critical rate and above it, every positive-part term vanishes. Below it, sufficiently small fixed $`z-1\gt0`$ gives a positive term. The left endpoint is excluded but has continuous objective value zero. The right endpoint $`z=2^{1/7}`$ is allowed. In particular $`E_q(0)=2\ln2/7`$.

Optimization needs no unsupported numerical global-search claim. For $`\rho\gt0`$, let $`g(t)=t-f(t)\gt0`$. Then $`g''(t)=-f''(t)\lt0`$, and the objective before taking its positive part is

```math
F(t)=\left(1+\frac\rho2\right)t-
\frac\rho2\frac{t^2}{g(t)},\qquad
0\lt t\le\frac{2\ln2}{7}.
```

The identity

```math
\left(\frac{t^2}{g(t)}\right)''
=\frac{2(g-tg')^2}{g^3}-\frac{t^2g''}{g^2}\gt0
```

proves strict concavity of $`F`$. In the subcritical case its global maximum is positive and unique: it is the right endpoint if its derivative there is nonnegative, and otherwise the unique interior derivative zero. For interval certification one may use

```math
\frac{dF}{d\lambda}
=2-\rho\frac{t^2 f'(t)-f(t)^2}{(t-f(t))^2},
\qquad\lambda=t/2.
\qquad\text{(D11)}
```

The [verifier](../verification/verify_guessing_decay.py) encloses all logarithms by a rational series with a bounded tail, encloses the algebraic endpoint by exact seventh-power comparisons, and brackets the unique derivative zero by certified signs. For $`q=99/100`$ it obtains the following exponent enclosures (rounded outward):

| $`\rho`$ | Maximizer | Certified interval for $`E_q(\rho)`$ |
|---:|---|---:|
| $`0`$ | Right endpoint | $`[0.198042051588,0.198042051589]`$ |
| $`0.01`$ | Right endpoint | $`[0.089980725766,0.089980725767]`$ |
| $`0.015`$ | Right endpoint | $`[0.035950062854,0.035950062855]`$ |
| $`0.02`$ | Unique interior point | $`[0.000050419501,0.000050419502]`$ |
| $`2/99`$ | Zero objective | $`0`$ |

These certify the displayed converse objective, not the best possible operational exponent.

## 5. Stronger finite profiles and the recovery comparison

Let $`P_z`$ be any established nondecreasing concave majorant of the actual branch score as a function of $`Y=z^r`$, continued constantly at one. The moment chapter's full polygon and the concave majorants of its integer-algebra caps are examples. The same inverse-ledger argument gives

```math
p\le P_z(z^{U_z}).
\qquad\text{(D12)}
```

To see this, a score larger than the right side would force an all-cut logarithmic moment floor strictly greater than $`U_z`$, contradicting (D5). For scores at or below $`P_z(1)`$, the inequality already holds by monotonicity. Formula (D12) keeps the chance term and the low-rank shape discarded in (D2), and can therefore be substantially stronger.

For the existing certified polygon at $`k=12,n=120,q=99/100,z=27/25,D=2000`$, the new verifier obtains:

| Prescribed $`m_T`$ | Simple bound (D4) | Existing finite polygon with (D12) | Specified recovery profile |
|---:|---:|---:|---:|
| $`0`$ | $`p\le0.157699360778`$ | $`p\le0.000244168541`$ | $`p\le0.000244173658`$ |
| $`171`$ | $`p\le0.366182779858`$ | $`p\le0.248198005521`$ | $`p\le0.293651830277`$ |
| $`200`$ | $`p\le0.418535080865`$ | $`p\le0.310300082242`$ | $`p\le0.367329468646`$ |
| $`243`$ | $`p\le0.509249724669`$ | $`p\le0.417900524354`$ | $`p\le0.485618841419`$ |

These are outward upper bounds from the discrete injection relaxation. The existing individual-gate bound still excludes $`p=2/5`$ at $`m_T=243`$; that stronger fact is not lost merely because this particular relaxation gives a weaker number. Optimizing over all applicable ledgers and profiles remains legitimate.

For integer algebra refinements, maximize each proved cap over the allowed pairs $`(s,c)`$ with $`r=2s+c`$ and $`s+c\le n`$, then form its concave majorant in $`z^r`$. This preserves arbitrary permitted encoders and the actual common endpoint decoder. A finite upper profile is not assumed attainable.

The verifier explicitly covers all 169 such pairs for $`0\le r\le24`$. Besides the output-dimension bound it includes the existing dimension-lifting bound $`p\le2^s\beta^k`$, where $`\beta=(1+1/\sqrt2)/2`$. It compares the direct cap $`2^{(s-k)/6}`$ with the specified recovery cap

```math
\begin{aligned}
F_{r,s}&=\min\{1,2^{r-2k},2^{s-k}\},\\
u_{r,s}&=\min_{0\le j\le k}(\tau_j+A_jF_{r,s}),\\
p&\le\frac{1+\sqrt{(d-1)(du_{r,s}-1)}}{d},\qquad d=2^k.
\end{aligned}
```

Here $`\tau_j,A_j`$ are the normalized DFW coefficients in the [recovery comparison](guessing_processing.md#5-comparison-with-existing-fidelity-and-recovery-reductions). The last line inverts the Renes quadratic on $`[1/d,1]`$; lower actual scores are already covered. Rank zero uses its exact replacement-channel cap $`1/d`$. Dyadic root enclosures make all finite decisions rational. At this fixed $`z`$, the complete integer direct hull and the hull combining both families equal the existing three-vertex polygon. The specified recovery hull has vertex ranks $`0,13,16,19,24`$ and gives the last table column. Thus no further integer improvement appears in this complete fixed-parameter profile comparison; other parameters or other proved inequalities remain separate questions.

The specified normalized DFW/Renes recovery caps can be used in exactly the same manner, branchwise before averaging. Alternatively, use their established averaged-channel fidelity floor and the complementary recovery moment theorem, then invert the same scalar ledger. Neither route is automatically stronger for every parameter. In particular, the earlier $`223/52`$ ceilings concern a specified finite comparison of recovery-based converse formulas at $`p=2/5`$; they do not determine a full decay exponent or any physical upper gate cost. A recovery argument whose fixed spectral cutoff was chosen for a positive constant score cannot simply reuse that cutoff when $`p`$ itself decays exponentially in $`k`$.

The simple exponent (D8) remains a guaranteed quantitative consequence. This chapter does not claim that full polygons, integer refinements, or every recovery theorem have the same optimal asymptotic exponent. The zero-budget finite example already illustrates why discarding the low-rank shape can lose information.

## 6. Scope and verification

The proof inputs are inherited from the envelope dynamics, the finite BB84 branch profile and the exact scalar injection ledger. The addition is their explicit finite inversion and the quantified exponential-success corollary, with a certified one-variable optimization. The verifier checks 864 exact scalar inversion cases, all 169 permitted integer algebra pairs at the displayed finite parameters, rigorous finite score enclosures, and rational-interval optimizer certificates; the universal claims rely on the written proofs.

The critical value $`2(1-q)/q`$ is a necessary-rate threshold below which this converse proves exponential decay under (D7). It is not an achievable processing capacity, an optimal coefficient, or evidence that any circuit near the threshold exists. The score verifies none of the noise, controller or chronology premises. No protected records, adaptive free lookup, LWE preparation wrapper, hardware gate-time noise model or altered challenge distribution has been introduced.
