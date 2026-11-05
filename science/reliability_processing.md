# Reliability requires continuing non-Clifford processing

[Memory model](model.md) · [Pauli envelopes](pauli_envelopes.md) · [Payload and interval laws](control_cost.md) · [Tanner construction](tanner_memory.md)

A small logical register still needs increasingly reliable protection when it is stored for a long time. The payload bound alone does not capture that cost. This chapter proves a second lower bound by controlling the probability that the backward Pauli envelope disappears completely. Combined with the payload bound, it matches the Tanner construction's processing order at fixed physical noise for every payload scaling, under a stated width constraint.

The complete-erasure survival mechanism and its use with Jensen's inequality are established in the cumulative-spacetime analysis of [Bharti, Haug and Tanggara, arXiv:2608.26272v1](https://arxiv.org/abs/2608.26272), Theorem 1 and Supplementary Information E. Their cost counts physical storage locations. The adaptation here uses the binary rank of a backward Pauli envelope, then converts its cumulative rank into prescribed T/T-dagger slots. Clifford operations and input-independent resets remain free in this model. The two resource counts are different.

## 1. Model, absorbing propagation and the finite statement

Use exactly the [existing memory model](model.md): $`n`$ noisy physical slots, $`D`$ independent product-depolarizing rounds with $`0<\gamma<1`$, a fixed Clifford+T/T-dagger/reset circuit, arbitrary free CPTP endpoints, and no protected controller or cross-time side channel. Every persistent data, syndrome, resource and work register is noisy and counted. The count $`m_T`$ includes every prescribed interior T/T-dagger slot. The identity/replacement flags below are a mathematical decomposition, not information supplied to the device.

Propagate one simultaneous Pauli envelope backward from the full physical output, initially $`W_n=\mathbb F_2^{2n}`$. For a T rotation about Pauli label $`w`$, choose

```math
V'=\begin{cases}
V,&[v,w]=0\text{ for every }v\in V,\\
V+\langle w\rangle,&\text{otherwise}.
\end{cases}
\qquad\text{(R1)}
```

The [rotation support formula](control_cost.md#3-rank-growth-charges-for-gates) proves this update valid and its rank growth at most one. In particular, $`V=0`$ stays zero. Clifford symplectic images, reset coordinate projections, sampled depolarizing intersections and unrecorded Pauli measurements also preserve zero and have the existing rank bounds. Blindly adjoining $`w`$ when $`V=0`$ would create an artificial nonzero envelope and must not be used here.

Write $`[x]_+=\max\{x,0\}`$ and, for $`0<\delta<1`$ and integer $`s\ge1`$, define

```math
\mathcal R_\gamma(s,\delta)=
\frac{\gamma s(1-\delta)}{\ln(1/\gamma)}
\left[\ln\frac{s(1-\delta)}{\delta}\right]_+.
\qquad\text{(R2)}
```

### Theorem R1 — Finite reliability budget

Suppose $`D\ge1`$, $`k\ge1`$, recovered Bell entanglement fidelity satisfies $`F_{\mathrm e}\ge1-\epsilon`$, and $`0<\epsilon<1-4^{-k}`$. Put

```math
\delta=\frac{\epsilon}{1-4^{-k}},\qquad b=1-\delta.
```

Every permitted memory obeys

```math
m_T\ge\left\lceil
\left[\mathcal R_\gamma(D,\delta)-2n\right]_+
\right\rceil.
\qquad\text{(R3)}
```

Combining the same additive ledger with the established all-cut payload floor gives the stronger joint statement

```math
m_T\ge\left\lceil\left[
\gamma Db\max\left\{2k,
\frac{[\ln(Db/\delta)]_+}{\ln(1/\gamma)}\right\}-2n
\right]_+\right\rceil.
\qquad\text{(R4)}
```

These are necessary counts with a finite initial-width reserve. Keep the other valid logarithmic, individual-gate and moment inequalities as well; (R3) need not dominate them at every finite parameter. A half-diamond error guarantee at most $`\epsilon`$ implies the fidelity premise. At $`D=0`$ there is no reliability charge from this theorem. Exact zero error with $`D\ge1`$ and finite width is impossible, as shown below.

## 2. Rank-dependent extinction probability

### Lemma R2 — A replacement layer kills rank $`r`$ with probability at least $`\gamma^r`$

Let $`V\le W_n`$ have positive binary rank $`r`$. Form an $`r\times2n`$ generator matrix with independent rows and choose $`r`$ independent coordinate columns. Let $`J`$ be the physical qubits containing those columns. Then $`|J|\le r`$, and restriction from $`V`$ to all coordinates of $`J`$ is injective.

If the next backward noise layer replaces every qubit in $`J`$ by the maximally mixed state, the [sampled noise update](pauli_envelopes.md#2-the-simultaneous-backward-envelope) gives

```math
V'=V\cap\ker\pi_E=\{0\},\qquad J\subseteq E.
```

Additional replacements cannot restore a surviving label. Conditional on the currently exposed suffix, the next layer's flags are independent Bernoulli variables. The probability of this sufficient event is $`\gamma^{|J|}\ge\gamma^r`$. The pivot set may be chosen as a deterministic function of the current envelope; it is not a physical adaptive operation or supplied flag. This proves the lemma.

This is a rank-dependent bound, rather than the event that all $`n`$ physical qubits are replaced. It applies to an arbitrary envelope, including those produced by non-Clifford gates and arbitrary fixed endpoint decoding. Rank one spanned by a weight-one Pauli gives the sharp single-layer value $`\gamma`$.

## 3. Survival bounds the cumulative rank

Fix a complete identity/replacement realization for the storage circuit, and traverse it backward using the absorbing updates. Immediately before backward noise layer $`j`$, let $`r_j`$ be the envelope rank and $`S_j=\mathbf 1_{\{r_j>0\}}`$. The current envelope is determined by suffix flags; the layer being crossed has not yet been sampled. Put

```math
A=\sum_{j=1}^D\Pr[S_j=1],\qquad
B=\sum_{j=1}^D\mathbb E r_j.
\qquad\text{(R5)}
```

A zero envelope means every pulled-back output observable is scalar. The sampled full physical storage channel is therefore a replacement channel. Composing it with any encoder and decoder still gives a logical replacement channel, whose Bell fidelity is exactly $`4^{-k}`$. Since branch probabilities are input-independent and Bell fidelity is linear in the channel, the probability $`p_0`$ of reaching zero by the end of the backward traversal satisfies

```math
1-\epsilon\le F_{\mathrm e}
\le1-(1-4^{-k})p_0,\qquad p_0\le\delta.
\qquad\text{(R6)}
```

Absorption implies $`\Pr[S_j=1]\ge1-p_0\ge b`$, and therefore $`A\ge Db`$. First-extinction events at distinct noise layers are disjoint. Extinction caused by an intervening free reset is an additional possible cause of $`p_0`$, not a negative contribution. Lemma R2 consequently gives

```math
\sum_{j=1}^D\mathbb E[S_j\gamma^{r_j}]
\le\Pr[\text{first extinction occurs at a noise layer}]
\le p_0\le\delta.
\qquad\text{(R7)}
```

Normalize the joint measure on layer indices and surviving suffix histories by $`A>0`$. Its mean rank is $`B/A`$, since rank is zero on the omitted histories. Jensen's inequality for the convex function $`x\mapsto e^{-x\ln(1/\gamma)}`$ yields

```math
\delta\ge A\exp\left[-\ln(1/\gamma)\frac BA\right].
```

Thus

```math
B\ge\frac{A}{\ln(1/\gamma)}
[\ln(A/\delta)]_+
\ge\frac{Db}{\ln(1/\gamma)}[\ln(Db/\delta)]_+.
\qquad\text{(R8)}
```

The second step is valid for all $`0<\delta<1`$: $`x[\ln(x/\delta)]_+`$ is zero on $`[0,\delta]`$ and has positive derivative on $`(\delta,\infty)`$. No independence between successive ranks or survival conditioning on future flags has been assumed. Current-layer independence is used only in Lemma R2 before that layer is exposed.

## 4. Convert cumulative rank into charged processing

At each noise layer, the [existing mean-rank contraction](pauli_envelopes.md#2-the-simultaneous-backward-envelope) gives expected rank loss at least $`\gamma\mathbb E r_j`$. Free primitives do not increase rank; each prescribed T/T-dagger slot increases it by at most one. Telescope over every primitive, including gaps next to the unrestricted endpoints. Starting from rank $`2n`$ gives

```math
\mathbb E r_{\rm end}\le2n+m_T-\gamma B,
\qquad m_T\ge\gamma B-2n.
\qquad\text{(R9)}
```

Combining (R8), nonnegativity and integrality proves (R3). The established Bell-fidelity rank floor is $`\mathbb E r_j\ge2kb`$ at every cut, including before each noise layer. Hence $`B\ge2kbD`$ as well. Taking the maximum of these two lower bounds on the same $`B`$ proves (R4). This does not identify rank with stored logical qubits or assert that the rank inequalities are a realizable optimal circuit.

If $`\epsilon=0`$, (R6) requires $`p_0=0`$. But the first backward noise layer completely replaces all $`n`$ slots with positive probability $`\gamma^n`$, forcing a zero envelope and a replacement logical channel. This contradiction proves the zero-error assertion for finite width. Extra ideal processing cannot restore a state that has become input-independent.

### Corollary R3 — Every contiguous interval

For an interval containing $`D_I\ge1`$ noise layers and $`m_{T,I}`$ prescribed T/T-dagger slots, let $`M_I`$ be a proved upper bound on the expected envelope rank entering the interval in backward order. Then

```math
m_{T,I}\ge\left\lceil
[\mathcal R_\gamma(D_I,\delta)-M_I]_+
\right\rceil.
\qquad\text{(R10)}
```

Always $`M_I=2n`$ is valid. Use the actual propagated global suffix envelope, which may depend on later noise, rather than restarting a newly chosen independent envelope. The global fidelity premise still bounds its zero probability at every cut by $`\delta`$. The first-extinction events within the interval have total probability at most $`\delta`$, so (R5)–(R9) apply with $`D_I,M_I,m_{T,I}`$. Boundary gates are assigned to the chosen primitive cuts, as in the existing interval law. At $`D_I=0`$, only the trivial nonnegative count follows.

## 5. Matching processing order for every payload scaling

### Corollary R4 — Fixed-noise matching under a width constraint

Fix $`0<\gamma\le1/20`$, $`0<\epsilon_0<1/2`$ and a constant $`C>0`$. Let $`D\to\infty`$, allow arbitrary $`k\ge1`$ and $`0<\epsilon\le\epsilon_0`$ to vary with $`D`$, and impose

```math
L=\ln(D/\epsilon),\qquad n\le C(k+L).
\qquad\text{(R11)}
```

There are constants $`c>0`$ and $`D_0`$, depending only on $`\gamma,C,\epsilon_0`$, such that every allowed memory with recovered Bell fidelity at least $`1-\epsilon`$ satisfies

```math
m_T\ge cD[k+\ln(D/\epsilon)]\qquad(D\ge D_0).
\qquad\text{(R12)}
```

**Proof.** For every $`k\ge1`$, $`\epsilon\le\delta\le4\epsilon/3`$ and $`b\ge b_0:=1-4\epsilon_0/3>0`$. Thus

```math
\ln(Db/\delta)\ge L+\ln(3b_0/4).
```

For sufficiently large $`D`$, uniformly over $`k,\epsilon`$ in the stated domain, the last expression is at least $`L/2`$. In (R4), therefore,

```math
\max\left\{2k,\frac{[\ln(Db/\delta)]_+}{\ln(1/\gamma)}\right\}
\ge\frac12\left(2k+\frac{L}{2\ln(1/\gamma)}\right)
\ge c_1(k+L)
```

for a fixed $`c_1>0`$. The expression before the outer positive part in (R4) is at least $`[\gamma b_0c_1D-2C]\,(k+L)`$. Choose $`D_0`$ so the bracket is at least $`\gamma b_0c_1D/2`$. This proves (R12).

The [Tanner construction](tanner_memory.md) achieves half-diamond error at most $`\epsilon`$, width $`O(k+L)`$ and $`m_T=O_\gamma(D[k+L])`$ in this same fixed-noise range. Its width constant depends on the chosen admissible family. Thus for a sufficiently large constant width factor permitting that construction, the necessary and sufficient processing order is

```math
m_T^{\rm opt}=\Theta_\gamma(D[k+\ln(D/\epsilon)]).
\qquad\text{(R13)}
```

Constants also depend on the fixed width allowance and code family. This describes optimal order within the stated width class, not an unrestricted optimization over arbitrarily large memories. In particular, fixed $`k`$ now has $`\Theta_\gamma(D\ln(D/\epsilon))`$ order; the earlier large-payload result is recovered when $`k`$ dominates the reliability logarithm. No additional $`\ln(n/k)=o(D)`$ condition is needed for this corollary.

## 6. Delayed-basis BB84 scores

Use one full-output envelope before any basis-dependent final readout. If the unconditional whole-string guessing score is at least $`p`$, with $`2^{-k}<p<1`$, then every zero-envelope branch has input-independent physical output. Its success after basis disclosure is at most $`2^{-k}`$; aborts can only lower it. Branch probabilities are independent of the hidden input and basis. Therefore

```math
p\le1-(1-2^{-k})p_0,\qquad
p_0\le\delta_p:=\frac{1-p}{1-2^{-k}}.
\qquad\text{(R14)}
```

Replace $`\delta`$ by $`\delta_p`$ in (R3) and (R10). This uses a separate score argument and does not identify guessing probability with Bell fidelity. For fixed $`k,\gamma,p>2^{-k}`$ and $`n=O(\ln D)`$, it gives $`m_T=\Omega_\gamma(D\ln D)`$ as $`D\to\infty`$. The conclusion may hold below the classical-storage score: noisy classical records also need continuing processing in the stated device model. It is not a certificate of quantum storage or a verification of the noise/controller premises. Perfect score $`p=1`$ is impossible at finite width with $`D\ge1`$ for the same complete-replacement reason.

The ideal finite-dimension BB84 guessing theorem is unchanged. Its new consequence is this reliability-sensitive conditional processing witness; no LWE preparation/processing compiler is imported.

## 7. Certified finite examples

For $`k=1`$, $`n=200`$, $`\gamma=\epsilon=1/100`$, exact outward logarithm enclosures give:

| Noise rounds $`D`$ | Reliability budget (R3) | Earlier individual-gate bound | Earlier scalar injection ledger |
|---:|---:|---:|---:|
| 100,000 | 2,989 | 2,439 | 1,983 |
| 1,000,000 | 38,422 | 24,503 | 19,923 |

Every entry is a rounded necessary T/T-dagger count under the stated Bell-fidelity premise. The last two columns use the named [control-cost formulas](control_cost.md), rather than every possible moment refinement. These parameters are not claimed achievable. The [certificate program](../verification/certify_reliability_examples.py) uses rational atanh-series remainders to enclose each logarithm and prove the integer rounding. The two fixed-width examples are not an assertion that a fixed finite memory can survive arbitrarily long.

## 8. Verification, attribution and limits

The [standard-library verifier](../verification/verify_reliability.py) checks small Pauli subspaces, exact replacement-mask hazards and mean contraction, absorbing gate support, full fixed schedules and first-extinction accounting. Exact logarithmic certificates provide finite examples. These checks supplement the universal proofs above and do not establish physical feasibility of the numerical examples.

The stochastic Pauli/envelope framework, CPTP replacement bottleneck, payload rank floor and T-growth ledger are inherited from this repository and its cited antecedents. The first-extinction/Jensen reliability mechanism is closely related to the complete-erasure spacetime argument of Bharti, Haug and Tanggara. The addition is its rank-dependent form under free Clifford/reset processing, the conversion to a charged T budget, the interval/BB84 transfer and its combination with the existing Tanner construction. The literature comparison is scoped; it is not an exhaustive priority certificate.

The finite $`-2n`$ reserve is essential. Exponentially large zero-T archives are not excluded by a result restricted to $`n=O(k+\ln(D/\epsilon))`$. The explicit coefficient contains $`1/\ln(1/\gamma)`$, so fixed-noise matching does not establish a uniform weak-noise optimum. Optimal numerical T coefficients, unrestricted width/work tradeoffs, noisy elementary processing gates, geometric locality and an achievable processing capacity remain open.
