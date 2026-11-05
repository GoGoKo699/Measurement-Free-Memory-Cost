# Memory model and scientific results

## The question

A quantum memory must remove errors while retaining the unknown state and its entanglement with other systems. Supplying clean qubits removes entropy. Processing the syndrome determines which correction to apply. This project asks how much non-Clifford processing is necessary when that processing must occur inside the noisy device, without a protected external classical controller.

The [qualitative need for non-Clifford processing](../references/RELATED_WORK.md) in this setting is established prior work. The results here quantify its continuing cost, prove a sharp order for the lifetime achievable with Clifford processing alone, and give explicit constructions on both sides of the space-versus-processing comparison.

The model deliberately permits very powerful ideal operations between noise layers. This makes its converse apply independently of geometric locality, while placing the constructions in an ideal control model rather than a detailed hardware model.

## 1. Systems, noise, and allowed control

There are at most $`n`$ physical qubits at any instant during storage. This includes all data, syndrome records, temporary work qubits, controller registers, and any resource that retains information between noisy rounds. Variable active width is represented by a fixed pool of $`n`$ slots; a slot may be reset before reuse.

A storage circuit contains exactly $`D\ge0`$ noise layers, with $`D`$ an integer. Each layer applies the product channel

```math
\mathcal N_\gamma^{\otimes n},
\qquad
\mathcal N_\gamma(X)
=(1-\gamma)X+\gamma\,\frac I2\mathrm{Tr}X.
\qquad\text{(M1)}
```

The noise is independent across qubits and time. We write

```math
q=1-\gamma,\qquad
\kappa=-\ln q,
\qquad 0\lt \gamma\lt 1.
```

The symbol $`\ln`$ denotes the natural logarithm; $`\log_2`$ has base two.

Equivalently, the single-qubit Pauli probabilities are $`1-3\gamma/4`$ for identity and $`\gamma/4`$ each for $`X,Y,Z`$. A computational-basis bit flips with probability $`\gamma/2`$. The noise-strength convention matters when comparing thresholds or batching intervals.

The fixed circuit between noise layers may contain:

| Operation or resource | Treatment |
|---|---|
| Arbitrary global Clifford unitary | Allowed; zero non-Clifford charge |
| Reset of a selected subsystem to a fixed state | Allowed; zero non-Clifford charge |
| Physical T or T-dagger gate | One unit of T count |
| Other explicitly specified internal operations | Require a proved rank-growth charge or an actual Clifford+T implementation |
| Intermediate external classical computation or storage | Not an available resource |
| Initial encoder and final recovery | Arbitrary fixed CPTP endpoint channels; outside the charged storage interval |

A reset on $`A`$ is the replacement channel

```math
\mathcal R_{\tau,A}(X)
=\tau_A\otimes\mathrm{Tr}_A X,
\qquad\text{(M2)}
```

with tensor factors restored to their original order. The replacement state is independent of the unknown input and of the retained register. It may be mixed, nonstabilizer, or entangled within the simultaneously prepared subsystem. An external source that retains correlations across different preparation times is a persistent memory and is not supplied for free.

All operations between successive noise layers are ideal in this model. Their elementary gate depth does not create additional noisy rounds automatically. An arbitrary global Clifford may therefore contain a large syndrome-extraction circuit. Likewise, all elementary non-Clifford gates inside an ideal recovery operation are counted even though the operation can be completed before the next noise layer.

There is no postselection and no data-dependent choice of the externally prescribed schedule. Coherent control using counted registers is allowed and must be implemented using the permitted gates. The proofs also allow unrecorded Pauli measurements as zero-charge channels, but no conclusion about arbitrary adaptive measurement-based controllers follows from that observation.

## 2. Transmission criterion and endpoints

Let $`K=2^k`$, with $`k\ge1`$. An arbitrary CPTP encoder $`\mathcal E`$ maps the logical input to the memory, $`\Phi_D`$ is its storage channel, and an arbitrary fixed CPTP decoder $`\mathcal D`$ returns $`k`$ logical qubits:

```math
\Lambda=\mathcal D\circ\Phi_D\circ\mathcal E.
\qquad\text{(M3)}
```

The reference system is untouched. No protected quantum side channel connects encoder and decoder.

Use the normalized Bell state

```math
|\Omega_K\rangle
=K^{-1/2}\sum_{j=1}^{K}|j\rangle|j\rangle
```

and its overlap entanglement fidelity

```math
F_{\mathrm e}(\Lambda)
=\langle\Omega_K|
(\mathrm{id}\otimes\Lambda)
(|\Omega_K\rangle\langle\Omega_K|)
|\Omega_K\rangle.
\qquad\text{(M4)}
```

This is the squared-overlap convention. A replacement logical channel has fidelity exactly $`K^{-2}=4^{-k}`$. The converse assumes only $`F_{\mathrm e}\ge1-\epsilon`$. The constructions prove the stronger worst-case, reference-stable guarantee

```math
\delta_\diamond(\Lambda)
:=\frac12\|\Lambda-\mathrm{id}_K\|_\diamond
\le\epsilon,
\qquad\text{(M5)}
```

which implies the fidelity condition.

Endpoint recovery can perform nonlinear syndrome decoding and use work registers after storage has finished. Such endpoint work supplies no noiseless information during the $`D`$ noisy rounds. If a nonlinear recovery is instead applied before storage ends, its processing belongs to the charged circuit. The [decoder-specific rank refinement](control_cost.md#7-the-final-decoder-and-the-meaning-of-the-terminal-rank) makes this boundary quantitative.

## 3. Principal results

### Reference-stable forgetting with Clifford operations and resets

For a storage circuit using only Clifford operations and replacement resets, and for any fixed physical input state $`\omega_0`$, set $`\tau=\Phi_D(\omega_0)`$. Then

```math
\frac12\|\Phi_D-\mathcal R_\tau\|_\diamond
\le\min\{1,2nq^D\}.
\qquad\text{(M6)}
```

Thus all input dependence, including correlations with an arbitrary reference, is lost to accuracy $`\delta`$ after

```math
D\ge \frac{\ln(2n/\delta)}{-\ln(1-\gamma)}.
\qquad\text{(M7)}
```

The proof in [Binary Pauli envelopes](pauli_envelopes.md) tracks the dimension of one binary linear envelope containing every pulled-back Pauli support. Resets project its coordinates; one noise layer contracts its expected rank by a factor $`q`$. A zero envelope identifies an exactly input-independent sampled channel.

### An extensive budget for non-Clifford processing

Suppose $`F_{\mathrm e}\ge1-\epsilon`$ and $`0\le\epsilon\lt 1-4^{-k}`$. Define

```math
a=2k\,\frac{1-\epsilon-4^{-k}}{1-4^{-k}}\gt 0.
\qquad\text{(M8)}
```

Every Clifford+T memory obeys the simple bound

```math
m_T\ge
a\left[D\kappa-\ln\!\left(\frac{2n}{a}\right)\right]_+,
\qquad [x]_+=\max\{x,0\}.
\qquad\text{(M9)}
```

The count includes both $`T`$ and $`T^\dagger`$. At fixed error bounded away from the trivial threshold and weak noise, this is an extensive cost of order $`k\gamma D`$ once $`D\kappa`$ exceeds the width allowance $`\ln(2n/a)`$ by a fixed factor greater than one. This is a necessary gate budget; the constructions below give sufficient resources in the same model.

[Control cost](control_cost.md) gives sharper integer-gate and finite-depth inequalities, an exact optimization of the corresponding scalar rank constraints, and a bound for every contiguous storage interval. It also replaces $`2n`$ by a smaller decoder-specific terminal rank when justified. A fixed Pauli-axis rotation has rank-growth charge at most one, while a general operation on $`s`$ qubits has a safe charge at most $`2s`$; those structural charges are not automatically physical T counts.

The improved fidelity floor (M8) follows from an elementary channel fact. A logical branch whose superoperator rank is at most $`2^r`$ has

```math
F_{\mathrm e}\le\min\{1,2^{r-2k}\}.
\qquad\text{(M10)}
```

Indeed, the superoperator trace equals $`K^2 F_{\mathrm e}`$, all eigenvalues of a CPTP map have modulus at most one, and the number of its nonzero eigenvalues, counted with algebraic multiplicity, is at most its rank. Applying the convex secant bound to each branch before averaging forces expected envelope rank at least $`a`$ at every cut. Every T gate can add at most one binary generator.

### Reliability-dependent continuing processing

For $`D\ge1`$ and $`0\lt\epsilon\lt1-4^{-k}`$, set $`\delta=\epsilon/(1-4^{-k})`$ and $`b=1-\delta`$. The [rank-extinction theorem](reliability_processing.md) strengthens the necessary processing budget to

```math
m_T\ge\left[
\gamma Db\max\left\{2k,
\frac{[\ln(Db/\delta)]_+}{\ln(1/\gamma)}\right\}-2n
\right]_+.
\qquad\text{(M14)}
```

The first term in the maximum retains the payload floor. The second follows from the conditional extinction probability $`\gamma^r`$ of a rank-$`r`$ envelope and the requirement that extinction occur with probability at most $`\delta`$. The proof uses an absorbing envelope update, so a zero envelope stays zero through every permitted primitive. It sums disjoint first-extinction events, applies Jensen's inequality to the surviving cut weights, and converts cumulative rank into T count using the same noise/gate balance. No noise flag is made available to the controller.

At finite width and $`D\ge1`$, exact zero Bell infidelity is impossible under positive depolarization. At $`D=0`$, there is no reliability-dependent interior processing requirement. The existing finite and interval bounds remain complementary; (M14) is not asserted to dominate each of them numerically.

### Sharp order of the zero-T quantum lifetime

For fixed $`k`$ and target fidelity $`4^{-k}\lt F_0\lt 1`$, the largest lifetime achievable with Clifford operations and resets alone satisfies

```math
D_{\max}(n)=\Theta(\gamma^{-1}\ln n)
\qquad\text{(M11)}
```

for sufficiently large width and $`0\lt \gamma\le1/16`$, with constants depending on the fixed target parameters.

The constructive direction in [Clifford quantum lifetime](clifford_lifetime.md) stores unknown quantum states, rather than only a classical bit. A stabilizer code produces a history of syndrome bits; many noisy copies of each bit are propagated through ternary trees. All majority voting and nonlinear decoding occur at the final endpoint. Every tree node remains a counted noisy qubit. Its exponential width and ideal global fanout are essential resources.

### Explicit five-qubit coherent recovery

For $`0\lt \gamma\le1/20`$ and $`0\lt \epsilon\lt 1/2`$, put

```math
\ell=1+\ln\frac{k(\gamma D+1)}{\epsilon},
\qquad \alpha=\log_2 5.
```

The [coherent recovery construction](coherent_recovery.md) achieves (M5) with

```math
n=O(k\ell^\alpha),
\qquad
m_T=O\bigl(k(\gamma D+1)\ell^\alpha\bigr).
\qquad\text{(M12)}
```

It concatenates the five-qubit code, batches weak noise, and applies an explicit syndrome correction using at most 42 T/T-dagger gates per decoding node. No additional work qubits are needed beyond the code block. The constants and finite integer parameter choices are given in the construction.

This five-qubit upper bound retains a polylogarithmic reliability overhead and its explicit range $`\gamma\le1/20`$.

### Constant-iteration Tanner recovery and matching order

Fix a good quantum Tanner-code family satisfying the hypotheses stated in [Tanner memory](tanner_memory.md). For $`0\lt\gamma\le1/20`$, $`k,D\ge1`$ and $`0\lt\epsilon\lt1/2`$, a fixed Clifford+T/reset memory achieves (M5) with

```math
n=O\!\left(k+\ln\frac{J}{\epsilon}\right),
\qquad m_T=O(nJ),
\qquad J=O(\gamma D+1).
\qquad\text{(M13)}
```

Here $`n`$ includes the code block and all syndrome, controller and temporary registers. The local decoder is a fixed Boolean circuit compiled coherently; every prescribed T/T-dagger slot is counted. Scratch is reset on a fixed schedule before use, so current syndrome extraction is ideal under the existing interior-processing assumption. No syndrome record is protected between noisy rounds. A fixed number of inner five-qubit concatenation levels reduces each outer position's fresh error to the outer decoder's tolerance; this level depends only on the fixed family, not on $`k,D,\epsilon`$. The bare Tanner option remains valid at its smaller direct noise range.

The code, constant-iteration decoder and repeated-round residual invariant are inherited from Gu et al. Their parallel runtime alone would not prove the charged gate bound. The added compilation counts local minimizations, comparisons, deterministic ties, updates and fixed routing. The correction channel is reference-stable, and correlated X/Z components of depolarizing noise are combined by a union bound.

Fix $`\gamma\in(0,1/20]`$, $`C\gt0`$ and $`0\lt\epsilon_0\lt1/2`$. Put $`L=\ln(D/\epsilon)`$. For every family with $`D\to\infty`$, $`k\ge1`$, $`0\lt\epsilon\le\epsilon_0`$ and $`n\le C(k+L)`$, (M14) implies $`m_T=\Omega_{\gamma,C,\epsilon_0}(D(k+L))`$, uniformly over $`k,\epsilon`$. Indeed $`b`$ is bounded below and $`\ln(Db/\delta)=L+O_{\epsilon_0}(1)`$; the maximum in (M14) controls $`k+L`$, and sufficiently large $`D`$ absorbs the $`2n`$ reserve. Equation (M13) attains this processing order at width $`O(k+L)`$. This closes the reliability-dominated case, including fixed $`k`$, as well as the payload-dominated case. Constants may be enormous and depend on the fixed positive noise strength; no optimal coefficient or uniform weak-noise law is inferred. The extended constructive range is proved by fixed-level inner protection, rather than inferred from an LDPC threshold.

## 4. How the results fit together

The conclusions distinguish two ways to support storage:

| Choice of internal processing | Established resource behavior | Main physical cost |
|---|---|---|
| Clifford operations and resets only | Optimal lifetime order is logarithmic in width at fixed logical size and fidelity | An exponential number of noisy syndrome records |
| Five-qubit coherent nonlinear recovery | Explicit width and gate count with a polylogarithmic reliability overhead | Continuing T/T-dagger processing |
| Constant-iteration Tanner recovery | Width $`O(k+\ln(J/\epsilon))`$ and work $`O(nJ)`$ for $`\gamma\le1/20`$ with fixed-level inner protection | Continuing charged decoding with counted resettable scratch |
| Any high-fidelity protocol in the model | Must satisfy the gate-budget and interval lower bounds | Information cannot cross every temporal cut without sufficient internal processing |

Fresh magic states alone do not expand the backward binary envelope. Their availability therefore does not eliminate the processing lower bound. This is a statement about operations inside the specified memory, rather than a magic-state consumption theorem for all forms of error correction.

All constructions permit ideal global processing. A local hardware design with faulty gates, measurement latency, or correlated noise must account for those additional features before importing an achievability conclusion.

### Delayed-basis BB84 scores as an operational premise

The [guessing-to-processing transfer](guessing_processing.md) applies the same ledger to a device that receives only $`k`$ BB84 challenge qubits. The secret string $`X`$ is never supplied; only its independent basis string is revealed after charged storage ends. For unconditional exact-string score $`p`$ and $`k\ge8`$, one full-output envelope with terminal rank $`2n`$ satisfies, at every primitive cut,

```math
R_t\ge a_k(p)=2k\left[\frac{p-2^{-k}}{1-2^{-k}}\right]_+.
```

Whenever $`a_k(p)\gt0`$, substitute this floor for $`a`$ in the gate and every-contiguous-interval bounds; the transfer also gives a positive all-size floor for $`1\le k\le7`$. The proof factors each sampled CPTP tail through a finite quantum-plus-classical Pauli algebra before any basis-dependent readout. All persistent controller records remain noisy and counted, and mathematical noise flags are not physically supplied. Established branchwise recovery reductions, with spectral normalization retained, reproduce both displayed floor formulas for every challenge size. The score does not verify the device-model premises or automatically import a cryptographic preparation wrapper.

### Exponential moments retain the noise law

The [moment-processing theorem](moment_processing.md) strengthens this ledger without changing the device model. For $`z\gt 1`$, define $`X_z=\log_z\mathbb E[z^r]`$. A product-noise layer gives $`X_z'\le\alpha_zX_z`$, with

```math
\alpha_z=\frac{\ln(1-q+qz^2)}{2\ln z}.
```

Each T/T-dagger gate adds at most one, and free Clifford operations and input-independent resets do not increase this potential. Direct quantum-plus-classical score bounds supply a floor for $`1\lt z\le2^{1/7}`$; a fixed recovery of the averaged full channel supplies a complementary floor (using $`F_0=1-\epsilon`$ for the Bell-fidelity premise) for $`1\lt z\le2`$. Once a positive floor $`a_z`$ is proved at every cut, replace the scalar ledger's $`a,q`$ by $`a_z,\alpha_z`$. The initial full-output potential is $`2n`$. An interval entry bound must control $`X_z`$ itself, rather than only mean rank.

For fixed positive guessing score and fixed $`q`$, along $`n\ge k`$, $`k,D\to\infty`$ and $`\ln(n/k)=o(D)`$, the discrete ledger implies

```math
m_T\ge(2-o(1))\,kD\,\frac{1-q}{q}.
```

This improves the fixed-score coefficient obtained from the convenient mean-rank floor; the weak-noise order remains $`k\gamma D`$. Finite examples and exact certificates are in the moment chapter. The conclusion is necessary and conditional on the same noise, controller and gate assumptions, with no physical-achievability claim.

The [quantitative decay corollary](guessing_decay.md) inverts this necessary ledger while retaining the finite initial-width reserve. Under the same growing-family conditions, a budget rate $`\limsup m_T/(kD)\le\rho`$ forces $`\liminf(-\ln p)/k\ge E_q(\rho)`$. The full generating-function version strengthens $`E_q`$ to $`\mathcal E_q`$ by propagating all moment arguments. It has $`\mathcal E_q(0)=\ln2`$ and is positive exactly when $`\rho\lt2(1-q)/q`$. This is a guarantee of exponential score decay below a necessary rate; it is not an achievable processing capacity.

For the stronger exponent, write

```math
\phi_j(\lambda)=\frac12\ln[1+q^j(e^{2\lambda}-1)],
\quad S_q(\lambda)=\sum_{j=1}^{\infty}\phi_j(\lambda),
\quad A(\lambda)=\min\left\{2\lambda,\frac{\ln2+5\lambda}{6},\ln2\right\}.
```

Then the same limiting hypotheses give

```math
\liminf\frac{-\ln p}{k}\ge\mathcal E_q(\rho),
\qquad\mathcal E_q(\rho)=
\max_{0\le\lambda\le\ln2}[A(\lambda)-\rho S_q(\lambda)]_+.
```

The [full proof and finite-width version](guessing_decay.md) retain the initial width and count every prescribed gate slot. The improvement comes from the noise-profile composition, without improving the branchwise guessing theorem.

## 5. Reading and verification routes

Read [Binary Pauli envelopes](pauli_envelopes.md) for the forgetting proof, then [Control cost](control_cost.md) for the fidelity bottleneck and gate budgets and [Reliability processing](reliability_processing.md) for the absorbing-extinction argument and all-payload fixed-noise matching law. The constructive chapters [Clifford quantum lifetime](clifford_lifetime.md), [five-qubit coherent recovery](coherent_recovery.md), and [Tanner memory](tanner_memory.md) establish the different ways to pay for storage. [Interpretation](interpretation.md) states the matching-order regime and the remaining questions about coefficients, finite constants and physical implementation.

[Related work](../references/RELATED_WORK.md) identifies inherited ingredients, the closest published theorem, and the distinctions from measurement-based or externally decoded models. [Verification](../verification/README.md) describes the finite checks and their scope. Finite computations test the explicit algebra and circuits; universal claims rely on the proofs.


