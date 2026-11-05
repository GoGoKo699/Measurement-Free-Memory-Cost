# Interpretation: width, storage time, and internal processing

The [model](model.md) separates three resources: noisy physical width, fresh-state preparation, and non-Clifford operations performed while the unknown state is stored. The main conclusion is that an unlimited supply of replacement states does not remove the need for continuing non-Clifford processing at modest width.

The [forgetting theorem](pauli_envelopes.md), [processing budget](control_cost.md), and the [syndrome archive](clifford_lifetime.md), [five-qubit recovery](coherent_recovery.md) and [Tanner recovery](tanner_memory.md) quantify this statement in one operational model.

## 1. Why the binary rank is the useful quantity

A Clifford maps a Pauli observable to another Pauli. A reset maps its binary label by deleting coordinates, multiplying the observable by a scalar that depends on the prepared state. These actions remain linear on labels even when a reset prepares a magic state.

Some coefficients vanish, and several labels can acquire the same image. The actual set of surviving labels can therefore be irregular. The proof encloses it in a binary linear space instead of assuming it remains a subspace. This enclosure is an upper bound on dependence: extra labels can weaken a converse, but cannot make it invalid.

Independent depolarization removes enough coordinates to contract the expected envelope dimension by $`q=1-\gamma`$ in every noisy round. This simultaneously controls all observables, rather than following a single expectation value or a separate path for each Pauli. If the envelope reaches dimension zero, the whole channel has become a replacement channel.

With internal non-Clifford gates, the envelope can grow. A T gate can add only one binary generator. A channel that preserves a logical Bell state needs a sufficiently large envelope at every temporal cut. The processing budget follows from combining these statements.

Binary rank, operator-space dimension, and logical qubits are different quantities:

| Quantity | Meaning | Use in the proof |
|---|---|---|
| $`r=\dim_{\mathbb F_2}V`$ | Number of independent binary Pauli labels in the enclosure | Contracts under noise and has a small gate-growth charge |
| $`2^r=\dim_{\mathbb C}\mathcal A_V`$ | Dimension of the span of all Paulis in that binary space | Bounds the rank of a logical channel passing through the cut |
| Quantum block size $`2^s`$ | Largest coherent matrix factor in the Pauli algebra | Exposes the limitation of a purely commuting, classical algebra |
| $`k`$ | Number of unknown logical qubits to preserve | Sets the fidelity-dependent floor on expected rank |

In particular, a large commuting algebra can retain much classical information while carrying no entanglement. The superoperator-rank argument gives the main scalar floor; the complementary Schmidt-number argument retains this quantum-versus-classical distinction.

## 2. What the initial width allowance buys

For the general decoder-independent bound, use the full terminal envelope $`M=2n`$. At fidelity $`F_{\mathrm e}\ge1-\epsilon`$ with $`0\le\epsilon\lt1-4^{-k}`$, define

```math
a=2k\,\frac{1-\epsilon-4^{-k}}{1-4^{-k}}.
```

Set $`\kappa=-\ln q`$. The simple necessary budget is

```math
m_T\ge
a\left[D\kappa-\ln(2n/a)\right]_+.
```

The first term grows with noisy storage duration. Physical width contributes the finite logarithmic allowance $`\ln(2n/a)`$. At fixed high fidelity, if $`D\kappa`$ exceeds that allowance by a fixed factor greater than one, the necessary T count is $`\Omega(k\gamma D)`$. At weak noise $`\kappa=\gamma+O(\gamma^2)`$. For finite encoding rate $`k/n=\Theta(1)`$, the allowance corresponds to only order $`1/\gamma`$ noisy rounds.

The exact finite-depth law in [Control cost](control_cost.md) improves the constants and gives a precise scalar optimization. Its optimum can first consume the initial rank reserve, then restore the minimum necessary rank before each subsequent noise layer. That is an optimization of necessary rank constraints. It does not construct a memory realizing those ranks.

The same argument applies to every contiguous interval. A protocol cannot satisfy an arbitrarily long interior interval solely by moving its charged operations before or after that interval.

The [delayed-basis BB84 application](guessing_processing.md) replaces the fidelity premise by an unconditional whole-string guessing score. The input is only the quantum challenge, its secret label is never disclosed, and its basis is disclosed after storage. For $`k\ge8`$, score $`p\gt2^{-k}`$ forces $`a_k(p)=2k(p-2^{-k})/(1-2^{-k})`$ at every cut, so the same gate and interval inequalities apply. A positive score excess can reflect retained classical information and does not alone demonstrate quantum storage. Existing branchwise fidelity/recovery reductions, with spectral normalization retained, reproduce the convenient mean-rank score floors, including the small-block polygon. The direct floor explains the operational connection with its explicit controller and chronology assumptions. The [moment theorem](moment_processing.md) gives a further consequence in the same model: actual product-noise probabilities contract the logarithm of an exponential rank moment, while each charged gate has the same additive growth bound. Combining this dynamics with a fixed recovery of the averaged channel strengthens the gate and interval ledger. For fixed positive guessing score and fixed $`q`$, along growing challenges and storage durations with $`\ln(n/k)=o(D)`$, its necessary rate approaches $`2k(1-q)/q`$ per round. This improves the fixed-score coefficient of the convenient mean-rank law while preserving the order $`k\gamma D`$ at weak noise. The entering interval reserve must bound the logarithmic moment; a mean-rank reserve cannot automatically be reused. The [inverse ledger](guessing_decay.md) additionally bounds the score itself at finite width and duration. In the growing-family limit it gives an explicit exponent $`E_q(\rho)`$, positive exactly for prescribed T rate $`\rho\lt2(1-q)/q`$. Propagating the full generating-function argument through successive noise layers strengthens that exponent to $`\mathcal E_q`$, with $`\mathcal E_q(0)=\ln2`$ and the same positivity threshold. This states how quickly success must vanish; it does not identify an achievable memory capacity.

### Reliability also constrains the time spent at each rank

The [reliability converse](reliability_processing.md) uses an absorbing choice of the same envelope. At binary rank $`r\gt0`$ and symplectic rank $`2s`$, at most $`r-s`$ physical pivot qubits identify every label. Replacing all of them annihilates the envelope, with conditional probability at least $`\gamma^{r-s}`$. Current replacement choices are independent of the later suffix defining the incoming envelope. A zero envelope stays zero under every adjoint primitive; in particular, a T update adds its axis only when some existing label anticommutes with it. Discarding the symplectic information gives the simpler rank-only consequence below.

Let $`\delta=\epsilon/(1-4^{-k})\in(0,1)`$ and $`b=1-\delta`$. Bell fidelity at least $`1-\epsilon`$ permits eventual extinction with probability at most $`\delta`$. First-extinction events at distinct noise layers are disjoint, so their total hazard is bounded by that same probability. Jensen's inequality, weighted by the surviving cuts, forces cumulative incoming rank at least

```math
\sum_{j=1}^D\mathbb E r_j
\ge\frac{Db}{\ln(1/\gamma)}[\ln(Db/\delta)]_+.
```

Noise loses at least $`\gamma\mathbb E r_j`$ expected rank at layer $`j`$, whereas a T/T-dagger gate adds at most one. The initial reserve is at most $`2n`$. Consequently

```math
m_T\ge\left[
\frac{\gamma Db}{\ln(1/\gamma)}[\ln(Db/\delta)]_+-2n
\right]_+.
```

The stronger finite theorem keeps the exact product of conditional survival probabilities, adds the existing all-cut lower bound on the coherent block size, and retains the exiting-rank floor in the gate ledger. These improvements change finite necessary counts without changing the fixed-noise asymptotic order. The complete-erasure/Jensen reliability mechanism has an established cumulative-spacetime precedent; its adaptation here uses backward Pauli rank and the T-growth ledger. For delayed-basis guessing, the reliability part applies with $`\delta=(1-p)/(1-2^{-k})`$ for every $`2^{-k}\lt p\lt1`$. Its coherent-block contribution must use the separate BB84 dimension bound, not the Bell-fidelity floor. The resulting bound can witness work retaining classical as well as quantum input dependence, without verifying the device assumptions.

## 3. Quantum storage at zero internal T count

The [Clifford-only construction](clifford_lifetime.md) extracts a stabilizer-syndrome history without performing corrections during storage. Each stored syndrome bit has several replicas, each with a growing ternary tree of noisy descendants. At the endpoint, majority reconstruction recovers the history, and successive syndrome differences identify the correctable error increments.

This construction preserves an unknown quantum state and its reference correlations. Only the syndrome record is copied. In each Pauli-noise branch, syndrome extraction acts on a definite syndrome sector and reveals no logical amplitudes.

All records remain physical noisy qubits, including old parents and their fresh children while they coexist at a branching boundary. Retired slots also remain noisy until reset for reuse. This counted width permits decoding to be deferred until the endpoint. The resulting lifetime order matches the converse:

```math
D_{\max}(n)=\Theta(\gamma^{-1}\ln n)
```

at fixed $`k`$ and target $`4^{-k}\lt F_0\lt1`$, for $`0\lt\gamma\le1/16`$ and sufficiently large width. The constants depend on the fixed target parameters; the exact width exponent is not optimized.

The endpoint recovery is unrestricted and may be computationally expensive. This is part of the operational statement: the memory must preserve recoverability during storage, while the final computation needed to recover the state is outside the interval.

## 4. A circuit that performs the correction during storage

The [coherent recovery construction](coherent_recovery.md) uses a different allocation of resources. A five-qubit Clifford decoder exposes four syndrome bits and one logical wire. Two quadratic Boolean polynomials specify the logical Pauli correction. Each polynomial has three quadratic terms, so six Toffoli-equivalent gates suffice. Exact seven-T decompositions give the stated upper bound of 42 T/T-dagger gates per decoding node.

Concatenating the code gives the conservative recursion $`p_{\ell+1}\le10p_\ell^2`$ for residual Pauli-label error probabilities. Batching order $`1/\gamma`$ weak-noise rounds between recoveries prevents the operation count from paying a constant correction cost in every arbitrarily short noise step. For integers $`k,D\ge1`$, $`0\lt\gamma\le1/20`$, and $`0\lt\epsilon\lt1/2`$, set $`\alpha=\log_2 5`$ and

```math
\ell=1+\ln\frac{k(\gamma D+1)}{\epsilon}.
```

The explicit construction has

```math
n=O(k\ell^\alpha),
\qquad
m_T=O\bigl(k(\gamma D+1)\ell^\alpha\bigr),
\qquad
\frac12\|\Lambda-\mathrm{id}_{2^k}\|_\diamond\le\epsilon.
```

All correction controls are counted quantum registers. The decoded syndrome wires are reset and reused for reencoding, so the construction needs no work qubits beyond the code block. It provides a physical circuit witness inside the ideal processing model, rather than an assumption that classical decoder running time equals coherent gate complexity.

This five-qubit witness retains logarithmic factors and conservative constants over its stated noise range.

The [Tanner construction](tanner_memory.md) uses a constant-rate, linear-distance code and a fixed number of iterations of the inherited parallel decoder. The iterations need only preserve a correctable residual region; exact return to the code space is deferred to the free endpoint. Every local minimization, maximum-weight choice, decomposition and tie is a fixed finite Boolean table. Reversible compilation and hardwired incidence give linear T count and linear workspace per correction batch. Syndrome and scratch are reset on a fixed schedule and remain part of the noisy width.

For a sufficiently small fixed code-dependent noise range, this gives

```math
n=O\!\left(k+\ln\frac{J}{\epsilon}\right),
\qquad m_T=O(nJ),
\qquad J=O(\gamma D+1),
\qquad\frac12\|\Lambda-\mathrm{id}_{2^k}\|_\diamond\le\epsilon.
```

Independent depolarizing noise makes each error component spatially independent within a batch. The X and Z components on a qubit remain correlated; the proof uses a union bound. A bounded-ratio family of code lengths permits arbitrary payloads by fixing unused logical qubits, without puncturing the code.

The construction is available throughout $`0\lt\gamma\le1/20`$ after a fixed number of inner five-qubit levels. Combining it with the reliability and payload converses now gives matching processing order for every payload at fixed positive noise. More precisely, fix $`\gamma\in(0,1/20]`$, $`C\gt0`$ and $`0\lt\epsilon_0\lt1/2`$. Along $`D\to\infty`$, with $`0\lt\epsilon\le\epsilon_0`$ and width $`n\le C[k+\ln(D/\epsilon)]`$, every permitted memory requires

```math
m_T=\Omega_{\gamma,C,\epsilon_0}\!\left(D[k+\ln(D/\epsilon)]\right).
```

The Tanner construction attains the same order at width $`O(k+\ln(D/\epsilon))`$. A fixed one-qubit payload therefore still requires order $`D\ln(D/\epsilon)`$ processing at such width. For payloads dominating the logarithm this reduces to the earlier order $`k\gamma D`$, with positive $`\gamma`$ fixed. The constants depend on the chosen noise strength, width cap and error ceiling as well as the fixed code family. They can be extremely large. The extension is proved by a coherent decode/correct, outer-recovery and reencode schedule, without an isotropic effective-noise assumption; its inner level is independent of payload, duration and target error. This does not settle the uniform weak-noise dependence or numerical optimum.

## 5. More gates do not eliminate the need for width

There is a simple obstruction even with unrestricted processing. In the stochastic decomposition of one noise layer, every one of the $`n`$ memory qubits is completely depolarized with probability $`\gamma^n`$. If this occurs at any round, the entire input dependence is lost. Subsequent operations and independent replacement states cannot recreate it.

The probability that no complete-register erasure occurs is $`(1-\gamma^n)^D`$. Every branch containing such an erasure is a replacement channel from the logical input to the final output, so its Bell fidelity is exactly $`K^{-2}=4^{-k}`$. Other branches have fidelity at most one. Therefore every protocol in this model satisfies

```math
F_{\mathrm e}
\le
4^{-k}+(1-4^{-k})(1-\gamma^n)^D.
\qquad\text{(I1)}
```

This argument permits arbitrary internal gates; it does not depend on the Clifford restriction or the T budget.

For $`D\ge1`$ and $`0\lt\epsilon\lt1-4^{-k}`$, requiring $`F_{\mathrm e}\ge1-\epsilon`$ in (I1) gives

```math
(1-\gamma^n)^D
\ge
1-\frac{\epsilon}{1-4^{-k}}.
```

Set $`\eta=\epsilon/(1-4^{-k})`$. Using $`\ln(1-x)\le-x`$ gives the necessary consequence

```math
n\ge\frac{\ln D-\ln[-\ln(1-\eta)]}{|\ln\gamma|}.
\qquad\text{(I2)}
```

When its right side is useful, this makes the width requirement explicit. At fixed $`\gamma`$ and small error, it grows logarithmically in $`D/\epsilon`$. The continuing-processing slope concerns feasible families with adequate width; fixed finite width still has the lifetime ceiling (I1).

The five-qubit construction uses a power of a logarithm in the storage horizon. With fixed-level inner protection throughout $`\gamma\le1/20`$, the Tanner construction's width instead has the additive form $`O(k+\ln(J_+/\epsilon))`$, with $`J_+=\max\{1,J\}`$. At fixed noise, the channel-dimension bound also forces $`n\ge k`$ when $`\epsilon\lt1/2`$. Together with (I2) and the Tanner upper bound, this proves $`n_{\min}=\Theta_\gamma(k+\ln(D/\epsilon))`$, uniformly over the stated payload and error range as $`D\to\infty`$. Thus the matching processing order holds within a sufficiently large fixed factor of the minimum possible width. The [full corollary](reliability_processing.md) makes the quantifiers explicit; the numerical and unrestricted width/work frontier remains open.

## 6. Endpoints, controllers, and magic states

An arbitrary final decoder may pull a single output Pauli back to many independent input Pauli labels. Its binary support rank can greatly exceed twice the number of logical output qubits. The general theorem therefore starts from $`2n`$, rather than silently using $`2k`$.

When the final decoder is Clifford followed by discarding unused wires, the terminal rank is $`2k`$. Any preceding nonlinear recovery then remains an internal charged operation. If that nonlinear recovery is moved outside the storage interval, its effect must instead be included in the decoder-support rank. This distinction explains why the constructions can use different endpoint strategies without contradicting the same converse.

An external noiseless classical syndrome decoder changes the model. It can store measurement results and evaluate nonlinear functions between noisy rounds without the counted device implementing those functions coherently. Ordinary measurement-based stabilizer correction is therefore not covered by a claim that it must consume T gates during memory.

Likewise, the theorem is not a magic-state consumption bound. Fresh nonstabilizer states are explicitly allowed as replacement states. Their scalar effect on the backward Pauli coefficients does not expand the binary label space. The charged resource is the internal operation that changes the allowed label dynamics.

## 7. Scientific reach and remaining questions

The [closest predecessor](../references/RELATED_WORK.md) established the qualitative Clifford-memory obstruction and a logarithmic-squared forgetting bound. The rank contraction supplies a single-logarithm bound stable under an arbitrary reference; the constructive archive establishes its quantum lifetime order. The channel-rank floor then turns a qualitative need for non-Clifford gates into a continuing, extensive requirement on their number and placement.

The five-qubit construction remains a simple finite-circuit witness over its explicit noise range. The Tanner compilation and rank-extinction converse give matching processing order for all payloads at fixed positive noise and width proportional to payload plus reliability logarithm. Gu et al. supply the good code and residual-error mechanism, and prior passive LDPC work supplies measurement-free correction with resettable ancillas. Bharti, Haug and Tanggara supply the complete-erasure/Jensen reliability mechanism for cumulative spacetime. The additions here are fixed-circuit accounting and a rank-based adaptation that charges continuing T/T-dagger work even though Clifford operations and resets are free. Optimal numerical coefficients, uniform weak-noise sharpness, useful finite constants and the full width/work frontier remain open.

All conclusions retain independent depolarizing noise, fixed schedules, counted persistent workspace, and ideal processing between noise layers. Correlated noise, gate faults inside the processing intervals, geometric locality, and an external controller require their own model and analysis.


