# Interpretation: width, storage time, and internal processing

The [model](model.md) separates three resources: noisy physical width, fresh-state preparation, and non-Clifford operations performed while the unknown state is stored. The main conclusion is that an unlimited supply of replacement states does not remove the need for continuing non-Clifford processing at modest width.

The [forgetting theorem](pauli_envelopes.md), [processing budget](control_cost.md), and [two](clifford_lifetime.md) [constructions](coherent_recovery.md) quantify this statement in one operational model.

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

The [delayed-basis BB84 application](guessing_processing.md) replaces the fidelity premise by an unconditional whole-string guessing score. The input is only the quantum challenge, its secret label is never disclosed, and its basis is disclosed after storage. For $`k\ge8`$, score $`p\gt2^{-k}`$ forces $`a_k(p)=2k(p-2^{-k})/(1-2^{-k})`$ at every cut, so the same gate and interval inequalities apply. A positive score excess can reflect retained classical information and does not alone demonstrate quantum storage. Existing branchwise fidelity/recovery reductions, with spectral normalization retained, reproduce the convenient mean-rank score floors, including the small-block polygon. The direct floor explains the operational connection with its explicit controller and chronology assumptions. The [moment theorem](moment_processing.md) gives a further consequence in the same model: actual product-noise probabilities contract the logarithm of an exponential rank moment, while each charged gate has the same additive growth bound. Combining this dynamics with a fixed recovery of the averaged channel strengthens the gate and interval ledger. For fixed positive guessing score and fixed $`q`$, along growing challenges and storage durations with $`\ln(n/k)=o(D)`$, its necessary rate approaches $`2k(1-q)/q`$ per round. This improves the fixed-score coefficient of the convenient mean-rank law while preserving the order $`k\gamma D`$ at weak noise. The entering interval reserve must bound the logarithmic moment; a mean-rank reserve cannot automatically be reused.

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

The logarithmic factors and conservative constants remain. This witness does not claim a practical resource estimate, a constant-rate code, or an optimal decoder.

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

The explicit coherent construction uses a power of a logarithm in the storage horizon. Closing that exponent gap, and closing the corresponding gate-count gap, are distinct optimality questions.

## 6. Endpoints, controllers, and magic states

An arbitrary final decoder may pull a single output Pauli back to many independent input Pauli labels. Its binary support rank can greatly exceed twice the number of logical output qubits. The general theorem therefore starts from $`2n`$, rather than silently using $`2k`$.

When the final decoder is Clifford followed by discarding unused wires, the terminal rank is $`2k`$. Any preceding nonlinear recovery then remains an internal charged operation. If that nonlinear recovery is moved outside the storage interval, its effect must instead be included in the decoder-support rank. This distinction explains why the two constructions can use different endpoint strategies without contradicting the same converse.

An external noiseless classical syndrome decoder changes the model. It can store measurement results and evaluate nonlinear functions between noisy rounds without the counted device implementing those functions coherently. Ordinary measurement-based stabilizer correction is therefore not covered by a claim that it must consume T gates during memory.

Likewise, the theorem is not a magic-state consumption bound. Fresh nonstabilizer states are explicitly allowed as replacement states. Their scalar effect on the backward Pauli coefficients does not expand the binary label space. The charged resource is the internal operation that changes the allowed label dynamics.

## 7. Scientific reach and remaining gap

The [closest predecessor](../references/RELATED_WORK.md) established the qualitative Clifford-memory obstruction and a logarithmic-squared forgetting bound. The rank contraction supplies a single-logarithm bound stable under an arbitrary reference; the constructive archive establishes its quantum lifetime order. The channel-rank floor then turns a qualitative need for non-Clifford gates into a continuing, extensive requirement on their number and placement.

The constructive coherent recovery shows how to meet that requirement up to a power of a logarithm, with explicit finite circuits and a guarantee stronger than the converse's Bell-fidelity premise. The remaining gap is stated at the level of the full operational model: better coherent decoding could improve the achieved resource law. A restriction of this particular concatenated construction would not by itself prove that such an improvement is impossible.

All conclusions retain independent depolarizing noise, fixed schedules, counted persistent workspace, and ideal processing between noise layers. Correlated noise, gate faults inside the processing intervals, geometric locality, and an external controller require their own model and analysis.


