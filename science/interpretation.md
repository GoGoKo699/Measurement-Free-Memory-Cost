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

For the general decoder-independent statement, the terminal envelope has rank at most $`2n`$. At fidelity $`F_{\mathrm e}\ge1-\epsilon`$, define

```math
a=2k\,\frac{1-\epsilon-4^{-k}}{1-4^{-k}}.
```

When $`a>0`$, the simple necessary budget is

```math
m_T\ge
a\left[D(-\ln q)-\ln(2n/a)\right]_+.
```

The first term grows with noisy storage duration. Physical width contributes the finite logarithmic allowance $`\ln(2n/a)`$. For fixed acceptable error and weak noise, long storage therefore requires a total T count of order $`k\gamma D`$. At finite encoding rate $`k/n=\Theta(1)`$, the allowance corresponds to only order $`1/\gamma`$ noisy rounds.

The exact finite-depth law in [Control cost](control_cost.md) improves the constants and gives a precise scalar optimization. Its optimum can first consume the initial rank reserve, then restore the minimum necessary rank before each subsequent noise layer. That is an optimization of necessary rank constraints. It does not construct a memory realizing those ranks.

The same argument applies to every contiguous interval. A protocol cannot satisfy an arbitrarily long interior interval solely by moving its charged operations before or after that interval.

## 3. Quantum storage at zero internal T count

The [Clifford-only construction](clifford_lifetime.md) extracts a stabilizer-syndrome history without performing corrections during storage. Each stored syndrome bit has several replicas, each with a growing ternary tree of noisy descendants. At the endpoint, majority reconstruction recovers the history, and successive syndrome differences identify the correctable error increments.

This construction preserves an unknown quantum state and its reference correlations. Only the syndrome record is copied. In each Pauli-noise branch, syndrome extraction acts on a definite syndrome sector and reveals no logical amplitudes.

All records remain physical noisy qubits. Their exponential growth is the resource that permits decoding to be deferred until the endpoint. The resulting lifetime order matches the converse:

```math
D_{\max}(n)=\Theta(\gamma^{-1}\ln n)
```

at fixed logical size and nontrivial fidelity target, within the stated weak-noise range. The constants and exact width exponent are not optimized.

The endpoint recovery is unrestricted and may be computationally expensive. This is part of the operational statement: the memory must preserve recoverability during storage, while the final computation needed to recover the state is outside the interval.

## 4. A circuit that performs the correction during storage

The [coherent recovery construction](coherent_recovery.md) uses a different allocation of resources. A five-qubit Clifford decoder exposes four syndrome bits and one logical wire. Two quadratic Boolean polynomials specify the logical Pauli correction. Each polynomial has three quadratic terms, so six Toffoli-equivalent gates suffice. Exact seven-T decompositions give the stated upper bound of 42 T/T-dagger gates per decoding node.

Concatenating the code suppresses the error quadratically per level. Batching order $`1/\gamma`$ weak-noise rounds between recoveries prevents the operation count from paying a constant correction cost in every arbitrarily short noise step. With $`\alpha=\log_2 5`$ and

```math
\ell=1+\ln\frac{k(\gamma D+1)}{\epsilon},
```

the explicit construction has

```math
n=O(k\ell^\alpha),
\qquad
m_T=O\bigl(k(\gamma D+1)\ell^\alpha\bigr).
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
\tag{I1}
```

This argument permits arbitrary internal gates; it does not depend on the Clifford restriction or the T budget.

For $`\epsilon<1-4^{-k}`$, requiring $`F_{\mathrm e}\ge1-\epsilon`$ in (I1) gives

```math
(1-\gamma^n)^D
\ge
1-\frac{\epsilon}{1-4^{-k}}.
```

At fixed $`\gamma`$ and small error, this requires at least logarithmic width in $`D/\epsilon`$. Thus the long-time slope of a gate-budget lower bound is a necessary resource statement for feasible families with adequate width. It does not establish indefinitely reliable storage at fixed finite width.

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
