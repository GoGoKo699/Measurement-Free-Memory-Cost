# Memory model and scientific results

## The question

A quantum memory must remove errors while retaining the unknown state and its entanglement with other systems. Supplying clean qubits removes entropy. Processing the syndrome determines which correction to apply. This project asks how much non-Clifford processing is necessary when that processing must occur inside the noisy device, without a protected external classical controller.

The [qualitative need for non-Clifford processing](../references/RELATED_WORK.md) in this setting is established prior work. The results here quantify its continuing cost, prove a sharp order for the lifetime achievable with Clifford processing alone, and give explicit constructions on both sides of the space-versus-processing comparison.

The model deliberately permits very powerful ideal operations between noise layers. This makes its converse apply independently of geometric locality, while placing the constructions in an ideal control model rather than a detailed hardware model.

## 1. Systems, noise, and allowed control

There are at most $`n`$ physical qubits at any instant during storage. This includes all data, syndrome records, temporary work qubits, controller registers, and any resource that retains information between noisy rounds. Variable active width is represented by a fixed pool of $`n`$ slots; a slot may be reset before reuse.

A storage circuit contains exactly $`D`$ noise layers. Each layer applies the product channel

```math
\mathcal N_\gamma^{\otimes n},
\qquad
\mathcal N_\gamma(X)
=(1-\gamma)X+\gamma\,\frac I2\mathrm{Tr}X.
\tag{M1}
```

The noise is independent across qubits and time. We write

```math
q=1-\gamma,\qquad
\kappa=-\ln q,
\qquad 0<\gamma<1.
```

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
\tag{M2}
```

with tensor factors restored to their original order. The replacement state is independent of the unknown input and of the retained register. It may be mixed, nonstabilizer, or entangled within the simultaneously prepared subsystem. An external source that retains correlations across different preparation times is a persistent memory and is not supplied for free.

All operations between successive noise layers are ideal in this model. Their elementary gate depth does not create additional noisy rounds automatically. An arbitrary global Clifford may therefore contain a large syndrome-extraction circuit. Likewise, all elementary non-Clifford gates inside an ideal recovery operation are counted even though the operation can be completed before the next noise layer.

There is no postselection and no data-dependent choice of the externally prescribed schedule. Coherent control using counted registers is allowed and must be implemented using the permitted gates. The proofs also allow unrecorded Pauli measurements as zero-charge channels, but no conclusion about arbitrary adaptive measurement-based controllers follows from that observation.

## 2. Transmission criterion and endpoints

Let $`K=2^k`$, with $`k\ge1`$. An arbitrary CPTP encoder $`\mathcal E`$ maps the logical input to the memory, $`\Phi_D`$ is its storage channel, and an arbitrary fixed CPTP decoder $`\mathcal D`$ returns $`k`$ logical qubits:

```math
\Lambda=\mathcal D\circ\Phi_D\circ\mathcal E.
\tag{M3}
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
\tag{M4}
```

This is the squared-overlap convention. The converse assumes only $`F_{\mathrm e}\ge1-\epsilon`$. The two constructions prove the stronger worst-case, reference-stable guarantee

```math
\delta_\diamond(\Lambda)
:=\frac12\|\Lambda-\mathrm{id}_K\|_\diamond
\le\epsilon,
\tag{M5}
```

which implies the fidelity condition.

Endpoint recovery can perform nonlinear syndrome decoding and use work registers after storage has finished. Such endpoint work supplies no noiseless information during the $`D`$ noisy rounds. If a nonlinear recovery is instead applied before storage ends, its processing belongs to the charged circuit. The [decoder-specific rank refinement](control_cost.md#7-the-final-decoder-and-the-meaning-of-the-terminal-rank) makes this boundary quantitative.

## 3. Principal results

### Reference-stable forgetting with Clifford operations and resets

For a storage circuit using only Clifford operations and replacement resets, and for any fixed physical input state $`\omega_0`$, set $`\tau=\Phi_D(\omega_0)`$. Then

```math
\frac12\|\Phi_D-\mathcal R_\tau\|_\diamond
\le\min\{1,2nq^D\}.
\tag{M6}
```

Thus all input dependence, including correlations with an arbitrary reference, is lost to accuracy $`\delta`$ after

```math
D\ge \frac{\ln(2n/\delta)}{-\ln(1-\gamma)}.
\tag{M7}
```

The proof in [Binary Pauli envelopes](pauli_envelopes.md) tracks the dimension of one binary linear envelope containing every pulled-back Pauli support. Resets project its coordinates; one noise layer contracts its expected rank by a factor $`q`$. A zero envelope identifies an exactly input-independent sampled channel.

### An extensive budget for non-Clifford processing

Suppose $`F_{\mathrm e}\ge1-\epsilon`$ and $`0\le\epsilon<1-4^{-k}`$. Define

```math
a=2k\,\frac{1-\epsilon-4^{-k}}{1-4^{-k}}>0.
\tag{M8}
```

Every Clifford+T memory obeys the simple bound

```math
m_T\ge
a\left[D\kappa-\ln\!\left(\frac{2n}{a}\right)\right]_+,
\qquad [x]_+=\max\{x,0\}.
\tag{M9}
```

The count includes both $`T`$ and $`T^\dagger`$. At fixed error bounded away from the trivial threshold, this is an extensive cost of order $`k\gamma D`$ once $`\gamma D`$ is sufficiently larger than $`\ln(n/k)`$. This is a necessary condition; it does not assert that spending those gates suffices.

[Control cost](control_cost.md) gives sharper integer-gate and finite-depth inequalities, an exact optimization of the corresponding scalar rank constraints, and a bound for every contiguous storage interval. It also replaces $`2n`$ by a smaller decoder-specific terminal rank when justified. A fixed Pauli-axis rotation has rank-growth charge at most one, while a general operation on $`s`$ qubits has a safe charge at most $`2s`$; those structural charges are not automatically physical T counts.

The improved fidelity floor (M8) follows from an elementary channel fact. A logical branch whose superoperator rank is at most $`2^r`$ has

```math
F_{\mathrm e}\le\min\{1,2^{r-2k}\}.
\tag{M10}
```

Indeed, the superoperator trace equals $`K^2 F_{\mathrm e}`$, all eigenvalues of a CPTP map have modulus at most one, and the number of its nonzero eigenvalues is at most its rank. Averaging the convex secant bound forces expected envelope rank at least $`a`$ at every cut. Every T gate can add at most one binary generator.

### Sharp order of the zero-T quantum lifetime

For fixed $`k`$ and fixed target fidelity # Memory model and scientific results

## The question

A quantum memory must remove errors while retaining the unknown state and its entanglement with other systems. Supplying clean qubits removes entropy. Processing the syndrome determines which correction to apply. This project asks how much non-Clifford processing is necessary when that processing must occur inside the noisy device, without a protected external classical controller.

The [qualitative need for non-Clifford processing](../references/RELATED_WORK.md) in this setting is established prior work. The results here quantify its continuing cost, prove a sharp order for the lifetime achievable with Clifford processing alone, and give explicit constructions on both sides of the space-versus-processing comparison.

The model deliberately permits very powerful ideal operations between noise layers. This makes its converse apply independently of geometric locality, while placing the constructions in an ideal control model rather than a detailed hardware model.

## 1. Systems, noise, and allowed control

There are at most $`n`$ physical qubits at any instant during storage. This includes all data, syndrome records, temporary work qubits, controller registers, and any resource that retains information between noisy rounds. Variable active width is represented by a fixed pool of $`n`$ slots; a slot may be reset before reuse.

A storage circuit contains exactly $`D`$ noise layers. Each layer applies the product channel

```math
\mathcal N_\gamma^{\otimes n},
\qquad
\mathcal N_\gamma(X)
=(1-\gamma)X+\gamma\,\frac I2\mathrm{Tr}X.
\tag{M1}
```

The noise is independent across qubits and time. We write

```math
q=1-\gamma,\qquad
\kappa=-\ln q,
\qquad 0<\gamma<1.
```

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
\tag{M2}
```

with tensor factors restored to their original order. The replacement state is independent of the unknown input and of the retained register. It may be mixed, nonstabilizer, or entangled within the simultaneously prepared subsystem. An external source that retains correlations across different preparation times is a persistent memory and is not supplied for free.

All operations between successive noise layers are ideal in this model. Their elementary gate depth does not create additional noisy rounds automatically. An arbitrary global Clifford may therefore contain a large syndrome-extraction circuit. Likewise, all elementary non-Clifford gates inside an ideal recovery operation are counted even though the operation can be completed before the next noise layer.

There is no postselection and no data-dependent choice of the externally prescribed schedule. Coherent control using counted registers is allowed and must be implemented using the permitted gates. The proofs also allow unrecorded Pauli measurements as zero-charge channels, but no conclusion about arbitrary adaptive measurement-based controllers follows from that observation.

## 2. Transmission criterion and endpoints

Let $`K=2^k`$, with $`k\ge1`$. An arbitrary CPTP encoder $`\mathcal E`$ maps the logical input to the memory, $`\Phi_D`$ is its storage channel, and an arbitrary fixed CPTP decoder $`\mathcal D`$ returns $`k`$ logical qubits:

```math
\Lambda=\mathcal D\circ\Phi_D\circ\mathcal E.
\tag{M3}
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
\tag{M4}
```

This is the squared-overlap convention. The converse assumes only $`F_{\mathrm e}\ge1-\epsilon`$. The two constructions prove the stronger worst-case, reference-stable guarantee

```math
\delta_\diamond(\Lambda)
:=\frac12\|\Lambda-\mathrm{id}_K\|_\diamond
\le\epsilon,
\tag{M5}
```

which implies the fidelity condition.

Endpoint recovery can perform nonlinear syndrome decoding and use work registers after storage has finished. Such endpoint work supplies no noiseless information during the $`D`$ noisy rounds. If a nonlinear recovery is instead applied before storage ends, its processing belongs to the charged circuit. The [decoder-specific rank refinement](control_cost.md#7-the-final-decoder-and-the-meaning-of-the-terminal-rank) makes this boundary quantitative.

## 3. Principal results

### Reference-stable forgetting with Clifford operations and resets

For a storage circuit using only Clifford operations and replacement resets, and for any fixed physical input state $`\omega_0`$, set $`\tau=\Phi_D(\omega_0)`$. Then

```math
\frac12\|\Phi_D-\mathcal R_\tau\|_\diamond
\le\min\{1,2nq^D\}.
\tag{M6}
```

Thus all input dependence, including correlations with an arbitrary reference, is lost to accuracy $`\delta`$ after

```math
D\ge \frac{\ln(2n/\delta)}{-\ln(1-\gamma)}.
\tag{M7}
```

The proof in [Binary Pauli envelopes](pauli_envelopes.md) tracks the dimension of one binary linear envelope containing every pulled-back Pauli support. Resets project its coordinates; one noise layer contracts its expected rank by a factor $`q`$. A zero envelope identifies an exactly input-independent sampled channel.

### An extensive budget for non-Clifford processing

Suppose $`F_{\mathrm e}\ge1-\epsilon`$ and $`0\le\epsilon<1-4^{-k}`$. Define

```math
a=2k\,\frac{1-\epsilon-4^{-k}}{1-4^{-k}}>0.
\tag{M8}
```

Every Clifford+T memory obeys the simple bound

```math
m_T\ge
a\left[D\kappa-\ln\!\left(\frac{2n}{a}\right)\right]_+,
\qquad [x]_+=\max\{x,0\}.
\tag{M9}
```

The count includes both $`T`$ and $`T^\dagger`$. At fixed error bounded away from the trivial threshold, this is an extensive cost of order $`k\gamma D`$ once $`\gamma D`$ is sufficiently larger than $`\ln(n/k)`$. This is a necessary condition; it does not assert that spending those gates suffices.

[Control cost](control_cost.md) gives sharper integer-gate and finite-depth inequalities, an exact optimization of the corresponding scalar rank constraints, and a bound for every contiguous storage interval. It also replaces $`2n`$ by a smaller decoder-specific terminal rank when justified. A fixed Pauli-axis rotation has rank-growth charge at most one, while a general operation on $`s`$ qubits has a safe charge at most $`2s`$; those structural charges are not automatically physical T counts.

The improved fidelity floor (M8) follows from an elementary channel fact. A logical branch whose superoperator rank is at most $`2^r`$ has

```math
F_{\mathrm e}\le\min\{1,2^{r-2k}\}.
\tag{M10}
```

Indeed, the superoperator trace equals $`K^2 F_{\mathrm e}`$, all eigenvalues of a CPTP map have modulus at most one, and the number of its nonzero eigenvalues is at most its rank. Averaging the convex secant bound forces expected envelope rank at least $`a`$ at every cut. Every T gate can add at most one binary generator.

### Sharp order of the zero-T quantum lifetime

For fixed $`k`$ and 4^{-k}<F_0<1`$, the largest lifetime achievable with Clifford operations and resets alone satisfies

```math
D_{\max}(n)=\Theta(\gamma^{-1}\ln n)
\tag{M11}
```

for sufficiently large width and $`0<\gamma\le1/16`$, with constants depending on the fixed target parameters.

The constructive direction in [Clifford quantum lifetime](clifford_lifetime.md) stores unknown quantum states, rather than only a classical bit. A stabilizer code produces a history of syndrome bits; many noisy copies of each bit are propagated through ternary trees. All majority voting and nonlinear decoding occur at the final endpoint. Every tree node remains a counted noisy qubit. Its exponential width and ideal global fanout are essential resources.

### Explicit coherent recovery within logarithmic factors

For $`0<\gamma\le1/20`$ and $`0<\epsilon<1/2`$, put

```math
\ell=1+\ln\frac{k(\gamma D+1)}{\epsilon},
\qquad \alpha=\log_2 5.
```

The [coherent recovery construction](coherent_recovery.md) achieves (M5) with

```math
n=O(k\ell^\alpha),
\qquad
m_T=O\bigl(k(\gamma D+1)\ell^\alpha\bigr).
\tag{M12}
```

It concatenates the five-qubit code, batches weak noise, and applies an explicit syndrome correction using at most 42 T/T-dagger gates per decoding node. No additional work qubits are needed beyond the code block. The constants and finite integer parameter choices are given in the construction.

This upper bound establishes a concrete comparison with (M9). The power-of-logarithm overhead has not been removed; neither optimal T count nor optimal physical width is asserted.

## 4. How the results fit together

The conclusions distinguish two ways to support storage:

| Choice of internal processing | Established resource behavior | Main physical cost |
|---|---|---|
| Clifford operations and resets only | Optimal lifetime order is logarithmic in width at fixed logical size and fidelity | An exponential number of noisy syndrome records |
| Charged coherent nonlinear recovery | Explicit width and gate count with a power-of-logarithm reliability overhead | Continuing T/T-dagger processing |
| Any high-fidelity protocol in the model | Must satisfy the gate-budget and interval lower bounds | Information cannot cross every temporal cut without sufficient internal processing |

Fresh magic states alone do not expand the backward binary envelope. Their availability therefore does not eliminate the processing lower bound. This is a statement about operations inside the specified memory, rather than a magic-state consumption theorem for all forms of error correction.

Both constructions permit ideal global processing. A local hardware design with faulty gates, measurement latency, or correlated noise must account for those additional features before importing an achievability conclusion.

## 5. Reading and verification routes

Read [Binary Pauli envelopes](pauli_envelopes.md) for the forgetting proof, then [Control cost](control_cost.md) for the fidelity bottleneck and gate budgets. The two constructive chapters, [Clifford quantum lifetime](clifford_lifetime.md) and [Coherent recovery](coherent_recovery.md), establish the different ways to pay for storage. [Interpretation](interpretation.md) compares their resources and explains the remaining optimality gap.

[Related work](../references/RELATED_WORK.md) identifies inherited ingredients, the closest published theorem, and the distinctions from measurement-based or externally decoded models. [Verification](../verification/README.md) describes the finite checks and their scope. Finite computations test the explicit algebra and circuits; universal claims rely on the proofs.
