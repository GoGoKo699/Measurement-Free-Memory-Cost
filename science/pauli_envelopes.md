# Binary Pauli envelopes and reference-stable forgetting

This chapter proves the envelope invariant and logarithmic forgetting theorem. The resource consequences for a memory that continues to transmit quantum information are proved in [control_cost.md](control_cost.md). Operational assumptions and the boundary of the charged storage interval are specified in [model.md](model.md).

## 1. Setting and conventions

The memory occupies at most $`n`$ qubits, including every persistent workspace or controller register. Between noise layers its fixed, input-independent circuit may apply Clifford unitaries and replacement resets. A replacement reset on a subsystem $`A`$ has the form


```math
\mathcal R_{\tau,A}(X)=\tau_A\otimes\mathrm{Tr}_A X,
\tag{E1}
```


with tensor factors returned to their original order. The prepared state $`\tau_A`$ is fixed and independent of the retained subsystem and the unknown logical input. It may be mixed, entangled within $`A`$, or nonstabilizer. The reset channel does not supply a hidden memory correlated with preparations at other times.

Each of $`D`$ noise layers applies, independently to every persistent qubit,


```math
\mathcal N_q(X)=qX+(1-q)\frac{I}{2}\mathrm{Tr}X,
\qquad 0<q<1.
\tag{E2}
```


The noise is independent across layers and locations. Write $`\gamma=1-q`$ and $`\kappa=-\ln q`$. Gates and resets between these layers are ideal operations in this model; gate depth is not silently converted into additional noise layers.

An arbitrary initial encoder and fixed final decoder may be placed outside the storage interval. There is no protected quantum bypass, free adaptive classical controller, or intermediate postselection. Any mutable classical information retained during storage must be represented in the counted noisy registers, with its processing implemented by the allowed circuit. Unrecorded Pauli measurements can be added without changing the proof below. General measurement and feedback operations require a separate resource accounting.

Let $`W_n=\mathbb F_2^{2n}`$, and label Hermitian tensor-product Paulis by $`v=(x,z)\in W_n`$, ignoring their overall phases. Their commutation relation is determined by


```math
[v,w]=x\cdot z'+z\cdot x'\pmod 2,
\qquad w=(x',z').
\tag{E3}
```


For an operator $`O`$, its Pauli support is the set of labels with nonzero coefficients in its Pauli expansion. For a binary subspace $`V\le W_n`$, write


```math
\mathcal A_V=\mathrm{span}_{\mathbb C}\{P_v:v\in V\}.
\tag{E4}
```


The representation of Clifford conjugation by a binary symplectic map is standard [Gottesman, 1998](https://arxiv.org/abs/quant-ph/9807006). The proofs below use that representation together with the linear label projection induced by a replacement reset.

## 2. The simultaneous backward envelope

Realize each factor of (E2) as an independent random choice between the identity channel, with probability $`q`$, and complete depolarization, with probability $`1-q`$. A complete realization $`b`$ specifies a CPTP circuit $`\Phi_b`$, and the original channel is the average $`\Phi=\mathbb E_b\Phi_b`$. The realization is only a mathematical decomposition; it is not a flag supplied to the memory or its decoder.

Propagate all terminal observables backward simultaneously. At a cut, an envelope is a binary subspace $`V`$ containing the Pauli support of every observable pulled back from the terminal output. The actual support of those observables need not itself be a binary subspace.

For the full physical output, initialize $`V=W_n`$. More generally, if a fixed decoder $`\mathcal D`$ maps the physical output to a system $`L`$, initialize


```math
V_{\mathrm{dec}}
=\mathrm{span}_{\mathbb F_2}
\bigcup_{P\text{ in an output Pauli basis}}
\mathrm{supp}_{\mathrm P}\!\left(\mathcal D^\dagger(P)\right),
\qquad M=\dim V_{\mathrm{dec}}\le 2n.
\tag{E5}
```


All statements below are valid with $`M=2n`$; the smaller value in (E5) is available when the decoder is specified. The Pauli-support rank of a general decoder is not bounded by twice the number of output qubits; see Section 7 of [control_cost.md](control_cost.md).

### Lemma E1 — Envelope invariant

The following backward updates preserve the envelope property.

1. **Clifford unitary.** For its binary symplectic action $`S`$, use $`V\mapsto SV`$. The dimension is unchanged.
2. **Replacement reset on $`A`$.** Let $`\pi_{A^c}`$ set the two coordinates of every qubit in $`A`$ to zero. Use $`V\mapsto\pi_{A^c}V`$. The dimension cannot increase.
3. **Complete depolarization on $`E`$.** Let $`\pi_E`$ retain the two coordinates of the qubits in $`E`$. Use $`V\mapsto V\cap\ker\pi_E`$. The dimension cannot increase.
4. **Unrecorded measurement of a Pauli $`P_w`$.** Use $`V\mapsto V\cap w^\perp`$, where $`w^\perp=\{v:[v,w]=0\}`$. The dimension cannot increase.

**Proof.** Clifford conjugation maps each Pauli to a phase times the Pauli with label $`Sv`$. For a reset, direct duality gives


```math
\mathcal R_{\tau,A}^\dagger(P_A\otimes P_B)
=\mathrm{Tr}(\tau_A P_A)I_A\otimes P_B.
\tag{E6}
```


The coefficient in (E6) can vanish, but every nonzero image has label $`\pi_{A^c}v`$. Ignoring zero coefficients only enlarges the envelope. This remains true for a joint reset to an arbitrary state on $`A`$. Complete depolarization kills precisely the Pauli terms with a nonidentity factor on an erased qubit. Finally, the unrecorded Pauli measurement channel has adjoint $`O\mapsto(O+P_wOP_w)/2`$, which preserves commuting Paulis and kills anticommuting ones. Each conclusion extends to linear combinations of Paulis. Induction over the backward circuit proves the invariant. $`\square`$

Replacement resets are essential to the scope of this argument. A general channel that maps individual Paulis to scalar multiples of other Paulis need not induce a linear map on their binary labels. Such a channel is not covered merely because its action is Pauli monomial.

### Lemma E2 — One noise layer contracts expected rank

Let $`V\le W_n`$ be any fixed subspace of dimension $`r`$. Independently include each physical qubit in $`E`$ with probability $`\gamma=1-q`$. Then


```math
\mathbb E_E\dim(V\cap\ker\pi_E)\le qr.
\tag{E7}
```


**Proof.** Choose a full-row-rank $`r\times2n`$ generator matrix $`G`$ for $`V`$. Choose $`r`$ linearly independent columns of $`G`$. Let $`G_E`$ be the restriction to the coordinate columns of the qubits in $`E`$, and let $`K_E`$ count the chosen independent columns whose qubits belong to $`E`$. Those columns remain independent inside $`G_E`$, so $`\mathrm{rank}G_E\ge K_E`$. Each chosen column is retained in $`G_E`$ with marginal probability $`\gamma`$. The two coordinate columns of a single qubit are erased together; independence of the chosen columns is not needed for the expectation. Thus


```math
\mathbb E\mathrm{rank}G_E
\ge\mathbb E K_E=\gamma r.
\tag{E8}
```


The kernel of the restricted coordinate map has dimension $`r-\mathrm{rank}G_E`$. Taking expectations proves (E7). $`\square`$

At a backward noise layer the current envelope is a function only of noise later in forward time. It is independent of the erasures sampled in the current layer. Consequently (E7) may be conditioned on the current envelope and iterated. Free operations from Lemma E1 can only decrease the rank between layers.

## 3. Exact replacement components

### Theorem E3 — Logarithmic forgetting, stable under an arbitrary reference

Let $`\Phi`$ be the storage channel, optionally composed with a fixed final decoder. Suppose the storage interval contains only the free primitives of Lemma E1 and $`D`$ noise layers. Initialize the envelope with dimension $`M`$, as in (E5), or take $`M=2n`$. For every fixed input density operator $`\omega_0`$, let $`\tau=\Phi(\omega_0)`$. Then


```math
\frac12\|\Phi-\mathcal R_\tau\|_\diamond
\le \min\{1,Mq^D\},
\qquad \mathcal R_\tau(X)=\mathrm{Tr}(X)\tau.
\tag{E9}
```


Equivalently, for every finite reference system $`R`$ and joint input state $`\omega_{MR}`$,


```math
\frac12\left\|
(\Phi\otimes\mathrm{id}_R)(\omega_{MR})
-\tau\otimes\omega_R
\right\|_1
\le \min\{1,Mq^D\}.
\tag{E10}
```


In particular, any two input states obey


```math
\|\Phi(\rho)-\Phi(\sigma)\|_1
\le 2\min\{1,Mq^D\}.
\tag{E11}
```


**Proof.** After propagating through the entire sampled storage circuit, write $`r_b=\dim V_b`$. Lemmas E1 and E2 imply


```math
\mathbb E_b r_b\le Mq^D,
\qquad
p:=\Pr(r_b>0)\le\min\{1,Mq^D\}.
\tag{E12}
```


The second inequality uses that a positive integer rank is at least one. If $`r_b=0`$, the adjoint image of every output observable is a scalar multiple of the identity. Therefore $`\Phi_b`$ is exactly a replacement channel: all input states have the same output, and trace preservation extends the formula to $`\Phi_b(X)=\mathrm{Tr}(X)\tau_b`$.

Separating the sampled circuits according to this event yields


```math
\Phi=(1-p)\mathcal R_{\tau_G}+p\Psi,
\tag{E13}
```


where $`\Psi`$ is a CPTP channel and the first term is the average of the replacement branches. The cases $`p=0`$ and $`p=1`$ are interpreted by dropping the absent term. Since $`\tau=\Phi(\omega_0)`$,


```math
\Phi-\mathcal R_\tau
=p\bigl(\Psi-\mathcal R_{\Psi(\omega_0)}\bigr).
\tag{E14}
```


The diamond distance between two channels is at most two, proving (E9) and hence (E10). Equation (E13) also gives $`\Phi(\rho)-\Phi(\sigma)=p[\Psi(\rho)-\Psi(\sigma)]`$, whose trace norm is at most $`2p`$. This proves (E11) without an extra factor of two. $`\square`$

For $`M>0`$, a sufficient depth for the right side of (E9) to be at most $`\delta\in(0,1)`$ is


```math
D\ge \frac{\ln(M/\delta)}{-\ln q}.
\tag{E15}
```


With $`M=2n`$, this is $`O(\gamma^{-1}\ln(n/\delta))`$ at weak noise. If $`M=0`$, the decoded channel is already a replacement channel. Arbitrary initial encoders and later decoders preserve the diamond estimate by composition. No factor depending on the dimension of a reference system is introduced.

## 4. Changing width and fresh systems

The fixed-width formulation covers a schedule with active width $`n_t\le n`$. Embed the schedule in $`n`$ slots. A newly admitted subsystem is prepared by a replacement reset of slots that were inactive; a discarded subsystem may be replaced by a fixed state and then ignored. Apply the noise to all $`n`$ slots. Noise on inactive slots has no operational effect because those slots are reset before admission and never used while inactive.

This padding argument requires every system that retains information across a noise layer to be among the noisy slots at that layer. A qubit cannot be kept in an exempt external store and later called fresh. Likewise, an external preparation source that holds a partner of an entangled state across several preparation times is a persistent system; it must be included if its correlations are available to the memory. Fixed joint preparations at one time are allowed by (E1).

The conservative choice $`n=\max_t n_t`$ yields every bound in this chapter. A more detailed width profile may strengthen an instance, but is not needed for the theorem.

## 5. Relation to the existing logarithmic-squared bound

Nelson, Rajakumar, and Gullans prove forgetting at depth $`O(\gamma^{-1}\log n\log(n/\delta))`$ for noisy Clifford circuits with resets. Their stochastic decomposition of depolarizing noise, backward action on Paulis, and replacement-branch reduction are the starting point for the present argument. Their published Theorem 1 and Observation 2 appear in [*Non-Clifford Gates Are Required for Long-Term Memory*, CCC 2026](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.CCC.2026.38), with the [full article](https://drops.dagstuhl.de/storage/00lipics/lipics-vol383-ccc2026/LIPIcs.CCC.2026.38/LIPIcs.CCC.2026.38.pdf) and [arXiv version](https://arxiv.org/abs/2510.08451) available online.

The additional ingredient here is the simultaneous binary linear envelope and its one-layer expected-rank contraction (E7). It replaces counting individual surviving Paulis. The replacement decomposition then yields the reference-stable statement (E9). The theorem permits arbitrary joint replacement states because their adjoint label action remains a coordinate projection.

The result concerns an explicitly bounded processing model. It does not assert that all implementations of quantum error correction require magic states or that a noiseless classical decoder is free during storage. The quantum-memory construction in [clifford_lifetime.md](clifford_lifetime.md) and the charged recovery construction in [coherent_recovery.md](coherent_recovery.md) address different ways to approach the resulting resource frontier.

## 6. Finite verification

The deterministic verification program in [verification/verify_pauli_algebra.py](../verification/verify_pauli_algebra.py) constructs dense two-qubit Pauli matrices independently of the circuit rank simulator. It checks every binary subspace of $`\mathbb F_2^4`$, the conditional expectations used in the companion chapter, the local rotation envelopes, and mixed backward circuits with nonstabilizer resets. Its numerical checks detect finite operator and sign errors; the dimension-independent assertions are established by the proofs above. The accompanying exact-rational verification checks the finite-depth scalar optimization used in [control_cost.md](control_cost.md).
