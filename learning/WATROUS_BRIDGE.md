# From Watrous to continuing memory processing

[Learning guide](README.md) · [30 milestones](MILESTONES.md) · [Memory model](../science/model.md)

## The physical question

A memory stores an unknown state while its data and syndrome records are repeatedly disturbed. Resetting qubits can remove entropy, but deciding how to correct a syndrome also requires processing. The question here is how much processing beyond Clifford operations is necessary when every persistent record belongs to the noisy device.

The model allows ideal global Clifford operations and replacement resets to arbitrary fixed states, including magic states. It charges every internal T or T-dagger gate. The initial encoder and final decoder may be arbitrary free CPTP channels. A noise round acts on every persistent qubit; elementary gate depth between rounds is not itself extra noisy time. The [exact model](../science/model.md) makes these choices part of the theorem.

## Source and edition

The single educational anchor is **John Watrous, *The Theory of Quantum Information*, Cambridge University Press (2018)**, DOI [10.1017/9781316848142](https://doi.org/10.1017/9781316848142). The author's [book page](https://cs.uwaterloo.ca/~watrous/TQI/) provides the January 2018 author copy and errata. Section numbers below follow that copy. It is for personal use and is not redistributed here.

The same anchor supports [Lattice Memory Certification](https://github.com/GoGoKo699/Lattice-Memory-Certification). Its route emphasizes discrimination and finite quantum storage; this route emphasizes channels, backward observables and a resource ledger. Neither route assigns the entire textbook.

## Anchor sections

| Watrous location | Foundation | Connection to this project |
|---|---|---|
| §§1.1.1–1.1.3 | Tensor products, operators, decompositions and norms | Pauli expansions, reference systems and binary algebra coordinates |
| §§1.2.1–1.2.2 | Convexity and probability | Independent noise choices, averaged rank and scalar bounds |
| §§2.1.1–2.1.3 | Classical and quantum registers, reductions and purification | Count data, syndrome and controller records without a protected bypass |
| §§2.2.1–2.2.3 | CPTP maps, Kraus representations and examples | Encoders, replacement resets, noisy storage and endpoint recovery |
| §§2.3.1–2.3.2 | Measurements | A final basis-dependent readout and its unconditional score |
| §3.1.2, including Theorem 3.10 | Ensemble discrimination and the square-root measurement | The inherited guessing theorem and recovery comparison |
| §§3.2–3.3 | Fidelity and channel distances | Bell-overlap fidelity, reference-stable error and forgetting |
| §§4.1.1 and 4.2.2 | Mixed-unitary channels, fixed points and spectra | Pauli twirls, conditional expectations and the channel-rank bottleneck |

Watrous supplies the finite-dimensional channel language. Binary symplectic Pauli algebra, stabilizer-code constructions, tensor-power hypercontractivity and the cited uncertainty/recovery inequalities are supplied in the scientific chapters and their primary references. Gottesman's work remains a scholarly source for stabilizer and Clifford methods, rather than a second educational anchor. The project does not attribute its envelope contraction or processing ledger to the textbook.

## Translate the conventions

Watrous's adjoint symbol $`*`$ is written $`\dagger`$ here. A channel maps states forward; its adjoint pulls observables backward. The equality $`\mathrm{Tr}[O\Phi(\rho)]=\mathrm{Tr}[\Phi^\dagger(O)\rho]`$ lets a final readout be studied at earlier cuts without changing any outcome probabilities.

The physical width is $`n`$, the quantum payload or BB84 challenge length is $`k`$, and the number of full noise rounds is $`D`$. The noise convention is

```math
\mathcal N_\gamma=(1-\gamma)\mathrm{id}+\gamma\mathcal R_{I/2},
\qquad q=1-\gamma,\qquad\kappa=-\ln q.
```

The total nonidentity Pauli-error probability is $`3\gamma/4`$; a computational-basis record flips with probability $`\gamma/2`$. These conventions cannot be exchanged silently. All persistent systems, including controller records and correlated resource registers, count toward $`n`$ and undergo the stated noise.

A Pauli label is $`v=(x,z)\in\mathbb F_2^{2n}`$. Multiplying Paulis adds their labels up to phase; Clifford conjugation acts by an invertible symplectic map. An envelope $`V`$ is a linear enclosure of all backward Pauli supports. The actual surviving support set need not be linear. Its binary dimension $`r=\dim V`$ differs from the complex dimension $`2^r`$ of its Pauli algebra.

Bell fidelity uses the squared-overlap convention $`F_{\rm e}=\langle\Omega|\rho|\Omega\rangle`$. Watrous's fidelity function is the unsquared root fidelity. The worst-case channel error is $`\tfrac12\|\Lambda-\mathrm{id}\|_\diamond`$. The constructions prove this stronger reference-stable criterion, whereas the main processing converse assumes only Bell-overlap fidelity.

## A two-qubit algebra example

Take a two-qubit cut with envelope $`V=\langle X_1,Z_1,Z_2\rangle`$. It has binary dimension three. The first two generators form an anticommuting pair; the last commutes with every generator. Thus

```math
\mathcal A_V\simeq M_2\otimes\mathrm{Diag}_2,
\qquad r=2s+c=3,
\qquad s=1,\quad c=1.
```

The trace-preserving conditional expectation onto this algebra dephases the second qubit in its computational basis. Its output can be represented by one quantum qubit and one classical bit. This is a genuine CPTP factorization, not a statement that three qubits of quantum information survive. The classical bit may retain input dependence, but it cannot carry a second coherent qubit.

Compare $`\langle Z_1,Z_2\rangle`$ and $`\langle X_1,Z_1\rangle`$. Both have binary rank two and complex algebra dimension four. The first is entirely classical; the second contains one coherent qubit. This explains why the science keeps the matrix-block size and total algebra dimension distinct before deriving a conservative rank-only bound.

## From backward observables to a processing ledger

Read [Pauli envelopes](../science/pauli_envelopes.md) for the complete proof. Clifford gates relabel an envelope without increasing its dimension; an input-independent replacement reset projects coordinates and also cannot increase dimension. A full independent depolarizing layer contracts its expected dimension by a factor at most $`q`$. The stochastic replacement choices used to prove this are channel decompositions, not measured flags delivered to a controller.

When the backward envelope becomes zero, the sampled tail loses every nonconstant input observable. It is therefore a replacement channel even for an input entangled with a reference. Averaging supplies the [diamond-norm forgetting law](../science/pauli_envelopes.md) and the logarithmic time allowance available without continuing charged processing.

The [control-cost chapter](../science/control_cost.md) turns preservation into a lower bound at every cut. The algebra of a rank-$`r`$ envelope has complex dimension $`2^r`$. A complete recovered logical CPTP channel through it has Bell fidelity at most $`\min\{1,2^{r-2k}\}`$. A convex secant followed by averaging yields a uniform expected-rank floor $`R_t\ge a`$. This floor constrains the entire circuit, including the cuts immediately before and after every primitive.

Each T or T-dagger gate adds at most one binary generator to a backward envelope. Summing logarithmic changes gives a continuing gate requirement after the finite width reserve. The exact discrete injection optimization gives a complementary bound, and the same ledger applies to every contiguous storage interval. Arbitrary endpoint decoding requires the conservative full-output reserve $`M=2n`$ unless a smaller decoder-specific envelope has actually been proved.

## From a delayed-basis score to the same ledger

The [guessing-to-processing account](../science/guessing_processing.md) supplies the full transfer and strongest inspected recovery comparison. The device receives only $`|X\rangle_\Theta`$, with $`X,\Theta\in\{0,1\}^k`$ uniform and independent. It receives $`\Theta`$ only after all charged storage ends. It never receives $`X`$. The score is the unconditional average probability of returning the entire string correctly.

At each cut, average earlier noise into a CPTP prefix. Fixing future replacement choices leaves a CPTP tail whose full-output adjoint range lies in the sampled Pauli algebra. That tail factors through $`s`$ quantum qubits and $`c`$ classical bits, with $`r=2s+c`$. The [Lattice guessing theorem](https://github.com/GoGoKo699/Lattice-Memory-Certification/blob/main/science/guessing-bound.md) gives a branch score at most $`\min\{1,2^{(s-k)/6}\}`$; ordinary uniform-input dimension counting also gives $`\min\{1,2^{s+c-k}\}`$.

Combining these bounds before averaging gives a rank floor even below one-half success. For $`k\ge8`$, a convenient form is

```math
R_t\ge a_k(p),
\qquad a_k(p)=2k\left[\frac{p-2^{-k}}{1-2^{-k}}\right]_+.
```

For example, $`k=12`$ and $`p=0.4`$ give $`a_k(p)\simeq9.5965`$. The smaller-$`k`$ bounds and the all-cut averaging proof are given in the authoritative chapter. A rank floor supplies a necessary gate budget only when the storage duration exceeds its width allowance.

Retaining spectral normalization in existing recovery inequalities reproduces the convenient mean-rank score floors when applied branchwise before averaging. The direct factorization explains the operational mechanism; these floors yield no additional gate or interval consequence over that recovery route. The test also does not establish that the assumed noise, counted controller or processing restrictions hold. Arbitrary uncharged feedback and postselection change the model, and the Lattice project's conditional LWE wrapper is not automatically inherited.

## Retaining the noise distribution

The [moment-processing proof](../science/moment_processing.md) uses the same counted device model. Choose independent coordinate pivots for a rank-$`r`$ Pauli envelope. At a noise layer, a surviving subspace can use only pivots on retained physical qubits. This gives, for $`1\lt z\le2`$,

```math
\mathbb E[z^{r'}\mid V]\le(1-q+qz^2)^{r/2}.
```

Jensen's inequality therefore makes $`X_z=\log_z\mathbb E[z^r]`$ contract by $`\alpha_z=\ln(1-q+qz^2)/(2\ln z)`$. A charged T-type gate still adds at most one. Applying the direct quantum-plus-classical branch score curve in the moment variable supplies one all-cut floor; fixing a recovery of the averaged channel supplies a complementary floor. The existing scalar optimization now applies to this potential. The source-to-project step is the pivot argument and its dynamical composition; it is not attributed to Watrous as a textbook theorem.

This strengthening keeps information lost by mean-rank concavification. It does not make the mathematical noise flags available to a controller. An interval entry bound must control $`X_z`$ itself; $`2n`$ is safe, while a bound only on $`\mathbb E[r]`$ need not be.

The [guessing-decay corollary](../science/guessing_decay.md) reverses the same scalar ledger: given width, duration and a prescribed gate budget, solve for the largest permitted score floor. Retaining the initial width gives a finite bound. Coasting for a vanishing fraction of a growing storage duration yields an explicit exponential rate for the score. A necessary-rate threshold of $`2(1-q)/q`$ results. The stronger full-profile route propagates $`\phi_j(\lambda)=\tfrac12\ln[1+q^j(e^{2\lambda}-1)]`$ and sums each prescribed slot's contributions across cuts before applying the score polygon. Its exponent reaches $`\ln2`$ at zero T rate and has the same positivity threshold; this converse does not construct a memory operating at that rate.

## Deferred and continuing recovery

The [Clifford-only construction](../science/clifford_lifetime.md) stores syndrome histories in growing trees of noisy records and defers nonlinear decoding to the unrestricted endpoint. It protects unknown quantum inputs, including their reference correlations, but pays exponential noisy width for logarithmic lifetime. Copying a syndrome record does not clone an unknown logical state.

The [coherent construction](../science/coherent_recovery.md) corrects internally. Each five-qubit decoding node uses at most 42 T/T-dagger gates without additional work qubits. Concatenation reduces logical error, and batching weak noise controls how often recovery is needed. The resulting width and work have a polylogarithmic reliability overhead relative to the extensive lower bound.

The [Tanner construction](../science/tanner_memory.md) replaces concatenation by an inherited good code and a constant-iteration parallel decoder. A finite local truth table must be synthesized as gates; a constant-time classical RAM statement is insufficient. Track all current-syndrome and temporary registers, compile nonlinear decisions into counted T/T-dagger slots, apply the controlled Pauli correction, and reset scratch on a fixed schedule. On each Pauli-noise branch the syndrome is independent of logical amplitudes, so discarding scratch preserves an entangled logical input. A deterministic residual invariant, Chernoff tails and a union bound over correlated X/Z components give half-diamond error at most epsilon. Bounded-ratio code lengths then yield width O(k+ln(J/epsilon)) and work O(nJ), with J=O(gamma D+1), for a sufficiently small direct noise range. A fixed number of inner five-qubit levels reduces the fresh error per outer position to that tolerance, extending the same asymptotic bounds to $`\gamma\le1/20`$. The inner level is independent of payload, duration and target error; its scratch and T slots are included. No isotropic effective-noise assumption is required.

Watrous supplies CPTP maps, partial traces and reference-stable channel comparison. Gu et al. supply the code and repeated-round decoder mechanism; passive LDPC and measurement-free correction are also prior ingredients. The fixed-circuit compilation and resource composition must be proved separately. Matching order when the payload dominates the reliability logarithm leaves numerical coefficients and hardware performance open.

The final [resource interpretation](../science/interpretation.md) also keeps the catastrophic whole-register replacement ceiling explicit. Gates cannot reconstruct a challenge after every persistent copy has been erased. The [milestone map](MILESTONES.md) ends by connecting that feasibility ceiling, the archive, both coherent constructions and the score-based corollary under one operational model.
