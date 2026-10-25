# Gottesman bridge to Measurement-Free Memory Cost

**Single educational anchor:** Daniel Gottesman, *An Introduction to Quantum Error Correction and Fault-Tolerant Quantum Computation*, [arXiv:0904.2557](https://arxiv.org/abs/0904.2557), preprint (2009), proceedings version (2010).

This source provides stabilizer, Clifford, syndrome, fault-tolerance, and concatenation foundations. The project's channel-rank, Pauli-envelope, processing-budget, and archive theorems are not attributed to this tutorial. Research attribution remains in [related work](../references/RELATED_WORK.md).

## Source-to-project route

| Topic in Gottesman | Relevant milestones | Project derivation |
|---|---|---|
| Quantum errors and correction criteria | M01–M04, M09 | [Model and reference-stable fidelity](../science/model.md) |
| Stabilizer codes, normalizer, syndrome extraction | M05–M08, M24–M26 | [Fully noisy syndrome archive](../science/clifford_lifetime.md) |
| Clifford group and Pauli transformations | M07, M13–M16 | [Binary envelope and forgetting](../science/pauli_envelopes.md) |
| Fault-tolerant gates and concatenation | M10–M11, M27–M28 | [Five-qubit coherent recovery](../science/coherent_recovery.md) |

## Translate the conventions

The stabilizer generator commutation signs form a syndrome; they do not reveal the unknown logical amplitudes. The project writes a Pauli label as `v=(x,z)` in a binary symplectic space. Clifford conjugation acts linearly on these labels. A *Pauli envelope* is a linear enclosure of all backward supports; the actual supports need not form a subspace.

The physical noise is `N_gamma=(1-gamma) id + gamma R_(I/2)`. The probability of a nonidentity Pauli error is `3 gamma/4`, while a computational-basis record flips with probability `gamma/2`. Every persistent data, syndrome, and workspace register is noisy. Ideal global Clifford processing and input-independent resets are allowed between noise layers; internal T and T-dagger gates are charged. The initial encoder and final decoder are separate free endpoints.

## The three missing bridges

**From Clifford dynamics to forgetting.** A noise layer contracts expected binary-envelope rank. Zero envelope means the sampled tail is a replacement channel, leading to reference-stable diamond-norm forgetting. Read [E1–E3](../science/pauli_envelopes.md).

**From algebra to required processing.** The binary-rank-`r` Pauli algebra has complex dimension `2^r`. Its rank bounds Bell fidelity; preserving logical entanglement therefore enforces a positive expected rank at every cut. Each T gate adds at most one generator, giving finite-time and interval processing lower bounds. Read [C1–C7](../science/control_cost.md).

**From syndrome correction to charged constructions.** One zero-T construction copies only syndrome records into a growing, noisy archive and defers nonlinear decoding to the endpoint. Another performs coherent five-qubit correction with at most 42 T/T-dagger gates per decoding node, followed by concatenation and batching. Read [Clifford lifetime](../science/clifford_lifetime.md) and [coherent recovery](../science/coherent_recovery.md).

These are exact statements only under the specified ideal processing model. No hardware gate-noise threshold, optimal T count, or matching constant-factor construction is claimed. Use the [30-milestone map](MILESTONES.md) to follow the complete story from prerequisites to resource frontier.
