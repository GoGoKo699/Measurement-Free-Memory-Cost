# Thirty milestones — Measurement-Free Memory Cost

The single educational anchor is **Watrous, The Theory of Quantum Information (2018)**, shared with Lattice Memory Certification. This reader-independent map allocates **12 background/preliminary milestones and 18 narrative/innovation milestones**. M01–M12 develop states, noise, records, controls and endpoints into the complete model. M13–M30 connect the envelope argument to the processing ledger, both constructions, the delayed-basis score corollary and the remaining resource gap. The phase label does not claim that inherited lemmas or corollaries are independent discoveries. Read alongside the [Watrous bridge](WATROUS_BRIDGE.md). Each milestone is an understanding target, not a personal test or progress record.

## M01 — Memory task and unknown inputs

**Phase:** Background / preliminaries. **Central question:** What must a memory preserve when its input is unknown and may be entangled? **Understanding target:** State the storage task as a recovered quantum channel and explain why preserving classical outcomes alone does not preserve an arbitrary input. **Dependencies and story role:** none; connect a local prerequisite to the complete memory task. **Conceptual connection:** An external reference tests whether the logical state and its quantum correlations survive together. **Source:** [model](../science/model.md), with foundational context in the [Watrous bridge](WATROUS_BRIDGE.md).

## M02 — Noisy physical width

**Phase:** Background / preliminaries. **Central question:** Which physical registers belong to the memory's noisy width? **Understanding target:** Count data, syndrome, work, controller and resource registers whenever they retain information between rounds, including slots later reset or reused. **Dependencies and story role:** M01; connect a local prerequisite to the complete memory task. **Conceptual connection:** Moving a record to an auxiliary register changes its location, not its exposure to noise. **Source:** [model](../science/model.md), with foundational context in the [Watrous bridge](WATROUS_BRIDGE.md).

## M03 — Depolarizing rounds

**Phase:** Background / preliminaries. **Central question:** How does one depolarizing round act on data and computational-basis records? **Understanding target:** Write the product noise channel, derive its Pauli probabilities and record-bit flip rate, and distinguish complete-depolarization strength from total nonidentity-error probability. **Dependencies and story role:** M02; connect a local prerequisite to the complete memory task. **Conceptual connection:** The same noise parameter enters quantum error accumulation and the classical archive's broadcast process. **Source:** [model](../science/model.md), with foundational context in the [Watrous bridge](WATROUS_BRIDGE.md).

## M04 — Bell fidelity and diamond error

**Phase:** Background / preliminaries. **Central question:** How do Bell fidelity and half-diamond error describe successful recovery? **Understanding target:** Define both criteria, derive the replacement-channel Bell overlap 4^(-k), and explain why a half-diamond guarantee implies the required Bell-fidelity bound. **Dependencies and story role:** M03; connect a local prerequisite to the complete memory task. **Conceptual connection:** The converse uses one entangled test input; the constructions protect all inputs and reference correlations. **Source:** [model](../science/model.md), with foundational context in the [Watrous bridge](WATROUS_BRIDGE.md).

## M05 — Stabilizer subspaces

**Phase:** Background / preliminaries. **Central question:** How does a commuting Pauli subspace define encoded quantum information? **Understanding target:** Relate independent stabilizer generators to the code space, identify logical Pauli representatives, and distinguish stabilizers from detectable errors. **Dependencies and story role:** M04; connect a local prerequisite to the complete memory task. **Conceptual connection:** Commuting syndrome observables can diagnose errors while leaving the logical amplitudes undetermined. **Source:** [coherent_recovery](../science/coherent_recovery.md), with foundational context in the [Watrous bridge](WATROUS_BRIDGE.md).

## M06 — Syndrome extraction

**Phase:** Background / preliminaries. **Central question:** Why can a syndrome be copied without copying the unknown logical state? **Understanding target:** Follow Clifford decoding, CNOT copying of syndrome coordinates and reencoding on a Pauli-error branch, identifying which registers contain logical information. **Dependencies and story role:** M05; connect a local prerequisite to the complete memory task. **Conceptual connection:** The archive stores error information extracted from the data, while the quantum payload remains in the code. **Source:** [clifford_lifetime](../science/clifford_lifetime.md), with foundational context in the [Watrous bridge](WATROUS_BRIDGE.md).

## M07 — Clifford symplectic action

**Phase:** Background / preliminaries. **Central question:** What structure does Clifford conjugation preserve on binary Pauli labels? **Understanding target:** Translate Pauli multiplication and commutation into binary addition and the symplectic form, then track these labels through a Clifford transformation. **Dependencies and story role:** M06; connect a local prerequisite to the complete memory task. **Conceptual connection:** A Pauli can spread across many qubits without creating additional independent binary directions. **Source:** [pauli_envelopes](../science/pauli_envelopes.md), with foundational context in the [Watrous bridge](WATROUS_BRIDGE.md).

## M08 — Replacement resets

**Phase:** Background / preliminaries. **Central question:** Why can replacement resets reduce but never increase backward-envelope rank? **Understanding target:** Compute the reset adjoint on a tensor-product Pauli and identify the coordinate projection, even when replacement coefficients vanish or the state is magic. **Dependencies and story role:** M07; connect a local prerequisite to the complete memory task. **Conceptual connection:** A magic replacement state can change Pauli coefficients while the surviving labels remain projections of existing directions. **Source:** [pauli_envelopes](../science/pauli_envelopes.md), with foundational context in the [Watrous bridge](WATROUS_BRIDGE.md).

## M09 — Distance and correctability

**Phase:** Background / preliminaries. **Central question:** How does code distance turn syndrome information into a correctable error estimate? **Understanding target:** Show why a low-weight error and a minimum-weight representative with the same syndrome differ by a stabilizer up to phase within the stated distance bound. **Dependencies and story role:** M08; connect a local prerequisite to the complete memory task. **Conceptual connection:** Correct recovery depends on the logical action of an error, rather than identifying its exact physical representative. **Source:** [clifford_lifetime](../science/clifford_lifetime.md), with foundational context in the [Watrous bridge](WATROUS_BRIDGE.md).

## M10 — Concatenation

**Phase:** Background / preliminaries. **Central question:** How does concatenating a small code suppress errors at successive levels? **Understanding target:** Draw the five-child decoding tree, explain why failure requires at least two erroneous children, and count physical leaves and internal recovery nodes. **Dependencies and story role:** M09; connect a local prerequisite to the complete memory task. **Conceptual connection:** Concatenation trades additional noisy width and repeated recovery work for a smaller residual logical-error probability. **Source:** [coherent_recovery](../science/coherent_recovery.md), with foundational context in the [Watrous bridge](WATROUS_BRIDGE.md).

## M11 — T count and magic supply

**Phase:** Background / preliminaries. **Central question:** What is charged when fresh magic states are free but internal processing is restricted? **Understanding target:** Separate the supply of replacement states from T/T-dagger gate count, and locate initial encoding, internal feedback and final decoding in the accounting. **Dependencies and story role:** M10; connect a local prerequisite to the complete memory task. **Conceptual connection:** A state resource and the coherent operation that processes syndrome information are different entries in this model. **Source:** [model](../science/model.md), with foundational context in the [Watrous bridge](WATROUS_BRIDGE.md).

## M12 — Complete operational model

**Phase:** Background / preliminaries. **Central question:** Which assumptions make the converse and both constructions comparable? **Understanding target:** Specify independent noise, ideal global processing, fixed schedules, counted persistent registers and unrestricted endpoints; distinguish noise rounds from elementary gate depth. **Dependencies and story role:** M11; connect a local prerequisite to the complete memory task. **Conceptual connection:** A resource comparison is meaningful only when every route uses the same noise and endpoint conventions. **Source:** [interpretation](../science/interpretation.md), with foundational context in the [Watrous bridge](WATROUS_BRIDGE.md).

## M13 — Backward Pauli supports

**Phase:** Narrative / innovation. **Central question:** How does an observable propagate backward through the noisy storage circuit? **Understanding target:** Track binary Pauli labels through Clifford conjugation, replacement resets and sampled complete-depolarization branches. Explain why the surviving supports need a linear enclosure. **Dependencies and story role:** M07 and M12; introduce the stochastic backward-propagation framework inherited from Nelson, Rajakumar and Gullans before the simultaneous rank estimate. **Conceptual connection:** Separate a sampled channel branch from a physical correction protocol. **Source:** [Pauli envelopes](../science/pauli_envelopes.md), [attribution](../references/RELATED_WORK.md) and the [Watrous bridge](WATROUS_BRIDGE.md).

## M14 — Expected-rank contraction

**Phase:** Narrative / innovation. **Central question:** Why does one noise layer contract the expected rank of a simultaneous Pauli envelope by at least a factor q? **Understanding target:** Select independent columns of a binary generator matrix and average over the independent qubit-survival variables. Identify where independence from the current noise layer enters. **Dependencies and story role:** M13; convert the backward-propagation framework into the rank balance used by the forgetting and processing bounds. **Conceptual connection:** Binary rank measures available operator directions rather than a stored classical record length. **Source:** [Expected-rank contraction](../science/pauli_envelopes.md) and [related work](../references/RELATED_WORK.md).

## M15 — Replacement branches

**Phase:** Narrative / innovation. **Central question:** Why does a zero backward envelope make a sampled channel a replacement channel? **Understanding target:** Use the adjoint action on all observables to show that the output is independent of the input, then separate zero-rank and surviving branches. **Dependencies and story role:** M14; connect the preceding proof ingredients toward the final resource comparison. **Conceptual connection:** The sampled branches are proof variables; their erasure history is unavailable to the memory or decoder. **Source:** [pauli_envelopes](../science/pauli_envelopes.md), with foundational context in the [Watrous bridge](WATROUS_BRIDGE.md).

## M16 — Reference-stable forgetting

**Phase:** Narrative / innovation. **Central question:** How does rank contraction imply forgetting even in the presence of a reference? **Understanding target:** Bound the probability of a nonzero envelope and convert the replacement mixture into a half-diamond estimate, using a valid terminal rank M. **Dependencies and story role:** M15; connect the preceding proof ingredients toward the final resource comparison. **Conceptual connection:** Forgetting the full input channel is stronger than comparing a chosen pair of unentangled input states. **Source:** [pauli_envelopes](../science/pauli_envelopes.md), with foundational context in the [Watrous bridge](WATROUS_BRIDGE.md).

## M17 — Pauli algebra blocks

**Phase:** Narrative / innovation. **Central question:** Which parts of a Pauli algebra can retain quantum rather than classical information? **Understanding target:** Separate anticommuting pairs from the commuting radical and identify the coherent matrix block, diagonal factor and total algebra dimension. **Dependencies and story role:** M16; connect the preceding proof ingredients toward the final resource comparison. **Conceptual connection:** Many commuting observable directions can retain classical labels while supporting little or no logical entanglement. **Source:** [control_cost](../science/control_cost.md), with foundational context in the [Watrous bridge](WATROUS_BRIDGE.md).

## M18 — Superoperator rank

**Phase:** Narrative / innovation. **Central question:** How does the tail channel's operator range constrain its recovered Bell fidelity? **Understanding target:** Bound superoperator rank by the Pauli algebra's dimension and derive the channel trace identity and rank inequality on the complete logical CPTP map. **Dependencies and story role:** M17; connect the preceding proof ingredients toward the final resource comparison. **Conceptual connection:** Superoperator rank concerns a linear map on operators, not the rank of an output state or Choi matrix. **Source:** [control_cost](../science/control_cost.md), with foundational context in the [Watrous bridge](WATROUS_BRIDGE.md).

## M19 — Fidelity-dependent rank floor

**Phase:** Narrative / innovation. **Central question:** Why must high recovered fidelity enforce an expected-rank floor at every cut? **Understanding target:** Apply the branchwise fidelity bound before averaging, use the secant inequality with the 4^(-k) baseline, and recover the positive floor a. **Dependencies and story role:** M18; connect the preceding proof ingredients toward the final resource comparison. **Conceptual connection:** A floor on average rank allows low-rank branches; it does not certify that a large envelope preserves a state. **Source:** [control_cost](../science/control_cost.md), with foundational context in the [Watrous bridge](WATROUS_BRIDGE.md).

## M20 — T-gate rank growth

**Phase:** Narrative / innovation. **Central question:** How many binary directions can one T or T-dagger gate add? **Understanding target:** Conjugate commuting and anticommuting Paulis through a Pauli-axis rotation and enclose the resulting supports by adding the rotation's label. **Dependencies and story role:** M19; connect the preceding proof ingredients toward the final resource comparison. **Conceptual connection:** Clifford gates and resets do not increase envelope rank; each counted T-type gate contributes at most one direction. **Source:** [control_cost](../science/control_cost.md), with foundational context in the [Watrous bridge](WATROUS_BRIDGE.md).

## M21 — Logarithmic budget

**Phase:** Narrative / innovation. **Central question:** How does a rank floor turn repeated contraction into a gate-count lower bound? **Understanding target:** Compare logarithmic rank loss at noise layers with bounded growth at counted gates, retaining the terminal-width allowance and the positive part. **Dependencies and story role:** M20; connect the preceding proof ingredients toward the final resource comparison. **Conceptual connection:** Physical width supplies a finite reserve; continuing work becomes necessary once the storage interval exhausts it. **Source:** [control_cost](../science/control_cost.md), with foundational context in the [Watrous bridge](WATROUS_BRIDGE.md).

## M22 — Finite-depth injection

**Phase:** Narrative / innovation. **Central question:** What finite-depth information is lost in the simpler logarithmic estimate? **Understanding target:** Solve the scalar injection relaxation by coasting on the initial reserve and then replenishing the floor; compare its certificate with integer-gate refinements. **Dependencies and story role:** M21; connect the preceding proof ingredients toward the final resource comparison. **Conceptual connection:** An optimal schedule for necessary scalar constraints is not automatically a realizable or optimal quantum circuit. **Source:** [control_cost](../science/control_cost.md), with foundational context in the [Watrous bridge](WATROUS_BRIDGE.md).

## M23 — Interval processing law

**Phase:** Narrative / innovation. **Central question:** Why must each interior storage interval satisfy its own processing constraint? **Understanding target:** Apply the rank balance to an interval with entering bound M_I, explain how later gates affect that envelope, and identify 2n as a safe bound. **Dependencies and story role:** M22; connect the preceding proof ingredients toward the final resource comparison. **Conceptual connection:** Work outside an interval cannot supply an unlimited rank reserve inside a long interval with no counted gates. **Source:** [control_cost](../science/control_cost.md), with foundational context in the [Watrous bridge](WATROUS_BRIDGE.md).

## M24 — Noisy syndrome archive

**Phase:** Narrative / innovation. **Central question:** How can a noisy branching archive retain syndrome history until final recovery? **Understanding target:** Track extraction, fanout, noise and slot reuse on the fixed schedule, including each root's first noisy block and every live or retired register. **Dependencies and story role:** M23; connect the preceding proof ingredients toward the final resource comparison. **Conceptual connection:** Archive width buys record redundancy; nonlinear majority decisions are deferred to the endpoint rather than performed during storage. **Source:** [clifford_lifetime](../science/clifford_lifetime.md), with foundational context in the [Watrous bridge](WATROUS_BRIDGE.md).

## M25 — Deferred quantum recovery

**Phase:** Narrative / innovation. **Central question:** How do recovered syndrome differences identify a correction for the final quantum state? **Understanding target:** Reconstruct consecutive error syndromes, choose representatives, and show that their cumulative correction restores the logical state on the input-independent good event. **Dependencies and story role:** M24; connect the preceding proof ingredients toward the final resource comparison. **Conceptual connection:** The final decoder can be nonlinear and costly while preserving the unknown input and its entanglement with a reference. **Source:** [clifford_lifetime](../science/clifford_lifetime.md), with foundational context in the [Watrous bridge](WATROUS_BRIDGE.md).

## M26 — Zero-T lifetime order

**Phase:** Narrative / innovation. **Central question:** When do the forgetting bound and the archive construction give matching lifetime order? **Understanding target:** Invert the archive's exponential width cost and compare it with reference-stable forgetting at fixed logical size, nontrivial target fidelity and the stated weak-noise range. **Dependencies and story role:** M25; connect the preceding proof ingredients toward the final resource comparison. **Conceptual connection:** Sharp logarithmic lifetime describes a width tradeoff; it does not provide polynomial-width storage for arbitrary duration. **Source:** [clifford_lifetime](../science/clifford_lifetime.md), with foundational context in the [Watrous bridge](WATROUS_BRIDGE.md).

## M27 — Five-qubit coherent correction

**Phase:** Narrative / innovation. **Central question:** How are the sixteen correctable syndromes turned into a charged recovery circuit? **Understanding target:** Derive the two Boolean correction polynomials, implement their six quadratic terms with exact seven-T decompositions, and account for the 42-T upper count and syndrome-wire reuse. **Dependencies and story role:** M08, M10 and M26; contrast online nonlinear correction with deferred archive decoding. **Conceptual connection:** The five-qubit code and seven-T synthesis are established ingredients; the explicit correction ledger supplies the construction used here. **Source:** [Coherent recovery](../science/coherent_recovery.md), [circuit checks](../verification/README.md) and [attribution](../references/RELATED_WORK.md).

## M28 — Batching and concatenation

**Phase:** Narrative / innovation. **Central question:** How do finite batching and concatenation parameters guarantee reliable storage? **Understanding target:** Derive the 10p_l^2 recursion from disjoint physical-error samples, then choose the integer batch size and concatenation level and compute the total width and T count. State the allowed k, D, gamma and epsilon domains. **Dependencies and story role:** M27; extend one exact recovery node to the complete storage interval. **Conceptual connection:** Independence concerns error labels, even when logical inputs are entangled. **Source:** [Concatenation and batching](../science/coherent_recovery.md) and the [Watrous bridge](WATROUS_BRIDGE.md).

## M29 — A delayed-basis guessing premise

**Phase:** Narrative / innovation. **Central question:** When can a whole-string guessing score replace the Bell-fidelity premise? **Understanding target:** Factor each full-output sampled tail through a CPTP quantum-plus-classical algebra, apply the finite guessing and dimension bounds before averaging, and transfer the resulting all-cut expected-rank floor to gate and contiguous-interval budgets with M=2n. **Dependencies and story role:** M12, M17, M19 and M23; connect the same resource ledger to a delayed measurement task. **Conceptual connection:** The device receives only the quantum challenge, X stays secret and bases arrive after storage. Noise flags are mathematical aids, and existing recovery bounds already yield the same large-block low-score consequence. **Source:** Watrous §§2.2–2.3 and 3.1.2; [guessing transfer and recovery comparison](../science/guessing_processing.md), [Watrous bridge](WATROUS_BRIDGE.md), and [Lattice application](https://github.com/GoGoKo699/Lattice-Memory-Certification/blob/main/science/noisy-processing.md).

## M30 — Proved resource bounds and remaining gap

**Phase:** Narrative / innovation. **Central question:** Which width and processing costs are necessary, which are achieved, and what still prevents indefinite storage? **Understanding target:** Compare the fidelity- and score-dependent processing laws, sharp fixed-payload zero-T lifetime, exponential-width archive and polylogarithmic coherent-recovery upper cost under the same noise and endpoint conventions; derive the catastrophic full-register replacement ceiling and explain why more gates cannot reverse complete loss of the input. **Dependencies and story role:** M12, M16, M21–M23, M26, M28 and M29; synthesize the complete scientific account. **Conceptual connection:** Necessary rank constraints and explicit constructions delimit a frontier without determining optimal physical T count. A score does not verify noise/controller premises or automatically import an LWE processing test. **Source:** Watrous channel foundations via the [bridge](WATROUS_BRIDGE.md); [interpretation](../science/interpretation.md), [finite resource comparison](../verification/Graph/README.md), [processing law](../science/control_cost.md), and [score scope](../science/guessing_processing.md).
