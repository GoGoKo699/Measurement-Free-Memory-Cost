# Finite width and processing comparison

The generator compares proved necessary bounds with the finite archive and five-qubit
storage constructions. It evaluates `k=1`, `gamma=epsilon=0.01`, and every integer
duration `D=1..1000`. The width panel shows the zero-T necessary width, the
noisy archive's sufficient width, and coherent recovery's width. The processing
panel evaluates the necessary T budget at that same coherent width and compares
it with the construction's upper count.

Run from the repository root:

```bash
python3 -m pip install -r verification/Graph/requirements.txt
python3 verification/Graph/generate_plots.py --check
```

The check regenerates into temporary directories, verifies all 1,000 finite
error certificates, and compares the vector PDF, CSV, and analytical certificate
fields with the committed [references](expected/). Pinned dependencies are in
[requirements.txt](requirements.txt). Add `--strict-provenance` only with the
original Python 3.12.14 dependency environment recorded in the resource
certificate to require identical runtime-version metadata too. The current
pins also reproduce the processing lower-bound PNG; the older resource
certificate records earlier Pillow and python-dateutil versions. Its numerical
fields and vector PDF are compared independently of that version metadata.

Run without `--check` to retain fresh output under `generated/`. The
[reproduction guide](../README.md#resource-comparison-plot) describes output
paths, overrides, reference protection, and the finite example.

The curves use the [processing law](../../science/control_cost.md),
[archive construction](../../science/clifford_lifetime.md), and
[coherent recovery](../../science/coherent_recovery.md). These plotted constructions retain their explicit finite constants. The [Tanner construction](../../science/tanner_memory.md), combined with the [reliability converse](../../science/reliability_processing.md), has matching processing order across payload sizes at fixed positive noise in its proved range and width proportional to payload plus reliability logarithm. Its potentially enormous family constants are not plotted as practical estimates.


The resource curve uses the mean-rank ledger stated in `science/control_cost.md`. The additional exponential-moment family in [moment processing](../../science/moment_processing.md) can strengthen a conditional lower bound; the rank-extinction reliability bound is another complementary converse. The plotted curve is not asserted to optimize these bounds.

## Finite processing lower bounds

[generate_figures.py](generate_figures.py) produces the separate
[processing lower-bound figure](expected/memory_bounds.pdf) at `k=1`, `n=200`,
`gamma=epsilon=0.01`, for durations from `D=10^4` to `10^6`. It compares the
refined finite reliability bound with the rank-only reliability and payload
bounds, using the same prescribed internal T/T-dagger count `m_T` throughout.
The fixed-noise Tanner order is an analytical comparison; its family constants
are not fitted or plotted here.

```bash
python3 -m pip install -r verification/Graph/requirements.txt
python3 verification/Graph/generate_figures.py --check
```

The command regenerates the PDF, PNG and [numerical rows](expected/memory_bounds.json)
in a temporary directory and compares them byte for byte with the committed
references. It independently evaluates the two labelled values with an
80-digit direct root, then checks them against the outward rational certificates
in `verification/expected/reliability_examples.json`. Those certificates are
recomputed by the main scientific suite. At `D=100000` and `D=1000000`, the
refined necessary counts are `3973` and `48238`; the corresponding rank-only
counts are `2989` and `38422`.

Run without `--check` to retain fresh files under
`verification/Graph/generated/processing-bounds/`, or select another destination
with `--output-dir PATH`. The generator protects the committed `expected/`
directory from output overrides. A check plus `--output-dir PATH` also retains
the checked fresh files for downstream use.

The source formulas are proved in [reliability processing](../../science/reliability_processing.md)
and [control cost](../../science/control_cost.md). These curves are necessary
counts, not achieved circuit costs or a physical feasibility certificate.
