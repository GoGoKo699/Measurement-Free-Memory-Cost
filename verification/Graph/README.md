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
[requirements.txt](requirements.txt). Add `--strict-provenance` on the original
Python 3.12.14 environment to require identical runtime-version metadata too.

Run without `--check` to retain fresh output under `generated/`. The
[reproduction guide](../README.md#resource-comparison-plot) describes output
paths, overrides, reference protection, and the finite example.

The curves use the [processing law](../../science/control_cost.md),
[archive construction](../../science/clifford_lifetime.md), and
[coherent recovery](../../science/coherent_recovery.md). These plotted constructions retain their explicit finite constants. The [Tanner construction](../../science/tanner_memory.md), combined with the [reliability converse](../../science/reliability_processing.md), has matching processing order across payload sizes at fixed positive noise in its proved range and width proportional to payload plus reliability logarithm. Its potentially enormous family constants are not plotted as practical estimates.


The resource curve uses the mean-rank ledger stated in `science/control_cost.md`. The additional exponential-moment family in [moment processing](../../science/moment_processing.md) can strengthen a conditional lower bound; the rank-extinction reliability bound is another complementary converse. The plotted curve is not asserted to optimize these bounds.
