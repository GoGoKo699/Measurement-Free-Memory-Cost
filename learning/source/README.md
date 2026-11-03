# Build the explanatory report

The report is a single self-contained LaTeX source. Its TikZ schematic is inline;
no external figure or bibliography file is required.

From the repository root, using pdfLaTeX from TeX Live 2025:

```bash
mkdir -p learning/report
pdflatex -interaction=nonstopmode -halt-on-error -output-directory learning/report learning/source/Measurement_Free_Memory_Cost_Report.tex
pdflatex -interaction=nonstopmode -halt-on-error -output-directory learning/report learning/source/Measurement_Free_Memory_Cost_Report.tex
```

The output is `learning/report/Measurement_Free_Memory_Cost_Report.pdf`. The
matching log and auxiliary files are written beside it and are not tracked.
Compile twice to resolve the contents and equation references. Standard packages
include AMS math, Latin Modern, geometry, microtype, booktabs, tabularx, TikZ,
fancyhdr, and hyperref.

The [verification workflow](../../.github/workflows/verification.yml) selects a
historic TeX Live 2025 distribution explicitly. It removes stale output before
compilation, checks the compiler and fresh log banners, and retains the built PDF,
matching source, log, compiler version, and source/PDF hashes as run artifacts.
This distinguishes a fresh build from the PDF stored in the repository.

The report explains the model and proofs continuously. The topic pages under
[`science/`](../../science/model.md) retain the scholarly attribution, full
claims, and authoritative proof routes.
