#!/usr/bin/env python3
"""Reproduce the finite processing lower-bound figure and its numerical checks.

This is a necessary-bound comparison, with no fitted construction constants.
Generation uses Matplotlib mathtext and requires no LaTeX installation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np

HERE = Path(__file__).resolve().parent
EXPECTED = HERE / "expected"
BLUE, ORANGE, GREEN, GRAY = "#0072B2", "#D55E00", "#009E73", "#666666"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.labelsize": 10, "axes.titlesize": 10,
    "xtick.labelsize": 9, "ytick.labelsize": 9,
    "legend.fontsize": 8.2, "axes.spines.top": False,
    "axes.spines.right": False, "axes.linewidth": 0.7,
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "savefig.facecolor": "white", "figure.facecolor": "white",
})


def refined_bound(D, k=1, n=200, gamma=0.01, epsilon=0.01):
    D = np.asarray(D, dtype=float)
    delta = epsilon / (1 - 4.0**(-k))
    b = 1 - delta
    a_r = 2*k*b
    a_s = k*max(0.0, (1-epsilon-2.0**(-k))/(1-2.0**(-k)))
    # Stable even for large D or small target error.
    L = b / np.log(1/gamma) * (-np.log(-np.expm1(np.log(b)/D)))
    raw = gamma*D*np.maximum(a_r, a_s+L)+a_r-2*n
    return np.ceil(np.maximum(0, raw))


def coarse_bound(D, k=1, n=200, gamma=0.01, epsilon=0.01):
    D = np.asarray(D, dtype=float)
    delta = epsilon / (1 - 4.0**(-k))
    b = 1 - delta
    rank = np.maximum(2*k, np.maximum(0, np.log(D*b/delta))/np.log(1/gamma))
    return np.ceil(np.maximum(0, gamma*D*b*rank-2*n))


def payload_bound(D, k=1, n=200, gamma=0.01, epsilon=0.01):
    D = np.asarray(D, dtype=float)
    a = 2*k*(1-epsilon-4.0**(-k))/(1-4.0**(-k))
    return np.ceil(a*np.maximum(0, D*(-np.log1p(-gamma))-np.log(2*n/a)))


def save(fig, target):
    target.parent.mkdir(parents=True, exist_ok=True)
    # Suppress time-bearing PDF metadata so regeneration is byte-stable.
    fig.savefig(target, bbox_inches="tight", pad_inches=0.035,
                metadata={"CreationDate": None, "ModDate": None,
                          "Creator": "generate_figures.py"})
    fig.savefig(target.with_suffix(".png"), dpi=220, bbox_inches="tight", pad_inches=0.035)
    plt.close(fig)


def memory_figure(output_dir):
    # A journal-sized quantitative panel. The matching-order summary belongs
    # in the caption and theorem, rather than a second, text-only panel.
    with plt.rc_context({
        "font.family": "STIXGeneral", "mathtext.fontset": "stix",
        "font.size": 8.5, "axes.labelsize": 9.3,
        "xtick.labelsize": 8.3, "ytick.labelsize": 8.3,
        "legend.fontsize": 7.9, "axes.linewidth": 0.65,
        "xtick.major.width": 0.65, "ytick.major.width": 0.65,
        "xtick.minor.width": 0.45, "ytick.minor.width": 0.45,
        "xtick.major.size": 3.0, "ytick.major.size": 3.0,
        "xtick.minor.size": 1.7, "ytick.minor.size": 1.7,
    }):
        fig, ax = plt.subplots(figsize=(3.35, 2.6), constrained_layout=True)
        D = np.unique(np.rint(np.geomspace(1e4, 1e6, 700))).astype(int)
        for values, color, ls, label in [
            (refined_bound(D), BLUE, "-", "Refined reliability"),
            (coarse_bound(D), ORANGE, "--", "Rank-only reliability"),
            (payload_bound(D), GRAY, ":", "Payload bound"),
        ]:
            ax.plot(D, np.where(values > 0, values, np.nan), color=color,
                    ls=ls, lw=1.4, label=label, solid_capstyle="round")
        ax.set(xscale="log", yscale="log", xlim=(1e4, 1.15e6),
               ylim=(10, 1e5), xlabel=r"Noisy rounds, $D$",
               ylabel=r"Lower bound on $m_T$")
        ax.grid(which="major", color="#E0E3E6", linewidth=0.45)
        ax.set_axisbelow(True)
        ax.legend(loc="upper left", frameon=False, borderaxespad=0.35,
                  handlelength=2.1, labelspacing=0.28, handletextpad=0.55)
        for x, y in [(100000, 3973), (1000000, 48238)]:
            ax.plot(x, y, marker="o", ms=3.5, mfc=BLUE,
                    mec="white", mew=0.45, zorder=5)
            offset = (-9, 9) if x < 1000000 else (-9, 1)
            ax.annotate(f"{y:,}", xy=(x, y), xytext=offset,
                        textcoords="offset points", ha="right", va="bottom",
                        fontsize=8.5, color=BLUE,
                        arrowprops={"arrowstyle": "-", "color": BLUE,
                                    "lw": 0.55, "shrinkA": 2, "shrinkB": 3})
        save(fig, Path(output_dir)/"memory_bounds.pdf")



def validate():
    """Compare stable plotting formulas with a direct 80-digit root evaluation."""
    mp.mp.dps = 80
    k, n = 1, 200
    gamma, eps = mp.mpf("0.01"), mp.mpf("0.01")
    delta = eps/(1-mp.mpf(4)**(-k)); b = 1-delta
    ar = 2*k*b
    ass = k*(1-eps-mp.mpf(2)**(-k))/(1-mp.mpf(2)**(-k))
    certified = json.loads((HERE.parent/"expected/reliability_examples.json").read_text())
    canonical = {int(x["parameters"]["D"]): x for x in certified["examples"]}
    examples = []
    for D in [100000, 1000000]:
        # Independent direct power instead of the plotting expm1 route.
        L = b/mp.log(1/gamma)*mp.log(1/(1-mp.root(b,D)))
        raw = gamma*D*max(ar,ass+L)+ar-2*n
        exact_ceiling = int(mp.ceil(max(mp.mpf(0),raw)))
        assert exact_ceiling == int(refined_bound(D))
        rank_only = int(coarse_bound(D))
        certificate = canonical[D]
        assert certificate["parameters"]["k"] == k
        assert certificate["parameters"]["n"] == n
        assert certificate["parameters"]["gamma"] == "1/100"
        assert certificate["parameters"]["epsilon"] == "1/100"
        assert exact_ceiling == certificate["sharpened_integer_gate_count"]
        assert rank_only == certificate["integer_gate_count"]
        enclosure = certificate["sharpened_joint_lower"]
        assert mp.mpf(enclosure["lower"]) <= raw <= mp.mpf(enclosure["upper"])
        examples.append({"D": D, "refined_raw": mp.nstr(raw,65),
                         "refined_integer": exact_ceiling,
                         "rank_only_integer": rank_only})
    D = np.unique(np.rint(np.geomspace(1e4, 1e6, 700))).astype(int)
    assert np.all(refined_bound(D) >= coarse_bound(D))
    rows = [[int(d), int(a), int(b), int(c)] for d, a, b, c in
            zip(D, refined_bound(D), coarse_bound(D), payload_bound(D))]
    return {"parameters": {"k": k, "n": n, "gamma": "0.01", "epsilon": "0.01"},
            "columns": ["D", "refined_reliability", "rank_only_reliability", "payload"],
            "curve_rows": rows, "labelled_examples": examples,
            "canonical_certificates": "verification/expected/reliability_examples.json",
            "scope": "Finite necessary bounds and displayed values. No physical feasibility, achievable cost or universal theorem is established by these checks."}


def generate(output_dir):
    result = validate()
    output_dir.mkdir(parents=True, exist_ok=True)
    memory_figure(output_dir)
    (output_dir/"memory_bounds.json").write_text(json.dumps(result, indent=2)+"\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true",
                        help="regenerate in a temporary directory and compare committed references")
    parser.add_argument("--output-dir", type=Path,
                        help="retain regenerated files outside the protected expected directory")
    args = parser.parse_args()
    if args.output_dir is not None:
        output = args.output_dir.resolve()
        if output == EXPECTED.resolve() or EXPECTED.resolve() in output.parents:
            parser.error("--output-dir cannot overwrite committed references")
    else:
        output = HERE/"generated/processing-bounds"
    if args.check:
        with tempfile.TemporaryDirectory(prefix="memory-bounds-") as tmp:
            fresh = Path(tmp)
            result = generate(fresh)
            checks = []
            for name in ["memory_bounds.pdf", "memory_bounds.png", "memory_bounds.json"]:
                actual = (fresh/name).read_bytes()
                expected = (EXPECTED/name).read_bytes()
                assert actual == expected, f"Generated {name} differs from the committed reference"
                checks.append({"path": str((EXPECTED/name).relative_to(HERE.parent.parent)),
                               "sha256": hashlib.sha256(actual).hexdigest()})
            if args.output_dir is not None:
                import shutil
                output.mkdir(parents=True, exist_ok=True)
                for item in fresh.iterdir():
                    shutil.copyfile(item, output/item.name)
        print(json.dumps({"status": "PASS", "checked_outputs": checks,
                          "labelled_examples": result["labelled_examples"],
                          "scope": result["scope"]}, indent=2))
    else:
        result = generate(output)
        print(json.dumps({"status": "PASS", "output_directory": str(output),
                          "labelled_examples": result["labelled_examples"],
                          "scope": result["scope"]}, indent=2))


if __name__ == "__main__":
    main()
