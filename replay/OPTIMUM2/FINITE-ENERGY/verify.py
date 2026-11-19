"""Run the assigned rational replays, rejection guards, and finite Lean 4.35.0 verification.

This driver verifies both independent coverings, all twenty-one guarded mutants, and Finite.lean.
"""

from pathlib import Path
import subprocess
import sys

root = Path(__file__).resolve().parent


def run(args: list[str]) -> subprocess.CompletedProcess[str]:
    p = subprocess.run(args, cwd=root, capture_output=True, text=True)
    print(p.stdout, end="")
    print(p.stderr, end="", file=sys.stderr)
    print("EXIT", p.returncode, "COMMAND", repr(args), flush=True)
    return p


mutants = {
    "replay.py": (
        "head-determinant",
        "stage1-cascade",
        "hankel-tail",
        "centered-cone",
        "far-imbalance",
        "large-energy",
        "coverage-gap",
        "strict-budget",
        "closure-pin",
        "decimal-constant",
        "upper-harnack",
        "reflected-tail",
        "first-tail-threshold",
    ),
    "independent.py": (
        "quadtree-area",
        "completion-sign",
        "matrix-gain",
        "missing-normal-mass",
        "tangent-intercept",
        "flat-admissibility",
        "psd-schur",
        "two-sided-bathtub",
    ),
}

for filename, names in mutants.items():
    if run([sys.executable, "-I", "-B", filename]).returncode != 0:
        raise SystemExit("REJECT: base " + filename)
    for name in names:
        p = run([sys.executable, "-I", "-B", filename, "--mutant", name])
        if p.returncode != 1 or "REJECT:" not in p.stdout + p.stderr:
            raise SystemExit("REJECT: unexpected mutant result " + filename + " " + name)

lean = [r"C:/Users/davis/.elan/bin/elan.exe", "run", "leanprover/lean4:v4.35.0", "lean"]
if run(lean + ["--threads=1", "--memory=500", "Finite.lean"]).returncode != 0:
    raise SystemExit("REJECT: finite Lean")
print("PASS: two exact bases; twenty-one guarded mutants; finite Lean 4.35.0")
