"""Run the SHARP replay once and then its nine deliberately corrupted variants.

From this directory, with the Python package python-flint available:

    python -B verify_sharp.py

(for example through   uv run --no-project --with python-flint python -B verify_sharp.py).

The base run of replay.py must exit with status 0 and print ALL CHECKS PASSED.
Each corrupted run must exit with status 1 and print a line that begins with
REJECT:. A run that fails for another reason, for instance a missing python-flint,
does not count as a rejection. This script exits with status 0 only if all ten
runs behave as required.

The script uses the standard library only. It is not one of the pinned members
of SOURCE-PINS.json and it does not change replay.py or cert_vectors.json.
"""
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

MUTANTS = [
    "identity",
    "cent_claim_up",
    "caus_claim_up",
    "floor_up",
    "kappa_up",
    "harnack_violate",
    "vector_tamper",
    "frac_up",
    "two_row_up",
]


def run(extra):
    proc = subprocess.run(
        [sys.executable, "-B", "replay.py", *extra],
        cwd=HERE,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=600,
    )
    return proc.returncode, proc.stdout, proc.stderr


def main():
    failures = []

    code, out, err = run([])
    ok = code == 0 and "ALL CHECKS PASSED" in out
    print(f"base replay: exit {code}: {'as required' if ok else 'FAILED'}")
    if not ok:
        failures.append("base")
        sys.stdout.write(out[-2000:])
        sys.stdout.write(err[-2000:])

    for name in MUTANTS:
        code, out, err = run(["--mutant", name])
        rejected = code == 1 and any(line.startswith("REJECT:") for line in out.splitlines())
        reason = ""
        if rejected:
            reason = next(line for line in out.splitlines() if line.startswith("REJECT:"))
        print(f"mutant {name}: exit {code}: {reason if rejected else 'NOT REJECTED'}")
        if not rejected:
            failures.append(name)
            sys.stdout.write(out[-1000:])
            sys.stdout.write(err[-1000:])

    if failures:
        print("FAIL:", ", ".join(failures))
        return 1
    print("PASS: base replay passes; nine corrupted runs rejected")
    return 0


if __name__ == "__main__":
    sys.exit(main())
