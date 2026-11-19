# Replay bundle for the finite-energy moment-norm ceilings

This directory accompanies the paper in `../main.tex`. It reproduces, in exact rational arithmetic, the finite part of the proof of

`k_B T_lin < C sqrt(rho)`, with `C_cent = sqrt(41/30)/(2*pi)` on the centered class and `C_caus = sqrt(45/32)/(2*pi)` on the causal class.

The result concerns linearized Eliashberg theory for a reduced multichannel model. It is a theorem about that model, not a statement about any material, and no replay output is a measurement.

It also reproduces the finite certificate behind Section 7 of the paper, which shows how far an edge-reserve proof of this kind can go. That part is independent of the first and is described under `OPTIMUM2/SHARP` below.

## Layout

The scripts pin six inherited proof files by SHA-256 and look for them at fixed sibling paths two levels above the script directory (`closure.json`, key `source_sha256`). The layout below therefore mirrors the original tree and must not be flattened:

```
replay/
  SHARP-Q2/RIGOR.md
  CENTERED-GAP/RIGOR.md
  GENERAL/RIGOR.md
  OPTIMUM/RIGOR.md
  OPTIMUM2/RIGOR.md
  note-centered/PHYS-Q2/RIGOR.md
  OPTIMUM2/SOFT-LIMIT/RIGOR.md   the written proof of the soft-limit corollary of Section 6
  OPTIMUM2/FINITE-ENERGY/
    replay.py        primary binary covering and exterior estimates
    independent.py   independent quadtree covering and 6x6 trial forms
    verify.py        driver: both replays, 21 guarded mutants, Lean
    Finite.lean      finite Lean 4 lemmas
    closure.json     the six inherited pins
    RIGOR.md         the written proof the certificates complement
  OPTIMUM2/SHARP/
    replay.py          Arb certificates for Section 7 (needs python-flint)
    cert_vectors.json  the two trial vectors of 160 entries and the weight description
    verify_sharp.py    runner: the base replay and nine corrupted runs (not pinned)
  SOURCE-PINS.json   sha256 and size of every member against its source
```

Every member listed in `SOURCE-PINS.json` is a byte-exact copy of its source. `SOURCE-PINS.json` records, for each member, its own SHA-256 and size and the SHA-256 of the source file it was copied from, together with a match flag. The one file of the bundle that is not a member is `OPTIMUM2/SHARP/verify_sharp.py`, a short runner written for this bundle; it changes none of the pinned bytes.

## Run

From `OPTIMUM2/FINITE-ENERGY`:

```sh
python -I -B verify.py
```

The proof scripts of that directory use only the Python standard library. The Lean step of `verify.py` expects a local `elan` installation at the hard-coded path and runs `C:/Users/davis/.elan/bin/elan.exe run leanprover/lean4:v4.35.0 lean --threads=1 --memory=500 Finite.lean`. `verify.py` offers no option to override that path, so replaying on another machine with the same Lean toolchain means editing the path in `verify.py`, which changes its pinned bytes.

Expected exit codes:

| command | exit |
|---|---|
| `replay.py` | 0, prints 99 causal leaves and 107 centered leaves, both at maximum depth 17 |
| `independent.py` | 0, prints 245 causal leaves and 287 centered leaves, both at maximum depth 9, and 72 trial forms |
| `replay.py --mutant NAME`, 13 names | 1, with a `REJECT:` line |
| `independent.py --mutant NAME`, 8 names | 1, with a `REJECT:` line |
| `Finite.lean` | 0 |
| `verify.py` | 0, final line `PASS: two exact bases; twenty-one guarded mutants; finite Lean 4.35.0` |

From `OPTIMUM2/SHARP`, with the Python package `python-flint`, which provides the Arb interval arithmetic, for example through `uv`:

```sh
uv run --no-project --with python-flint python -B verify_sharp.py
```

The runner itself uses only the standard library, so any interpreter that can import `flint` will do. The replay was exercised with python-flint 0.9.0 on Python 3.13.9, and no other versions were tried. To run the steps one at a time, call `python -B replay.py` and `python -B replay.py --mutant NAME` with `NAME` one of

```
identity  cent_claim_up  caus_claim_up  floor_up  kappa_up
harnack_violate  vector_tamper  frac_up  two_row_up
```

| command | exit |
|---|---|
| `replay.py` in `OPTIMUM2/SHARP` | 0, ends with `ALL CHECKS PASSED` |
| `replay.py --mutant NAME`, the 9 names above | 1, with a `REJECT:` line |
| `verify_sharp.py` | 0, final line `PASS: base replay passes; nine corrupted runs rejected` |

`verify_sharp.py` is not called by `verify.py`, and neither script touches the pinned files of the other.

## What is and is not checked

The two coverings certify the rectangle inequalities on the compact domain `[1/3, 3] x [0, 1/2]` and the exterior inequalities outside it. The 72 six-by-six forms are a supplementary finite sample and not the certificate. Lean checks polynomial identities and rational budgets after clearing denominators. It does not formalize the infinite operator, the covering loop, the Harnack inequalities or the Poisson rearrangement; those steps are proved in `OPTIMUM2/FINITE-ENERGY/RIGOR.md` and in the paper.

The `OPTIMUM2/SHARP` replay evaluates the Rayleigh quotients of the two weight sequences of Section 7 in 192-bit Arb ball arithmetic, with the infinite normal-row sums bounded above by their telescoped tails. It checks the four lower bounds (1) and the decrease of both sequences for every row up to 4000, and it compares the resulting enclosures with the stated lower bounds, with the two-row quotients and with `(2 pi c)^2`, where `c = 0.1827262477723295176168720` is the upper endpoint of an enclosure of `C_inf` written into the script. That enclosure is the gamma = 2 row of the certificate record "Certified Eliashberg gamma-model linear-instability constants: replay certificate" (Zenodo, https://doi.org/10.5281/zenodo.24366341). Its enclosures are conditional on the analytic tail estimates in that record's written proof. This bundle does not recompute it, and the replay takes it as given. The reduction of Lemmas 7 and 8, the membership of the weights beyond row 4000 and the tail estimate itself are proved in Section 7 of the paper and are not replayed.

## Historical notes and review status

All fifteen pinned files above are retained byte-for-byte from the reviewed proof workspace so that the SHA-256 checks in `closure.json` and `SOURCE-PINS.json` remain intact. The bundled `RIGOR.md` files therefore retain internal working-note labels (`PROVED`, `CERTIFIED`, `NUMERICAL`, `OPEN`, `KILLED`, `[INFERENCE]`) and relative workspace paths from that tree, and the docstring of `OPTIMUM2/SHARP/replay.py` keeps its label `CERTIFIED` and its pointers to `RESULT.md`, `RIGOR.md` and `GAMMA/RESULT.md` of the proof workspace, which are not in this bundle; those tags are historical working-note annotations rather than claims of peer review, external material properties, or measurements. Where the pre-finalization tightest-leaf table in `OPTIMUM2/FINITE-ENERGY/RIGOR.md` differs from `../main.tex` Table 2, the values in `../main.tex` Table 2 match the bundled `replay.py` output. Separate AI agents reviewed the underlying notes, certificates, and manuscript; this preprint and supplement have not been peer reviewed.

`SOURCE-PINS.json` is a transport check. Nothing in this bundle uploads, sends or publishes anything.
