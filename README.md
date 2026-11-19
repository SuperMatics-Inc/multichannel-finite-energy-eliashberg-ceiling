# Moment-norm ceilings √(41/30)/(2π) and √(45/32)/(2π) for multichannel Eliashberg branches

Davis Rens  
SuperMatics Inc.  
ORCID: [0000-0001-5839-5197](https://orcid.org/0000-0001-5839-5197)

For reciprocal multichannel linearized Eliashberg models with bounded, energy-dependent channel densities of states, the paper proves

`k_B T_lin < C sqrt(rho)`.

Here `rho` is the norm of the glue moment matrix, `C = sqrt(41/30)/(2*pi)` on the centered normal class, and `C = sqrt(45/32)/(2*pi)` on the causal Harnack class specified in the manuscript. The constants are upper bounds, not claimed sharp. The result concerns a reduced model's normal-branch linear instability. It does not bound a general many-body transition temperature and makes no material, synthesis, or measurement claim.

Section 7 proves a limitation of the edge-wise reserve method on the paper's weight classes. Its best possible constant exceeds 0.185411 on the centered class and 0.188201 on the Harnack relaxation. These are method floors, not physical transition-temperature lower bounds. The second floor concerns the relaxation, not actual causal branches. Full multichannel self-consistency may evade this obstruction; the sharp scalar endpoint for unequal channel densities of states remains open.

## Paper and replay

- Paper: [10.5281/zenodo.24380829](https://doi.org/10.5281/zenodo.24380829).
- Replay supplement: [10.5281/zenodo.24380827](https://doi.org/10.5281/zenodo.24380827).
- Manuscript source: `paper/main.tex`, `paper/refs.bib`.
- Written analytic proof: `replay/OPTIMUM2/FINITE-ENERGY/RIGOR.md`.

Preserve the nested `replay/` layout. From `replay/OPTIMUM2/FINITE-ENERGY`, run:

```sh
python -I -B verify.py
```

These checkers use only the Python standard library. Success prints `PASS: two exact bases; twenty-one guarded mutants; finite Lean 4.35.0` and exits 0. The primary covering has 99 causal and 107 centered leaves; the independent covering has 245 and 287 leaves and 72 supplementary trial forms. All 21 corrupted variants must exit 1 at a named guard. Trial forms are samples, not universal certificates.

The byte-pinned driver calls `C:/Users/davis/.elan/bin/elan.exe run leanprover/lean4:v4.35.0 lean --threads=1 --memory=500 Finite.lean`. It has no skip or path-override flag. CI installs elan and the exact toolchain, then places a small shell adapter at that literal path relative to the replay working directory on Linux. The adapter raises only this command's memory cap from 500 to 2048 MB; the Linux toolchain exceeded the original 500 MB kernel allowance. All other arguments, the toolchain, and the pinned driver and proof bytes remain unchanged. No proof is skipped. On Windows that same path is absolute. On other machines, either provide the literal path as CI does or run the two Python checkers and Lean separately, explicitly recording any runtime difference from the admitted driver.

For the separate method-barrier certificate, install `python-flint==0.8.0`, then run from `replay/OPTIMUM2/SHARP`:

```sh
python -I -B verify_sharp.py
```

The base Arb interval replay at 192 bits must exit 0 and print `ALL CHECKS PASSED`. The wrapper also runs nine corrupted variants. Each must exit 1 with a named `REJECT:` guard; an import error does not count. The wrapper exits 0 only when all ten runs behave as required.

`replay/SOURCE-PINS.json` checks all 15 frozen members; six inherited proof pins are also checked internally. `ADMITTED.sha256` records the pre-publication manuscript and replay snapshot, whose manuscript predates the three release-link metadata additions and bibliography line-wrapping adjustment. `MANIFEST.sha256` hashes every release source file except itself. The frozen replay README retains a historical `../main.tex` path; the released manuscript is under `paper/`.

Build the paper with `tectonic --keep-logs main.tex` from `paper/`. CI checks a warning-free 15-page build and the ORCID, paper DOI, supplement DOI and repository links. The released PDF is deposited separately; the ZIP is the exact single-root source commit created with `git archive`.

## Review status and limitations

This preprint has not been peer reviewed. AI agents generated the proofs, certificates, code, and manuscript drafts and edits under D. Rens' direction and review. Separate AI agents performed the independent proof, certificate and manuscript reviews. The final manuscript delta review accepted the package on 19 November 2026. Independent AI review is not journal peer review.

The analytic theorem is proved in the manuscript under its explicit hypotheses. Exact rational and Arb interval replays certify finite inequalities and closure data. Lean checks finite polynomial identities and rational budgets; it does not formalize the infinite operator, Harnack or Poisson inequalities, the width cone, soft-limit locking, or the method-barrier infinite-tail argument. Historical proof-note labels remain unchanged for byte custody.

The paper and manuscript are licensed under CC BY 4.0. The replay code, Lean file, proof notes and repository infrastructure are licensed under MIT. Do not enable Zenodo's GitHub integration for this repository; the linked supplement already provides its DOI.
