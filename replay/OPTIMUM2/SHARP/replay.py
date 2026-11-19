"""SHARP replay: re-verifies every CERTIFIED number in RESULT.md / RIGOR.md.

Run:
  UV_OFFLINE=1 UV_CACHE_DIR=D:/dev/.cache/uv uv run --no-project --with python-flint \
      python -B replay.py                 # must exit 0
  ... python -B replay.py --mutant NAME   # every mutant must exit 1

Certificates (exact rationals + 192-bit Arb balls, outward rounding):
  E1  exact bracket identities 2(1-19/60)=41/30, 2(1-19/64)=45/32.
  E2  centered instance: b_n = c_n = (2/pi) atan(W/w_n), W = 33/4, t = 1/10^6.
      Membership in the FINITE-ENERGY centered class (box DOS at widths
      sigma_m = w_m, i.e. the Proposition-1 lower envelope b_m = phi_{w_m}(b_0)),
      and the four Harnack bounds (checked for m < L, analytic beyond).
  E3  causal-relaxation instance: b_n = c_n = max(b_0/w_n, w_n b_0 - (w_n - 1)),
      b_0 = 7/8, t = 1/10^6; all four Harnack bounds checked exactly for m < L
      and by the closed form for m >= 3.
  E4  for each instance a rigorous lower bound R_lo on the Rayleigh quotient
      v^T A v / v^T v of the frozen-weight diagonal edge operator (X = Y), with
      the infinite normal-row sums bounded above by their exact telescoped tail.
      Then the FINITE-ENERGY edge reserve (eq. edge_reserve, r = 1) forces
      2(1 - d_*) >= R_lo, i.e. every edge-reserve constant is >= sqrt(R_lo)/(2 pi).
  E5  R_lo > kappa_hi = (2 pi C_inf_hi)^2 with the CERTIFIED upper endpoint
      C_inf_hi = 0.1827262477723295176168720 (GAMMA/RESULT.md, gamma = 2 row).
  E6  two-row witnesses: trial vectors supported on Matsubara rows n = 0, 1 only
      already give R_lo > kappa_hi for both instances.
  E7  consistency with FINITE-ENERGY's certified reserves, and the certified
      share of each open bracket [C_inf, C_class] lying below the edge-route floor.
"""
import json, os, sys
from fractions import Fraction as Fr

from flint import arb, ctx, fmpq

ctx.prec = 192
HERE = os.path.dirname(os.path.abspath(__file__))
MUT = sys.argv[sys.argv.index("--mutant") + 1] if "--mutant" in sys.argv else None
L = 4000  # explicit truncation of the normal-row sums; exact tail beyond

# claimed certified values (decimal rationals); RESULT.md quotes these
CLAIM = {"cent": Fr(1357186, 10**6), "caus": Fr(1398311, 10**6)}
FLOOR = {"cent": Fr(185411, 10**6), "caus": Fr(188201, 10**6)}  # lower bounds on sqrt(R)/(2 pi)
CINF_HI = "0.1827262477723295176168720"
if MUT == "cent_claim_up":
    CLAIM["cent"] += Fr(1, 10**4)
if MUT == "caus_claim_up":
    CLAIM["caus"] += Fr(1, 10**4)
if MUT == "floor_up":
    FLOOR["cent"] += Fr(1, 10**5)
if MUT == "kappa_up":
    CINF_HI = "0.1856"


def check(cond, msg):
    if not cond:
        print("REJECT:", msg)
        sys.exit(1)
    print("ok:", msg)


def q(x):
    return arb(fmpq(x.numerator, x.denominator))


# ---------------- E1: exact identities ----------------
r1 = Fr(41, 31) if MUT == "identity" else Fr(41, 30)
check(2 * (1 - Fr(19, 60)) == r1, "2(1-19/60) = 41/30")
check(2 * (1 - Fr(19, 64)) == Fr(45, 32), "2(1-19/64) = 45/32")

with open(os.path.join(HERE, "cert_vectors.json")) as f:
    DATA = json.load(f)
t = Fr(DATA["t"])
tA = q(t)
N0 = DATA["N0"]
J = L + N0 + 4
k = [None] + [1 / (tA + j * j) for j in range(1, J)]  # k_0 = 1/t cancels exactly (X = Y, b = c)
w = [2 * n + 1 for n in range(J)]
PI = arb.pi()


def weights(spec):
    """b_0..b_L as arb balls (b_L used as the monotone tail majorant)."""
    if spec["kind"] == "box":
        W = q(Fr(spec["W"]))
        return [2 / PI * (W / w[n]).atan() for n in range(L + 1)], None
    b0 = Fr(spec["b0"])
    if MUT == "harnack_violate":
        exact = [b0] + [max(b0 / w[n], w[n] * b0 - (w[n] - 1)) for n in range(1, L + 1)]
        exact[1] -= Fr(1, 100)
    else:
        exact = [max(b0 / w[n], w[n] * b0 - (w[n] - 1)) for n in range(L + 1)]
    return [q(x) for x in exact], exact


def membership(name, spec, b, exact):
    if spec["kind"] == "box":
        b0 = b[0]
        Wt = (PI * b0 / 2).tan()  # = W, so phi_{w_m}(b_0) = (2/pi) atan(W/w_m) = b_m
        check(b0 > 0 and b0 < 1, f"{name}: 0 < b_0 < 1")
        ok = True
        for m in range(1, L + 1):
            phi = 2 / PI * (Wt / w[m]).atan()
            ok &= (b[m] - phi).contains(0) or abs(b[m] - phi) < arb(10) ** -50
            wm = w[m]
            ok &= bool(b[m] >= b0 / wm) and bool(b[m] <= wm * b0)
            ok &= bool(b[m] >= wm * b0 - (wm - 1)) and bool(b[m] <= (b0 + wm - 1) / wm)
            ok &= bool(b[m] <= b[m - 1])
        check(ok, f"{name}: b_m = phi_(w_m)(b_0) (centered box at sigma_m = w_m), Harnack bounds, monotone, m <= {L}")
    else:
        b0 = exact[0]
        ok = True
        for m in range(1, L + 1):
            wm = w[m]
            bm = exact[m]
            ok &= max(b0 / wm, wm * b0 - (wm - 1), Fr(0)) <= bm <= min(wm * b0, (b0 + wm - 1) / wm, Fr(1))
            ok &= bm <= exact[m - 1]
        # m >= 3: b_m = b0/w_m (closed form); the four bounds reduce to 1/w <= w and b0 <= 1
        ok &= all(exact[m] == b0 / w[m] for m in range(3, L + 1))
        check(ok, f"{name}: four Harnack bounds and monotone decrease, exactly for m <= {L}; closed form b0/w_m for m >= 3")


def rayleigh_lower(b, v):
    n0 = len(v)
    va = [arb(x) for x in v]  # floats -> exact dyadic balls
    if MUT == "vector_tamper":
        va = [x * (1 + arb(i % 7) / 50) for i, x in enumerate(va)]
    sb = [b[n].sqrt() for n in range(n0)]
    sw = [arb(w[n]).sqrt() for n in range(n0)]
    u = [va[n] * sb[n] / sw[n] for n in range(n0)]
    quad = arb(0)
    for n in range(n0):
        row = arb(0)
        for m in range(n0):
            if m != n:
                row += (k[abs(n - m)] + k[n + m + 1]) * u[m]
        quad += u[n] * row
    for n in range(n0):
        D = arb(0)
        for m in range(L):
            if m != n:
                D += b[m] * (k[abs(n - m)] - k[n + m + 1])
        tail = sum((k[j] for j in range(L - n, L + n + 1)), arb(0))
        D += b[L] * tail  # b_m <= b_L for m >= L, p_nm >= 0, sum_{m>=L} p_nm = sum_{j=L-n}^{L+n} k_j
        diag = (2 * b[n] * k[2 * n + 1] - D) / w[n]
        quad += diag * va[n] * va[n]
    nv = sum((x * x for x in va), arb(0))
    return quad / nv


# E6 trial vectors on rows 0,1 (top eigenvectors of the leading 2x2 block, rounded; exact dyadics)
TWO_ROW = {"cent": [0.94548797, 0.32565701], "caus": [0.95361671, 0.30102355]}
CLAIM2 = {"cent": Fr(13329, 10**4), "caus": Fr(13870, 10**4)}
if MUT == "two_row_up":
    CLAIM2["caus"] = Fr(13871, 10**4)
kappa_hi = (2 * PI * arb(CINF_HI)) ** 2
print("kappa_hi =", kappa_hi.str(20))
for name in ("cent", "caus"):
    inst = DATA["instances"][name]
    b, exact = weights(inst["spec"])
    membership(name, inst["spec"], b, exact)
    R = rayleigh_lower(b, inst["v"])
    print(f"{name}: R ball =", R.str(25))
    check(bool(R > q(CLAIM[name])), f"{name}: R_lo > {CLAIM[name]} = {float(CLAIM[name])}")
    check(bool(q(CLAIM[name]) > kappa_hi), f"{name}: claimed R_lo > kappa_hi (edge route cannot reach C_inf)")
    fl = q(CLAIM[name]).sqrt() / (2 * PI)
    check(bool(fl > q(FLOOR[name])), f"{name}: edge-route floor sqrt(R_lo)/(2pi) > {float(FLOOR[name])}")
    dmax = 1 - CLAIM[name] / 2
    print(f"{name}: every edge reserve satisfies d_* <= 1 - R/2 < 1 - R_claim/2 = {dmax} = {float(dmax):.8f}")
    # E6: two Matsubara rows (n = 0, 1) already defeat every edge reserve at the C_inf level
    R2 = rayleigh_lower(b, TWO_ROW[name])
    print(f"{name}: two-row R ball =", R2.str(25))
    check(bool(R2 > q(CLAIM2[name])), f"{name}: two-row R_lo > {CLAIM2[name]} = {float(CLAIM2[name])}")
    check(bool(q(CLAIM2[name]) > kappa_hi), f"{name}: two-row claim > kappa_hi (rows 0,1 suffice for the obstruction)")

# consistency with FINITE-ENERGY certified reserves (no contradiction claimed or found)
check(CLAIM["cent"] < Fr(41, 30), "cent: R_lo < 41/30 = 2(1-19/60) (consistent with FINITE-ENERGY centered reserve)")
check(CLAIM["caus"] < Fr(45, 32), "caus: R_lo < 45/32 = 2(1-19/64) (consistent with FINITE-ENERGY causal reserve)")
check(CLAIM["cent"] < 2 * (1 - Fr(317, 1000)), "cent: R_lo < 683/500 = 2(1-317/1000) (consistent with certified reserve d_cent,*)")
check(CLAIM["caus"] < 2 * (1 - Fr(119, 400)), "caus: R_lo < 281/200 = 2(1-119/400) (consistent with certified reserve d_caus,*)")

# share of each open bracket [C_inf, C_class] that lies below the edge-route floor
CINF_LO = arb("0.1827262477723295078215179")
C_CLASS = {"cent": q(Fr(41, 30)).sqrt() / (2 * PI), "caus": q(Fr(45, 32)).sqrt() / (2 * PI)}
FRAC = {"cent": Fr(805, 1000), "caus": Fr(911, 1000)}
if MUT == "frac_up":
    FRAC["cent"] = Fr(806, 1000)
for name in ("cent", "caus"):
    print(f"{name}: C_class =", C_CLASS[name].str(20))
    frac = (q(FLOOR[name]) - arb(CINF_HI)) / (C_CLASS[name] - CINF_LO)
    check(bool(frac > q(FRAC[name])), f"{name}: (floor - C_inf)/(C_class - C_inf) > {float(FRAC[name])} (share unreachable by any edge reserve)")
print("ALL CHECKS PASSED")
sys.exit(0)
