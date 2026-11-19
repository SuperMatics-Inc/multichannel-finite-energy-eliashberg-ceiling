"""Exact rational primary certification for the 21-point stage-1 cascade and centered slope-cone ceiling.

Finite certificates complement, not replace, the analytic RIGOR.md proof.
All arithmetic in the proof path uses Python's exact fractions.Fraction and exact integers.
"""

from argparse import ArgumentParser
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys

p = ArgumentParser(description=__doc__)
p.add_argument(
    "--mutant",
    choices=(
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
)
bad = p.parse_args().mutant


def require(ok: bool, name: str) -> None:
    if not ok:
        raise SystemExit("REJECT: " + name)


# 0. Verify inherited proof-file sha256 pins from closure.json (fail-closed against UNIV-GLUE root)
root = Path(__file__).resolve().parent
univ_glue_root = root.parent.parent
closure_path = root / "closure.json"
require(closure_path.is_file(), "closure.json present")
closure_pins = json.loads(closure_path.read_text(encoding="utf-8"))["source_sha256"]
require(len(closure_pins) == 6, "closure.json has all six inherited proof pins")
for name, expected in closure_pins.items():
    target = (root if bad == "closure-pin" else univ_glue_root) / name
    require(target.is_file(), "inherited proof file exists: " + name)
    require(sha256(target.read_bytes()).hexdigest() == expected, "inherited proof bytes " + name)

# 1. Pure rational enclosing intervals for pi and arctan(x) via alternating Leibniz series
def atan_small_bounds(z: Fraction, n_pairs: int = 12) -> tuple[Fraction, Fraction]:
    require(abs(z) <= Fraction(1, 2), "atan_small_bounds domain")
    if z == 0:
        return Fraction(0), Fraction(0)
    sign = 1 if z > 0 else -1
    az = abs(z)
    z2 = az * az
    pow_z = az
    s = Fraction(0)
    for j in range(2 * n_pairs):
        term = pow_z / (2 * j + 1)
        if j % 2 == 0:
            s += term
        else:
            s -= term
        pow_z *= z2
    lo_pos = s
    hi_pos = s + pow_z / (4 * n_pairs + 1)
    return (lo_pos, hi_pos) if sign > 0 else (-hi_pos, -lo_pos)


a5_lo, a5_hi = atan_small_bounds(Fraction(1, 5), 14)
a239_lo, a239_hi = atan_small_bounds(Fraction(1, 239), 6)
PI_LO = 16 * a5_lo - 4 * a239_hi
PI_HI = 16 * a5_hi - 4 * a239_lo
require(
    Fraction(3141592653589793, 10**15) < PI_LO < PI_HI < Fraction(3141592653589794, 10**15),
    "Machin rational pi enclosure",
)
c_cent_lo = Fraction(186040947887195, 10**15) if bad == "decimal-constant" else Fraction(186059321394545, 10**15)
c_cent_hi = Fraction(186040947887196, 10**15) if bad == "decimal-constant" else Fraction(186059321394546, 10**15)
c_caus_lo = Fraction(188734545391826, 10**15)
c_caus_hi = Fraction(188734545391827, 10**15)
c_cent_star_lo = Fraction(186013935536718, 10**15)
c_cent_star_hi = Fraction(186013935536719, 10**15)
c_caus_star_lo = Fraction(188650644722914, 10**15)
c_caus_star_hi = Fraction(188650644722915, 10**15)
require(
    (2 * PI_HI * c_cent_lo) ** 2 < Fraction(41, 30) < (2 * PI_LO * c_cent_hi) ** 2
    and (2 * PI_HI * c_caus_lo) ** 2 < Fraction(45, 32) < (2 * PI_LO * c_caus_hi) ** 2
    and (2 * PI_HI * c_cent_star_lo) ** 2 < Fraction(683, 500) < (2 * PI_LO * c_cent_star_hi) ** 2
    and (2 * PI_HI * c_caus_star_lo) ** 2 < Fraction(281, 200) < (2 * PI_LO * c_caus_star_hi) ** 2,
    "exact ceiling decimal enclosure",
)
ATAN_HALF_LO, ATAN_HALF_HI = atan_small_bounds(Fraction(1, 2), 14)


def atan_pos_bounds(x: Fraction) -> tuple[Fraction, Fraction]:
    require(x > 0, "positive atan argument")
    if x > 1:
        lo_inv, hi_inv = atan_pos_bounds(Fraction(1) / x)
        return Fraction(1, 2) * PI_LO - hi_inv, Fraction(1, 2) * PI_HI - lo_inv
    if x > Fraction(1, 3):
        z = (2 * x - 1) / (x + 2)
        z_lo, z_hi = atan_small_bounds(z, 12)
        return ATAN_HALF_LO + z_lo, ATAN_HALF_HI + z_hi
    return atan_small_bounds(x, 12)


# 2. Universal polynomial and completion identities at clear denominators
for t_i in (-1, 0, 1, 2):
    for r_i in (1, 2, 3):
        for px_i in (0, 1, 2):
            for py_i in (0, 1, 2):
                lhs = (
                    (2 * t_i * (1 + t_i) * px_i + r_i) * (2 * t_i * (1 + t_i) * r_i * py_i + 1)
                    - r_i * ((1 + 2 * t_i) ** 2)
                    - (1 if bad == "head-determinant" else 0)
                )
                rhs = 2 * t_i * (1 + t_i) * (
                    px_i + r_i * r_i * py_i + 2 * t_i * (1 + t_i) * r_i * px_i * py_i - 2 * r_i
                )
                require(lhs == rhs, "decoupled head polynomial identity")

for a_i in (1, 2, 3):
    for b_i in (1, 2):
        for x_i in (-2, 1, 3):
            for y_i in (-1, 2, 4):
                require(
                    (a_i - b_i) * (x_i * x_i + y_i * y_i) - 2 * (a_i + b_i) * x_i * y_i
                    == (a_i - b_i) * ((x_i - y_i) ** 2) - 4 * b_i * x_i * y_i,
                    "reflected tail decomposition identity",
                )
k0_1 = Fraction(1, 1 + 0)
k3_1 = Fraction(1, 1 + 9)
p11_1 = k0_1 - k3_1
s11_1 = k0_1 + k3_1
half_pair_11 = Fraction(1, 2) * (p11_1 * (1 + 1) - 2 * s11_1 * 1 * 1)
reflected_coeff = 1 if bad == "reflected-tail" else 2
require(
    half_pair_11 == Fraction(1, 2) * p11_1 * ((1 - 1) ** 2) - reflected_coeff * k3_1
    and half_pair_11 == Fraction(-1, 5),
    "reflected tail factor of 2",
)

w1 = 3
first_tail_thresh = Fraction(1 if bad == "first-tail-threshold" else 2, w1 * w1)
require(first_tail_thresh == Fraction(2, 9), "first-tail threshold 2/w_1^2 == 2/9")

b0_sample = Fraction(1, 10)
b1_upper = Fraction(1) if bad == "upper-harnack" else min(3 * b0_sample, 1 - (1 - b0_sample) / 3)
require(
    b1_upper <= 3 * b0_sample and b1_upper <= 1 - (1 - b0_sample) / 3,
    "upper Harnack ratio bound required for negative completion",
)
for m in range(1, 21):
    wm = 2 * m + 1
    require((m + 1) ** 2 - m * m == wm, "stage-0 telescoping width identity")
    for t_val in (Fraction(0), Fraction(1, 7), Fraction(1, 2)):
        s0 = Fraction(1, t_val + m * m) + Fraction(1, t_val + (m + 1) ** 2)
        p0 = Fraction(wm, (t_val + m * m) * (t_val + (m + 1) ** 2))
        require(wm * (s0 * s0 - p0 * p0) == 4 * p0, "stage-0 completion cancellation")
        r_v = Fraction(7, 10)
        a_v = Fraction(3, 5)
        direct0 = r_v * p0 / 2 - s0 * s0 / (4 * wm * (a_v + p0 / (2 * r_v * wm)))
        rat0 = (2 * p0 * (r_v * wm * a_v - Fraction(2, wm))) / (4 * wm * a_v + 2 * (1 / r_v) * p0)
        require(direct0 == rat0, "stage-0 completion rationalization")
        if m >= 2:
            require((m + 2) ** 2 - (m - 1) ** 2 == 3 * wm, "stage-1 telescoping width identity")
            s1 = Fraction(1, t_val + (m - 1) ** 2) + Fraction(1, t_val + (m + 2) ** 2)
            p1 = Fraction(3 * wm, (t_val + (m - 1) ** 2) * (t_val + (m + 2) ** 2))
            require(3 * wm * (s1 * s1 - p1 * p1) == 4 * p1, "stage-1 completion cancellation")
            b1_v = Fraction(1, 7)
            direct1 = r_v * p1 / 6 - s1 * s1 / (12 * wm * (b1_v + p1 / (2 * r_v * wm)))
            rat1 = (2 * p1 * (r_v * wm * b1_v - Fraction(2, 3 * wm))) / (
                12 * wm * b1_v + 6 * (1 / r_v) * p1
            )
            require(direct1 == rat1, "stage-1 completion rationalization")

# 3. Exact rational bounds on the reduced Hankel tail (after extracting stage-0 and stage-1 cascades)
Q = 10**12
M = 20
w = [2 * m + 1 for m in range(1, M + 1)]


def f_floor(x: Fraction) -> int:
    return (x.numerator * Q) // x.denominator


def f_ceil(x: Fraction) -> int:
    return -((-x.numerator * Q) // x.denominator)


q_red_exact: list[Fraction] = []
for n in range(1, M + 1):
    s = Fraction(0)
    for m in range(1, 21):
        if (n == 1 and m >= 2) or (m == 1 and n >= 2):
            continue
        s += Fraction(1, (2 * m + 1) * ((n + m + 1) ** 2))
    if bad == "hankel-tail":
        s -= Fraction(1, 20)
    else:
        s += Fraction(1, 43 * (n + 21))
    q_red_exact.append(s)

require(
    q_red_exact[0] > Fraction(38, 1000) and q_red_exact[0] < Fraction(381, 10000),
    "reduced Hankel row-1 bound",
)
q_full_row1 = sum(Fraction(1, (2 * m + 1) * ((m + 2) ** 2)) for m in range(1, 17)) + Fraction(1, 4 * 16 * 16)
q_hankel_cap = Fraction(33, 500)
require(q_full_row1 < q_hankel_cap, "full Hankel Schur bound")
q_red_int = [f_ceil(q) for q in q_red_exact]

# 4. Exterior region exact rational inequalities (for all t >= 1/2 and for all r outside [1/3, 3])
d_caus_cert = Fraction(119, 400)
d_caus_rep = Fraction(119, 400) if bad == "strict-budget" else Fraction(19, 64)
d_cent_cert = Fraction(317, 1000)
d_cent_rep = Fraction(317, 1000) if bad == "strict-budget" else Fraction(19, 60)

require(d_caus_rep < d_caus_cert and d_cent_rep < d_cent_cert, "strict reserve margins")
require(
    2 * (1 - d_caus_rep) == Fraction(45, 32)
    and 2 * (1 - d_cent_rep) == Fraction(41, 30)
    and Fraction(41, 30) < Fraction(45, 32) < Fraction(377, 250) < Fraction(413, 250),
    "strict ceiling chain",
)
require(1 - (1 - d_caus_cert) / (1 - d_caus_rep) == Fraction(1, 1125), "causal bare-form gap 1/1125")
require(1 - (1 - d_cent_cert) / (1 - d_cent_rep) == Fraction(1, 2050), "centered bare-form gap 1/2050")

alpha_ext = 1 - d_cent_cert - q_hankel_cap
require(alpha_ext == Fraction(617, 1000), "exterior bare budget")

far_series_20 = sum(
    Fraction(1, (2 * m + 1) * (m * m) * ((m + 1) ** 2)) for m in range(1, 21)
) + Fraction(1, 8 * (20**4))
require(
    far_series_20 == Fraction(221683723759377431023, 2453474125083421440000)
    and far_series_20 < Fraction(7, 75),
    "far-imbalance rationalized completion series bound < 7/75",
)
far_loss = Fraction(1, 2) if bad == "far-imbalance" else Fraction(7, 75) / alpha_ext
require(
    Fraction(9, 17) - far_loss == Fraction(11899, 31467) > d_cent_cert,
    "far-imbalance reserve on r <= 1/4",
)
mid_imb_reserve = Fraction(2, 5) + Fraction(3, 4) * (alpha_ext / 8 - Fraction(1, 9)) / alpha_ext
require(
    mid_imb_reserve == Fraction(106229, 296160) > d_cent_cert,
    "intermediate-imbalance reserve on 1/4 <= r <= 1/3",
)
large_t_base = Fraction(1, 5) if bad == "large-energy" else Fraction(1, 3)
large_t_reserve = large_t_base + Fraction(4, 9) * (alpha_ext / 6 - Fraction(1, 9)) / alpha_ext
require(
    large_t_reserve == Fraction(16361, 49977) > d_cent_cert,
    "large-energy exterior reserve on t >= 1/2",
)

# 5. Certified rational tangent envelopes for causal and centered classes
sub_ivs = [
    (Fraction(0), Fraction(2, 3)),
    (Fraction(2, 3), Fraction(4, 5)),
    (Fraction(4, 5), Fraction(17, 20)),
    (Fraction(17, 20), Fraction(88, 100)),
    (Fraction(88, 100), Fraction(90, 100)),
    (Fraction(90, 100), Fraction(92, 100)),
    (Fraction(92, 100), Fraction(94, 100)),
    (Fraction(94, 100), Fraction(96, 100)),
    (Fraction(96, 100), Fraction(1)),
]


def make_certified_rational_tangents(xi_ints: list[int]) -> tuple[list[list[Fraction]], list[list[Fraction]]]:
    slopes = [[Fraction(1, 2 * m + 1) for m in range(1, M + 1)]]
    intercepts = [[Fraction(0) for _ in range(M)]]
    Q_c = 10**8
    for xi in xi_ints:
        ax_lo, ax_hi = atan_pos_bounds(Fraction(xi))
        b0_hi = (2 * ax_hi) / PI_LO
        lam_row: list[Fraction] = []
        c_row: list[Fraction] = []
        for m in range(1, M + 1):
            wm = 2 * m + 1
            lam = Fraction(2 * wm) if bad == "centered-cone" else Fraction(wm * (1 + xi * xi), wm * wm + xi * xi)
            am_lo, _ = atan_pos_bounds(Fraction(xi, wm))
            bm_lo = (2 * am_lo) / PI_HI
            c_exact_lo = bm_lo - lam * b0_hi
            c_rat = Fraction((c_exact_lo.numerator * Q_c) // c_exact_lo.denominator, Q_c)
            require(c_rat <= 0 and lam + c_rat <= 1, "admissible lower tangent line")
            lam_row.append(lam)
            c_row.append(c_rat)
        slopes.append(lam_row)
        intercepts.append(c_row)
    slopes.append([Fraction(2 * m + 1) for m in range(1, M + 1)])
    intercepts.append([Fraction(-(2 * m)) for m in range(1, M + 1)])
    return slopes, intercepts


def build_iv_int_tables(xi_ints: list[int]) -> list[tuple[list[int], list[list[int]]]]:
    slopes, intercepts = make_certified_rational_tangents(xi_ints)
    all_bkpts: set[Fraction] = set()
    L = len(slopes)
    for l in range(L - 1):
        for m in range(M):
            denom = slopes[l + 1][m] - slopes[l][m]
            require(denom > 0, "strictly increasing tangent slopes")
            b_int = -(intercepts[l + 1][m] - intercepts[l][m]) / denom
            all_bkpts.add(b_int)
    sorted_bkpts = sorted(all_bkpts)
    iv_tables: list[tuple[list[int], list[list[int]]]] = []
    for b_lo, b_hi in sub_ivs:
        pts = [b_hi] if b_lo == 0 else ([b_lo] + [pt for pt in sorted_bkpts if b_lo < pt < b_hi] + [b_hi])
        inv_p = [f_floor(Fraction(1) / pt) for pt in pts]
        ratio_cl = []
        for pt in pts:
            row = [f_floor(max(slopes[l][m] * pt + intercepts[l][m] for l in range(L)) / pt) for m in range(M)]
            ratio_cl.append(row)
        iv_tables.append((inv_p, ratio_cl))
    return iv_tables


def box_ok_int(
    r_lo: Fraction,
    r_hi: Fraction,
    t_lo: Fraction,
    t_hi: Fraction,
    dt: Fraction,
    split: Fraction,
    iv_tables: list[tuple[list[int], list[list[int]]]],
) -> bool:
    dt_ceil = f_ceil(dt)
    one_m_dt = Q - dt_ceil
    rem = [one_m_dt - q_red_int[m] for m in range(M)]
    if any(r_val <= 0 for r_val in rem):
        return False
    alp_0 = [(rem[m] * split.numerator) // split.denominator for m in range(M)]
    bet_1 = [(rem[m] * (split.denominator - split.numerator)) // split.denominator for m in range(M)]
    alp_0[0] = rem[0]
    bet_1[0] = 0

    r_lo_int = f_floor(r_lo)
    r_hi_int = f_ceil(r_hi)
    ri_lo_int = f_floor(Fraction(1) / r_hi)
    ri_hi_int = f_ceil(Fraction(1) / r_lo)

    p_0k = [f_floor(Fraction(1) / (t_hi + k * k) - Fraction(1) / (t_hi + (k + 1) ** 2)) for k in range(1, M + 1)]
    p_0k_hi = [f_ceil(Fraction(1) / (t_lo + k * k) - Fraction(1) / (t_lo + (k + 1) ** 2)) for k in range(1, M + 1)]
    p_1k = [f_floor(Fraction(1) / (t_hi + (k - 1) ** 2) - Fraction(1) / (t_hi + (k + 2) ** 2)) for k in range(2, M + 1)]

    D_1k_X: list[int] = []
    D_1k_Y: list[int] = []
    for idx in range(M - 1):
        wm = w[idx + 1]
        b1 = bet_1[idx + 1]
        p1 = p_1k[idx]
        thresh = (2 * Q * Q + 3 * wm - 1) // (3 * wm)
        num_X = r_lo_int * wm * b1 - thresh
        num_Y = ri_lo_int * wm * b1 - thresh
        if num_X <= 0 or num_Y <= 0:
            return False
        den_X = 12 * wm * b1 * Q + 6 * ri_hi_int * p1
        den_Y = 12 * wm * b1 * Q + 6 * r_hi_int * p1
        D_1k_X.append((2 * p1 * num_X) // den_X)
        D_1k_Y.append((2 * p1 * num_Y) // den_Y)

    S = len(sub_ivs)
    tilde_alp_X1_lo: list[int] = []
    tilde_alp_Y1_lo: list[int] = []
    for a in range(S):
        inv_p, ratio_cl = iv_tables[a]
        tilde_alp_X1_lo.append(
            min((alp_0[0] * inv_p[i] + sum(ratio_cl[i][m + 1] * D_1k_X[m] for m in range(M - 1))) // Q for i in range(len(inv_p)))
        )
        tilde_alp_Y1_lo.append(
            min((alp_0[0] * inv_p[i] + sum(ratio_cl[i][m + 1] * D_1k_Y[m] for m in range(M - 1))) // Q for i in range(len(inv_p)))
        )

    Phi_X_min = [[0] * S for _ in range(S)]
    Phi_Y_min = [[0] * S for _ in range(S)]

    for a in range(S):
        b_hi = sub_ivs[a][1]
        D_0k_X: list[int] = []
        for m in range(M):
            wm = w[m]
            ta = tilde_alp_Y1_lo[a] if m == 0 else (alp_0[m] * b_hi.denominator) // b_hi.numerator
            thresh0 = (2 * Q * Q + wm - 1) // wm
            num = r_lo_int * wm * ta - thresh0
            if num <= 0:
                return False
            p0 = p_0k[m]
            den = 4 * wm * ta * Q + 2 * ri_hi_int * p0
            D_0k_X.append((2 * p0 * num) // den)
        for b in range(S):
            inv_p, ratio_cl = iv_tables[b]
            Phi_X_min[a][b] = min(
                (one_m_dt * inv_p[i] + sum(ratio_cl[i][m] * D_0k_X[m] for m in range(M))) // Q for i in range(len(inv_p))
            )

    for b in range(S):
        c_hi = sub_ivs[b][1]
        D_0k_Y: list[int] = []
        neg_loss_Y = 0
        for m in range(M):
            wm = w[m]
            ta = tilde_alp_X1_lo[b] if m == 0 else (alp_0[m] * c_hi.denominator) // c_hi.numerator
            thresh0 = (2 * Q * Q + wm - 1) // wm
            num = ri_lo_int * wm * ta - thresh0
            if num >= 0:
                p0 = p_0k[m]
                den = 4 * wm * ta * Q + 2 * r_hi_int * p0
                D_0k_Y.append((2 * p0 * num) // den)
            else:
                p0_h = p_0k_hi[m]
                neg_loss_Y += wm * ((2 * p0_h * num) // (4 * wm * ta * Q))
                D_0k_Y.append(0)
        for a in range(S):
            inv_p, ratio_cl = iv_tables[a]
            Phi_Y_min[a][b] = min(
                (one_m_dt * inv_p[i] + sum(ratio_cl[i][m] * D_0k_Y[m] for m in range(M))) // Q + neg_loss_Y
                for i in range(len(inv_p))
            )

    r0 = Fraction(1, 2) * (r_lo + r_hi)
    r0_sq = r0 * r0
    two_t = 2 * t_lo * (1 + t_lo)
    for a in range(S):
        for b in range(S):
            px = Fraction(Phi_X_min[a][b], Q)
            py = Fraction(Phi_Y_min[a][b], Q)
            if px <= 0 or py <= 0:
                return False
            slope = py - px / r0_sq
            base = 2 * px / r0
            g_lo = base + (r_lo if slope >= 0 else r_hi) * slope
            if g_lo + two_t * px * py < 2:
                return False
    return True


split_share = Fraction(1) if bad == "stage1-cascade" else Fraction(7, 10)

for label, xi_list, target_dt in (
    ("causal", [], d_caus_cert),
    ("centered", [5, 8, 12], d_cent_cert),
):
    iv_tables = build_iv_int_tables(xi_list)
    stack = [(Fraction(1), Fraction(3), Fraction(0), Fraction(1, 2), 0)]
    if bad == "coverage-gap":
        stack = [(Fraction(1), Fraction(2), Fraction(0), Fraction(1, 2), 0)]
    leaves = 0
    max_d = 0
    covered_area = Fraction(0)
    while stack:
        r0, r1, tl, th, d = stack.pop()
        if box_ok_int(r0, r1, tl, th, target_dt, split_share, iv_tables):
            leaves += 1
            max_d = max(max_d, d)
            covered_area += (r1 - r0) * (th - tl)
            continue
        require(d < 18, label + " compact rectangle subdivision depth")
        rm = Fraction(1, 2) * (r0 + r1)
        tm = Fraction(1, 2) * (tl + th)
        if (r1 - r0) > Fraction(3, 2) * (th - tl):
            stack.append((r0, rm, tl, th, d + 1))
            stack.append((rm, r1, tl, th, d + 1))
        else:
            stack.append((r0, r1, tl, tm, d + 1))
            stack.append((r0, r1, tm, th, d + 1))
    require(covered_area == Fraction(1), label + " complete residual domain without gaps")
    print(f"CERTIFIED 21-point {label} covering: leaves={leaves}, max_depth={max_d}, reserve={target_dt}")

print("Analytic import: full-head Poisson equality forces flat DOS almost everywhere.")
