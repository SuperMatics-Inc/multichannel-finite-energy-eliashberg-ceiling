"""Independent 4-way full-domain quadtree covering and actual-weight finite matrix PSD verification.

No import of the primary replay; all arithmetic is exact fractions.Fraction and exact integers.
"""

from argparse import ArgumentParser
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

p = ArgumentParser(description=__doc__)
p.add_argument(
    "--mutant",
    choices=(
        "quadtree-area",
        "completion-sign",
        "matrix-gain",
        "missing-normal-mass",
        "tangent-intercept",
        "flat-admissibility",
        "psd-schur",
        "two-sided-bathtub",
    ),
)
bad = p.parse_args().mutant


def demand(ok: bool, name: str) -> None:
    if not ok:
        raise SystemExit("REJECT: " + name)


root = Path(__file__).resolve().parent
univ_glue_root = root.parent.parent
for rel_name, expected in json.loads((root / "closure.json").read_text(encoding="utf-8"))["source_sha256"].items():
    target = univ_glue_root / rel_name
    demand(sha256(target.read_bytes()).hexdigest() == expected, "inherited proof bytes " + rel_name)


# 1. Independent rational atan and pi bounds (12-pair / 8-pair Leibniz enclosure)
def atan_leibniz(z: Fraction, pairs: int) -> tuple[Fraction, Fraction]:
    demand(0 <= z <= Fraction(1, 2), "atan_leibniz nonneg domain")
    if z == 0:
        return Fraction(0), Fraction(0)
    z2 = z * z
    cur = z
    acc = Fraction(0)
    for k in range(2 * pairs):
        term = cur / (2 * k + 1)
        acc = acc + term if (k % 2 == 0) else acc - term
        cur *= z2
    return acc, acc + cur / (4 * pairs + 1)


a5_l, a5_h = atan_leibniz(Fraction(1, 5), 12)
a239_l, a239_h = atan_leibniz(Fraction(1, 239), 8)
PI_L = 16 * a5_l - 4 * a239_h
PI_H = 16 * a5_h - 4 * a239_l
AH_L, AH_H = atan_leibniz(Fraction(1, 2), 16)


def atan_pos(x: Fraction) -> tuple[Fraction, Fraction]:
    demand(x > 0, "atan_pos argument")
    if x > 1:
        li, hi = atan_pos(Fraction(1) / x)
        return Fraction(1, 2) * PI_L - hi, Fraction(1, 2) * PI_H - li
    if x > Fraction(1, 3):
        u = (2 * x - 1) / (x + 2)
        lu, hu = atan_leibniz(abs(u), 14)
        return (AH_L + lu, AH_H + hu) if u >= 0 else (AH_L - hu, AH_H - lu)
    return atan_leibniz(x, 14)


for m in range(1, 21):
    wm = 2 * m + 1
    r_v = Fraction(7, 10)
    t_v = Fraction(1, 11)
    a_v = Fraction(3, 5)
    p0 = Fraction(wm, (t_v + m * m) * (t_v + (m + 1) ** 2))
    s0 = Fraction(1, t_v + m * m) + Fraction(1, t_v + (m + 1) ** 2)
    sgn = -1 if bad == "completion-sign" else 1
    direct0 = r_v * p0 / 2 - s0 * s0 / (4 * wm * (a_v + p0 / (2 * r_v * wm)))
    rat0 = (2 * p0 * (r_v * wm * a_v - Fraction(2, wm))) / (4 * wm * a_v + sgn * 2 * (1 / r_v) * p0)
    demand(direct0 == rat0, "independent stage-0 completion rationalization")
    if m >= 2:
        p1 = Fraction(3 * wm, (t_v + (m - 1) ** 2) * (t_v + (m + 2) ** 2))
        s1 = Fraction(1, t_v + (m - 1) ** 2) + Fraction(1, t_v + (m + 2) ** 2)
        b1_v = Fraction(1, 7)
        direct1 = r_v * p1 / 6 - s1 * s1 / (12 * wm * (b1_v + p1 / (2 * r_v * wm)))
        rat1 = (2 * p1 * (r_v * wm * b1_v - Fraction(2, 3 * wm))) / (12 * wm * b1_v + sgn * 6 * (1 / r_v) * p1)
        demand(direct1 == rat1, "independent stage-1 completion rationalization")

# 2. Independent 4-way quadtree covering of the FULL [1/3, 3] x [0, 1/2] domain (area = 4/3)
Q = 10**12
M = 20
w = [2 * m + 1 for m in range(1, M + 1)]


def q_floor(x: Fraction) -> int:
    return (x.numerator * Q) // x.denominator


def q_ceil(x: Fraction) -> int:
    return -((-x.numerator * Q) // x.denominator)


q_red: list[int] = []
for n in range(1, M + 1):
    s = Fraction(1, 43 * (n + 21))
    for m in range(1, 21):
        if (n == 1 and m >= 2) or (m == 1 and n >= 2):
            continue
        s += Fraction(1, (2 * m + 1) * ((n + m + 1) ** 2))
    q_red.append(q_ceil(s))

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


def make_tables(xi_list: list[int]) -> list[tuple[list[int], list[list[int]]]]:
    slopes = [[Fraction(1, 2 * m + 1) for m in range(1, M + 1)]]
    intercepts = [[Fraction(0) for _ in range(M)]]
    Q_c = 10**8
    for xi in xi_list:
        _, ax_hi = atan_pos(Fraction(xi))
        b0_hi = (2 * ax_hi) / PI_L
        l_row: list[Fraction] = []
        c_row: list[Fraction] = []
        for m in range(1, M + 1):
            wm = 2 * m + 1
            lam = Fraction(wm * (1 + xi * xi), wm * wm + xi * xi)
            am_lo, _ = atan_pos(Fraction(xi, wm))
            bm_lo = (2 * am_lo) / PI_H
            c_lo = bm_lo - lam * b0_hi
            c_rat = Fraction(1, 10) if bad == "tangent-intercept" else Fraction((c_lo.numerator * Q_c) // c_lo.denominator, Q_c)
            demand(c_rat <= 0 and lam + c_rat <= 1, "certified tangent intercept nonpositive")
            l_row.append(lam)
            c_row.append(c_rat)
        slopes.append(l_row)
        intercepts.append(c_row)
    slopes.append([Fraction(2 * m + 1) for m in range(1, M + 1)])
    intercepts.append([Fraction(-(2 * m)) for m in range(1, M + 1)])

    bkpts: set[Fraction] = set()
    L = len(slopes)
    for l in range(L - 1):
        for m in range(M):
            bkpts.add(-(intercepts[l + 1][m] - intercepts[l][m]) / (slopes[l + 1][m] - slopes[l][m]))
    s_bkpts = sorted(bkpts)
    tables: list[tuple[list[int], list[list[int]]]] = []
    for b_lo, b_hi in sub_ivs:
        pts = [b_hi] if b_lo == 0 else ([b_lo] + [pt for pt in s_bkpts if b_lo < pt < b_hi] + [b_hi])
        inv_p = [q_floor(Fraction(1) / pt) for pt in pts]
        ratio_cl = [[q_floor(max(slopes[l][m] * pt + intercepts[l][m] for l in range(L)) / pt) for m in range(M)] for pt in pts]
        tables.append((inv_p, ratio_cl))
    return tables


def cell_ok_full(
    r_lo: Fraction,
    r_hi: Fraction,
    t_lo: Fraction,
    t_hi: Fraction,
    dt: Fraction,
    split: Fraction,
    tables: list[tuple[list[int], list[list[int]]]],
) -> bool:
    one_m_dt = Q - q_ceil(dt)
    rem = [one_m_dt - q_red[m] for m in range(M)]
    alp_0 = [(rem[m] * split.numerator) // split.denominator for m in range(M)]
    bet_1 = [(rem[m] * (split.denominator - split.numerator)) // split.denominator for m in range(M)]
    alp_0[0] = rem[0]
    bet_1[0] = 0

    r_lo_i = q_floor(r_lo)
    r_hi_i = q_ceil(r_hi)
    ri_lo_i = q_floor(Fraction(1) / r_hi)
    ri_hi_i = q_ceil(Fraction(1) / r_lo)

    p_0k = [q_floor(Fraction(1) / (t_hi + k * k) - Fraction(1) / (t_hi + (k + 1) ** 2)) for k in range(1, M + 1)]
    p_0k_hi = [q_ceil(Fraction(1) / (t_lo + k * k) - Fraction(1) / (t_lo + (k + 1) ** 2)) for k in range(1, M + 1)]
    p_1k = [q_floor(Fraction(1) / (t_hi + (k - 1) ** 2) - Fraction(1) / (t_hi + (k + 2) ** 2)) for k in range(2, M + 1)]

    D_1k_X: list[int] = []
    D_1k_Y: list[int] = []
    for idx in range(M - 1):
        wm = w[idx + 1]
        b1 = bet_1[idx + 1]
        p1 = p_1k[idx]
        thresh = (2 * Q * Q + 3 * wm - 1) // (3 * wm)
        num_X = r_lo_i * wm * b1 - thresh
        num_Y = ri_lo_i * wm * b1 - thresh
        if num_X <= 0 or num_Y <= 0:
            return False
        D_1k_X.append((2 * p1 * num_X) // (12 * wm * b1 * Q + 6 * ri_hi_i * p1))
        D_1k_Y.append((2 * p1 * num_Y) // (12 * wm * b1 * Q + 6 * r_hi_i * p1))

    S = len(sub_ivs)
    t_alp_X1 = [
        min((alp_0[0] * tables[a][0][i] + sum(tables[a][1][i][m + 1] * D_1k_X[m] for m in range(M - 1))) // Q for i in range(len(tables[a][0])))
        for a in range(S)
    ]
    t_alp_Y1 = [
        min((alp_0[0] * tables[a][0][i] + sum(tables[a][1][i][m + 1] * D_1k_Y[m] for m in range(M - 1))) // Q for i in range(len(tables[a][0])))
        for a in range(S)
    ]

    Phi_X = [[0] * S for _ in range(S)]
    Phi_Y = [[0] * S for _ in range(S)]
    sign_den = -1 if bad == "completion-sign" else 1

    for a in range(S):
        b_hi = sub_ivs[a][1]
        D_0k_X: list[int] = []
        neg_X = 0
        for m in range(M):
            wm = w[m]
            ta = t_alp_Y1[a] if m == 0 else (alp_0[m] * b_hi.denominator) // b_hi.numerator
            thresh0 = (2 * Q * Q + wm - 1) // wm
            num = r_lo_i * wm * ta - thresh0
            if num >= 0:
                p0 = p_0k[m]
                den = 4 * wm * ta * Q + sign_den * 2 * ri_hi_i * p0
                demand(den > 0, "positive stage-0 completion denominator")
                D_0k_X.append((2 * p0 * num) // den)
            else:
                p0_h = p_0k_hi[m]
                neg_X += wm * ((2 * p0_h * num) // (4 * wm * ta * Q))
                D_0k_X.append(0)
        for b in range(S):
            inv_p, ratio_cl = tables[b]
            Phi_X[a][b] = min(
                (one_m_dt * inv_p[i] + sum(ratio_cl[i][m] * D_0k_X[m] for m in range(M))) // Q + neg_X for i in range(len(inv_p))
            )

    for b in range(S):
        c_hi = sub_ivs[b][1]
        D_0k_Y: list[int] = []
        neg_Y = 0
        for m in range(M):
            wm = w[m]
            ta = t_alp_X1[b] if m == 0 else (alp_0[m] * c_hi.denominator) // c_hi.numerator
            thresh0 = (2 * Q * Q + wm - 1) // wm
            num = ri_lo_i * wm * ta - thresh0
            if num >= 0:
                p0 = p_0k[m]
                den = 4 * wm * ta * Q + sign_den * 2 * r_hi_i * p0
                demand(den > 0, "positive stage-0 completion denominator")
                D_0k_Y.append((2 * p0 * num) // den)
            else:
                p0_h = p_0k_hi[m]
                neg_Y += wm * ((2 * p0_h * num) // (4 * wm * ta * Q))
                D_0k_Y.append(0)
        for a in range(S):
            inv_p, ratio_cl = tables[a]
            Phi_Y[a][b] = min(
                (one_m_dt * inv_p[i] + sum(ratio_cl[i][m] * D_0k_Y[m] for m in range(M))) // Q + neg_Y for i in range(len(inv_p))
            )

    r0 = Fraction(1, 2) * (r_lo + r_hi)
    r0_sq = r0 * r0
    two_t = 2 * t_lo * (1 + t_lo)
    for a in range(S):
        for b in range(S):
            px = Fraction(Phi_X[a][b], Q)
            py = Fraction(Phi_Y[a][b], Q)
            if px <= 0 or py <= 0:
                return False
            slope = py - px / r0_sq
            g_lo = 2 * px / r0 + (r_lo if slope >= 0 else r_hi) * slope
            if g_lo + two_t * px * py < 2:
                return False
    return True


for label, xi_list, target_dt in (
    ("causal", [], Fraction(119, 400)),
    ("centered", [5, 8, 12], Fraction(317, 1000)),
):
    tables = make_tables(xi_list)
    queue = [(Fraction(1, 3), Fraction(1), Fraction(0), Fraction(1, 2), 0), (Fraction(1), Fraction(3), Fraction(0), Fraction(1, 2), 0)]
    area = Fraction(0)
    leaves = 0
    max_d = 0
    while queue:
        rl, rh, tl, th, d = queue.pop()
        if cell_ok_full(rl, rh, tl, th, target_dt, Fraction(7, 10), tables):
            area += (rh - rl) * (th - tl)
            leaves += 1
            max_d = max(max_d, d)
            continue
        demand(d < 12, label + " 4-way quadtree closes")
        rm = Fraction(1, 2) * (rl + rh)
        tm = Fraction(1, 2) * (tl + th)
        quads = [
            (rl, rm, tl, tm, d + 1),
            (rm, rh, tl, tm, d + 1),
            (rl, rm, tm, th, d + 1),
        ]
        if bad != "quadtree-area":
            quads.append((rm, rh, tm, th, d + 1))
        queue.extend(quads)
    demand(area == Fraction(4, 3), label + " full [1/3, 3] x [0, 1/2] area 4/3")
    print(f"CERTIFIED independent 4-way {label} cover: leaves={leaves}, max_depth={max_d}, area={area}")


# 3. Actual-weight 6x6 trial forms with 21 retained normal weights on each side
def psd(mat: list[list[Fraction]]) -> bool:
    n = len(mat)
    a = [row[:] for row in mat]
    for k in range(n):
        piv = a[k][k]
        if piv < 0:
            return False
        if piv == 0:
            if any(a[k][j] != 0 for j in range(k + 1, n)):
                return False
            continue
        for i in range(k + 1, n):
            fac = (a[i][k] if bad == "psd-schur" else a[k][i]) / piv
            for j in range(i, n):
                a[i][j] -= fac * a[k][j]
                if bad != "psd-schur":
                    a[j][i] = a[i][j]
    return True


indef_6x6 = [[Fraction(0) for _ in range(6)] for _ in range(6)]
for idx in range(6):
    indef_6x6[idx][idx] = Fraction(3, 2) if idx in (3, 4) else Fraction(1)
indef_6x6[0][3] = indef_6x6[3][0] = Fraction(1)
indef_6x6[0][4] = indef_6x6[4][0] = Fraction(1)
v_test = [Fraction(-2), Fraction(0), Fraction(0), Fraction(1), Fraction(1), Fraction(0)]
quad_val = sum(v_test[i] * indef_6x6[i][j] * v_test[j] for i in range(6) for j in range(6))
demand(quad_val == Fraction(-1) and not psd(indef_6x6), "PSD Schur complement rejects indefinite 6x6")

a2_lo, a2_hi = atan_pos(Fraction(2))
a1_lo, a1_hi = atan_pos(Fraction(1))
a23_lo, a23_hi = atan_pos(Fraction(2, 3))
b0_two_sided_lo = (2 * a2_lo) / PI_H
bm_sigma_half_lo = (2 * a1_lo) / PI_H
bm_sigma_half_hi = (2 * a1_hi) / PI_L
phi_w1_hi = (2 * a23_hi) / PI_L
claimed_lower = b0_two_sided_lo if bad == "two-sided-bathtub" else phi_w1_hi
demand(
    claimed_lower <= bm_sigma_half_lo and bm_sigma_half_hi < b0_two_sided_lo,
    "two-sided width cone [1/w_m, w_m] bathtub lower bound",
)

c_slopes = [[Fraction(1, 2 * m + 1) for m in range(1, M + 1)]]
c_inter = [[Fraction(0) for _ in range(M)]]
Q_c = 10**6
for xi_val in (5, 8, 12):
    _, ax_hi = atan_pos(Fraction(xi_val))
    b0_hi = (2 * ax_hi) / PI_L
    l_r = []
    c_r = []
    for m in range(1, M + 1):
        wm = 2 * m + 1
        lam = Fraction(wm * (1 + xi_val * xi_val), wm * wm + xi_val * xi_val)
        am_lo, _ = atan_pos(Fraction(xi_val, wm))
        bm_lo = (2 * am_lo) / PI_H
        c_lo = bm_lo - lam * b0_hi
        c_r.append(Fraction((c_lo.numerator * Q_c) // c_lo.denominator, Q_c))
        l_r.append(lam)
    c_slopes.append(l_r)
    c_inter.append(c_r)
c_slopes.append([Fraction(2 * m + 1) for m in range(1, M + 1)])
c_inter.append([Fraction(-(2 * m)) for m in range(1, M + 1)])

count = 0
d_cert = Fraction(319, 1000)
for r in (Fraction(1, 4), Fraction(1, 3), Fraction(7, 10), Fraction(1), Fraction(3), Fraction(4)):
    for t in (Fraction(1, 1000), Fraction(1, 10), Fraction(1)):
        for code in range(4):
            b_sides: list[list[Fraction]] = []
            for side in range(2):
                head = Fraction((code + side + 5) % 9, 8)
                prof = [head]
                for m in range(1, 21):
                    wm = 2 * m + 1
                    low = max(c_slopes[l][m - 1] * head + c_inter[l][m - 1] for l in range(len(c_slopes)))
                    high = min(wm * head, Fraction(wm - 1 + head, wm), Fraction(1))
                    prof.append(low if code % 2 == 0 else (low + high) / 2)
                b_sides.append(prof)
            mat = [[Fraction(0) for _ in range(6)] for _ in range(6)]
            for side in range(2):
                scale = 1 / r if side == 0 else r
                for n in range(3):
                    i = side * 3 + n
                    mat[i][i] = (1 - d_cert) * scale * (2 * n + 1) * b_sides[side][n]
                    for m in range(21):
                        if bad == "missing-normal-mass" and m > 0:
                            continue
                        p_nm = Fraction(1, t + (n - m) ** 2) - Fraction(1, t + (n + m + 1) ** 2)
                        mat[i][i] += p_nm * b_sides[side][n] * b_sides[1 - side][m] / 2
            for n in range(3):
                for m in range(3):
                    s_nm = Fraction(1, t + (n - m) ** 2) + Fraction(1, t + (n + m + 1) ** 2)
                    gain = s_nm * b_sides[0][n] * b_sides[1][m] / (1 if bad == "matrix-gain" else 2)
                    mat[n][3 + m] -= gain
                    mat[3 + m][n] -= gain
            demand(psd(mat), "actual-weight finite form with 21 retained normal weights")
            count += 1

head_w = Fraction(1)
second_w = Fraction(1, 2) if bad == "flat-admissibility" else Fraction(1)
demand(second_w >= 3 * head_w - 2, "near-tight full-head Poisson admissibility")
print(f"CERTIFIED supplementary 6x6 trial forms (21 retained normal weights): count={count}")
print("Analytic import: full-head Poisson equality forces flat DOS almost everywhere.")
