# Two-stage frequency-cascade and centered slope-cone bounds for finite-energy multichannel Eliashberg branches

## Verdict and scope

**PROVED, with computer-assisted finite rational inequalities:** on the reciprocal reduced Eliashberg model at arbitrary finite glue energies,

1. **Centered finite-channel class (`C_cent`):** for every centered even-DOS branch solving the matched positive normal equation (or lying in the two-sided centered width cone $w_m^{-1} \le s_{m,x}/s_{0,x} \le w_m$),
   \[
   k_B T_{\rm lin} < C_{\rm cent}\sqrt{\rho},
   \qquad C_{\rm cent} = \frac{\sqrt{41/30}}{2\pi}.
   \]
   **NUMERICAL / CERTIFIED enclosure:** $C_{\rm cent} \approx 0.18605932139454540$ (certified in exact rational arithmetic in `replay.py` to lie in $[186059321394545/10^{15},\ 186059321394546/10^{15}]$ using Machin's formula with alternating Leibniz remainder bounds). Compared to the prior sharp scalar constant $C_\infty = 0.18272624777232952$ and the preceding multipoint bound $C_{21} = \sqrt{377/250}/(2\pi) \approx 0.19544330695411738$, this reduces the relative excess $(C - C_\infty)/C_\infty$ from $6.960\%$ to $1.824\%$ and closes $73.79\%$ of the open bracket $[C_\infty, C_{21}]$ ($89.35\%$ of the published centered bracket $[C_\infty, U]$, where $U = \sqrt{6466219/3576050}/(2\pi) \approx 0.2140146706$).

2. **General causal finite-channel class (`C_caus`):** for every causal branch obeying the four-inequality Harnack weight law $0 \le b_{m,x} \le 1$, $b_{0,x}/w_m \le b_{m,x} \le w_m b_{0,x}$, and $w_m b_{0,x} - (w_m - 1) \le b_{m,x} \le (b_{0,x} + w_m - 1)/w_m$ (including shifted branches with a channel-wise causal radial propagator),
   \[
   k_B T_{\rm lin} < C_{\rm caus}\sqrt{\rho},
   \qquad C_{\rm caus} = \frac{\sqrt{45/32}}{2\pi}.
   \]
   **NUMERICAL / CERTIFIED enclosure:** $C_{\rm caus} \approx 0.18873454539182638$ (certified in exact rational arithmetic in `replay.py` to lie in $[188734545391826/10^{15},\ 188734545391827/10^{15}]$), reducing the general causal excess over $C_\infty$ from $6.960\%$ to $3.288\%$ and closing $52.75\%$ of $[C_\infty, C_{21}]$.

3. **Single-channel scalar finite-energy case (settled by prior proof):** for a single channel ($d = 1$) at arbitrary finite glue energies and arbitrary bounded energy-dependent DOS $0 \le N(\xi) \le N_{\max}$, the sharp bound $k_B T_{\rm lin} < C_\infty \sqrt{\rho}$ is already **PROVED** in the scalar finite-band paper (`PUBLISH/finite-repo`, Zenodo 24364879). For coupled channels ($d \ge 2$) with energy-independent (flat) DOS, $k_B T_{\rm lin} \le C_\infty \sqrt{\rho}$ at arbitrary glue energies is already **PROVED** in `../../note-centered/PHYS-Q2/RIGOR.md`.

**CERTIFIED:** two independent exact-rational coverings (`replay.py` using binary subdivision on $[1, 3] \times [0, 1/2]$ and `independent.py` using 4-way quadtree subdivision on the full $[1/3, 3] \times [0, 1/2]$ domain without $r \leftrightarrow r^{-1}$ symmetry reduction) certify the strictly larger edge reserves
\[
d_{\rm cent,*} = \frac{317}{1000} > d_{\rm cent} = \frac{19}{60},
\qquad
d_{\rm caus,*} = \frac{119}{400} > d_{\rm caus} = \frac{19}{64}.
\]
The bare-form stability gaps at $k_B T = C_{\rm cent}\sqrt{\rho}$ and $k_B T = C_{\rm caus}\sqrt{\rho}$ are $1 - (1 - d_{\rm cent,*})/(1 - d_{\rm cent}) = 1/2050$ and $1 - (1 - d_{\rm caus,*})/(1 - d_{\rm caus}) = 1/1125$, respectively (corresponding to the unrounded reserve thresholds $C_{\rm cent,*} = \sqrt{683/500}/(2\pi) \in [186013935536718/10^{15},\ 186013935536719/10^{15}]$ and $C_{\rm caus,*} = \sqrt{281/200}/(2\pi) \in [188650644722914/10^{15},\ 188650644722915/10^{15}]$, also certified in `replay.py`).

**OPEN:** whether the sharp coupled finite-energy supremum $S$ over the centered finite-channel class equals $C_\infty$ or lies strictly inside $[C_\infty, C_{\rm cent}] \approx [0.18272625, 0.18605932]$ remains **OPEN**. Our initial **NUMERICAL** scan (`search_finite_energy.py`, `SEARCH.json`) across 50 single- and two-channel Einstein and Lorentzian finite-energy rows with flat, box, shell, and notched centered DOS profiles found no pairing eigenvalue $\ge 1$ at $k_B T = C_\infty\sqrt{\rho}$ (maximum gain $0.999491612731$ in the soft flat benchmark; maximum non-flat finite-energy gain $0.99945037$ in wide boxes approaching the flat limit). Neither a finite numerical search nor a form upper bound proves $S = C_\infty$.
## 1. Inherited model, normalization, and known single-channel / flat-DOS results

**PROVED under these inherited hypotheses (`../../OPTIMUM/RIGOR.md`, `../../OPTIMUM2/RIGOR.md`, `../../note-centered/PHYS-Q2/RIGOR.md`):**

Let $(X_{\rm ch}, \eta)$ be a finite channel space with positive DOS-height weights $\eta_x > 0$, fixed measurable DOS profiles $0 \le N_x(\xi) \le N_{\max,x}$, energy-independent radial vertices, and reciprocal nonnegative measures $d\mathfrak{m}_{xy}(E) = d\mathfrak{m}_{yx}(E) \ge 0$ supported on $E \in (0, \infty)$ with total moment matrix
\[
\mathcal{M}_{xy} = \int_{(0, \infty)} d\mathfrak{m}_{xy}(E),
\]
having finite positive operator norm $\rho = \|\mathcal{M}\|_{\ell^2(\eta) \to \ell^2(\eta)} \in (0, \infty)$. At each tested temperature $\tau = k_B T > 0$, the normal frequency rows (including zero transfer) have a finite bound $R(\tau) < \infty$ (no temperature-independent row bound is assumed), and each spectral energy $E > 0$ has dimensionless squared energy $t = (E / (2\pi\tau))^2 > 0$ and Matsubara transfer kernel $k_j(t) = (t + j^2)^{-1}$. An actual positive reflected normal branch solves the matched imaginary normal equation; a larger normal kernel is permitted when its excess over the pairing transfer is odd-folded nonnegative, a bounded pairing kernel dominated entrywise by the nonnegative attractive envelope cannot increase the onset, and a positive semidefinite static pairing-only subtraction is harmless on the matched branch (`PUBLISH/optimum/source/main.tex` §2). When a symmetry application is invoked, the normal propagator and attractive envelope preserve the chosen sector.

For each channel $x$ and Matsubara index $n \ge 0$, let $w_n = 2n + 1$ and let $b_{n,x} \in [0, 1]$ be the height-normalized DOS Poisson average at the positive reflected normal width $s_{n,x} \ge w_n$. Write $\ell_{n,x}=\Re Z_x(i\pi\tau w_n)$ for the radial energy shift; $\ell_{n,x}=0$ on the centered branch:
\[
b_{n,x} = \int_{-\infty}^{\infty} \frac{N_x(\xi)}{N_{\max,x}} \frac{\pi\tau s_{n,x}}{(\pi\tau s_{n,x})^2 + (\xi-\ell_{n,x})^2} \frac{d\xi}{\pi} \in [0, 1].
\]
On the causal shifted radial-normal class, write $F_x=g_x\circ Z_x$, where $g_x$ is the normalized DOS Herglotz transform, with a real subtraction if needed, and $Z_x$ is the causal upper-half-plane self-map representing the radial propagator. Its imaginary part $B_x(z)=\Im F_x(z)$ is harmonic and satisfies $0\le B_x(z)\le1$ throughout the upper half-plane: the normalized DOS is bounded by one, and its Poisson kernel integrates to one. Consequently both $B_x$ and $1-B_x$ are nonnegative harmonic functions; this does not assert that the holomorphic function $1-F_x$ is Herglotz. Applying Harnack's inequality to these two harmonic functions between $z_0=i\pi\tau$ and $z_m=i\pi\tau w_m$, with zero-function cases handled by continuity, proves the four-inequality causal Harnack cone $\mathcal{C}_{\rm caus}$ for $b_{n,x}=B_x(i\pi\tau w_n)$:
\begin{equation}
0 \le b_{m,x} \le 1,
\qquad
\frac{b_{0,x}}{w_m} \le b_{m,x} \le w_m b_{0,x},
\qquad
w_m b_{0,x} - (w_m - 1) \le b_{m,x} \le \frac{b_{0,x} + w_m - 1}{w_m}
\qquad (m \ge 1).
\label{eq:harnack}
\end{equation}
On the centered even-DOS class ($N_x(\xi) = N_x(-\xi) \in [0, N_{\max,x}]$), the self-consistent widths satisfy the two-sided width cone $w_m^{-1} \le s_{m,x}/s_{0,x} \le w_m$ (`../../CENTERED-GAP/RIGOR.md`, `../../SHARP-Q2/RIGOR.md`). The ratio of the two centered Poisson kernels is then in $[w_m^{-1},w_m]$; integrating both $N_x/N_{\max,x}$ and its complement gives all four inequalities in \eqref{eq:harnack}. Thus the centered Poisson class $\mathcal{C}_{\rm cent}$ is contained in $\mathcal{C}_{\rm caus}$. No actual-branch restriction is added relative to `OPTIMUM`; the synthetic causal relaxation explicitly requires all four Harnack bounds.

For any finite-support real test profile $v_{n,x}$ and normalized amplitudes $g_{n,x} = \sqrt{w_n b_{n,x}} v_{n,x}$, the folded stability quadratic form is
\begin{equation}
Q(v) = \sum_x \eta_x \|g_x\|_{\ell^2}^2 + \frac{1}{2(2\pi\tau)^2} \sum_{x,y} \eta_x \eta_y \int_{(0, \infty)} I_{xy}(t)\, d\mathfrak{m}_{xy}(E),
\label{eq:folded_form}
\end{equation}
where for one reciprocal channel pair $(x, y)$, one spectral energy $t > 0$, channel weights $b_n = b_{n,x}$, $c_m = b_{m,y}$, and real trial coordinates $X_n = v_{n,x}$, $Y_m = v_{m,y}$,
\[
I_{xy}(t) = \sum_{n,m \ge 0} b_n c_m \left[ p_{nm}(t)\bigl(X_n^2 + Y_m^2\bigr) - 2 S_{nm}(t) X_n Y_m \right],
\]
with $p_{nm}(t) = k_{|n-m|}(t) - k_{n+m+1}(t) > 0$ and $S_{nm}(t) = k_{|n-m|}(t) + k_{n+m+1}(t) > 0$. For any edge allocation ratio $r > 0$, define the bare edge norm
\[
B_r(X, Y) = \sum_{n \ge 0} w_n \bigl(r^{-1} b_n X_n^2 + r c_n Y_n^2\bigr) = \sum_{n \ge 0} \bigl(a_n^2 + d_n^2\bigr),
\]
in terms of the unit-bare-norm coordinates
\[
a_n = r^{-1/2} \sqrt{w_n b_n}\, X_n,
\qquad
d_m = r^{1/2} \sqrt{w_m c_m}\, Y_m.
\]
If an edge reserve $d_* \in (0, 1)$ satisfies
\begin{equation}
B_r(X, Y) + \frac{1}{2} I_{xy}(t) \ge d_* B_r(X, Y)
\qquad \forall\, r > 0,\ t > 0,
\label{eq:edge_reserve}
\end{equation}
for all admissible weight sequences $(b_n), (c_m)$, then choosing $h = (I - \mathcal{M}/R')^{-1}\mathbf{1} > 0$ with $R' \downarrow \rho$ and $r = h_x / h_y$ yields
\[
Q(v) \ge \left(1 - \frac{2(1 - d_*)\rho}{(2\pi\tau)^2}\right) \sum_x \eta_x \|g_x\|_{\ell^2}^2,
\]
so every temperature satisfying $(2\pi\tau)^2 > 2(1 - d_*)\rho$ is strictly stable (`../../OPTIMUM/RIGOR.md` §1, `../../OPTIMUM2/RIGOR.md` §2).

**What is already settled (`PROVED` in prior repository artifacts; not re-derived here):**
- In the single-channel scalar case ($d = 1$), `PUBLISH/finite-repo` (Zenodo 24364879) proves $k_B T_{\rm lin} < C_\infty \sqrt{\rho}$ for all finite-energy glue measures and all bounded energy-dependent DOS profiles. Why does that proof not immediately close $d \ge 2$ with unequal DOS profiles? Because for $d = 1$, $X_n = Y_n$ and $b_n = c_n$, so the ordinary-transfer term $\frac{1}{2} b_n b_m p_{nm}(X_n - X_m)^2$ has identical weights $b_n b_m$ on the diagonal $X_n^2 + X_m^2$ and cross term $2 X_n X_m$, matching the flat-DOS comparison operator up to a diagonal multiplier $\le 1$. When $d \ge 2$ and channel $x$ has a narrow DOS ($b_{n,x} \ll 1$) while channel $y$ has a wide DOS ($c_{m,y} \approx 1$), the edge weight $b_{n,x} c_{m,y}$ suppresses the bare normal mass $w_n b_{n,x} X_n^2$ of channel $x$ at the same linear rate $b_{n,x}$ as the cross-channel pairing term $b_{n,x} c_{m,y} S_{nm} X_n Y_m$.
- In `../../OPTIMUM2/RIGOR.md`, a single-stage head-to-tail completion ($n = 0$ against $m = 1, \dots, 20$) with independent Bernoulli head masks achieved $d_* = 247/1000$ ($C_{21} = \sqrt{377/250}/(2\pi)$). Why did it stop at $247/1000$? Two distinct mechanisms limited the single-stage bound:
  1. **Bernoulli mask decoupling slack:** reducing $(b_0, c_0)$ to $\{0, 1\}^2$ before adding the head-to-tail completions evaluated the square-completion denominator at the worst-case own-head weight $b_0 = 1$ ($\alpha + p_m/(2 r w_m)$ instead of $\alpha/b_0 + p_m/(2 r w_m)$), discarding the strong stabilization of small $b_0$ while simultaneously allowing $c_0$ to vary independently in the head reserve $1 - c_0 + e_2 c_0$.
  2. **Unassisted $m = 1$ tail bottleneck:** at small $t$ and imbalance $r \approx 1/3$, the $m = 1$ completion $D_1(r, t)$ has threshold $r \alpha = 2/w_1^2 = 2/9$ (equivalently $r \alpha / 2 = 1/w_1^2 = 1/9$). With static tail allocation $\alpha \le 1 - d_* - q$, as soon as $d_*$ approaches $0.30$, $\alpha \le 0.634$, so at $r = 1/3$ one has $r \alpha \le 0.2113 < 2/9 = 0.2222\ldots$, making $D_1(1/3, t) < 0$. In the actual stability form, however, the $m = 1$ mode is itself stabilized by the positive tail-tail ordinary-transfer squares $\frac{1}{2} b_1 c_k p_{1,k}(X_1 - Y_k)^2$ for $k \ge 2$.

---

## 2. Two-stage frequency cascade ($m = 1 \leftrightarrow k \ge 2$, then $n = 0 \leftrightarrow m \ge 1$)

**PROVED.** Let $M = 20$ and $w_n = 2n + 1$. In unit-bare-norm coordinates $a_n = r^{-1/2}\sqrt{w_n b_n}\,X_n$ and $d_m = r^{1/2}\sqrt{w_m c_m}\,Y_m$ (when $b_n = 0$ or $c_m = 0$, the corresponding pairing terms vanish and the bare reserve is trivially $1 \ge d_*$), the quadratic form $B_r(X, Y) + \frac{1}{2}I_{xy}(t)$ decomposes into four disjoint sectors:

1. **Head-head sector ($n = m = 0$):**
   \[
   Q_{00} = \left(1 + \frac{r c_0 p_{00}}{2}\right) a_0^2 + \left(1 + \frac{r^{-1} b_0 p_{00}}{2}\right) d_0^2 - \sqrt{b_0 c_0}\, S_{00} a_0 d_0,
   \]
   where $p_{00} = 2\beta = [t(1 + t)]^{-1}$, $S_{00} = 2(\beta + u)$, $\beta = [2t(1 + t)]^{-1}$, and $u = (1 + t)^{-1}$.

2. **Stage-0 head-tail sectors ($n = 0, m \ge 1$ and $m = 0, n \ge 1$):**
   Writing $p_{0,m} = k_m - k_{m+1} = \frac{w_m}{(t + m^2)(t + (m+1)^2)}$ and $s_{0,m} = k_m + k_{m+1}$,
   \[
   Q_{0,{\rm tail}} = \sum_{m \ge 1} \left[ \frac{r c_m p_{0,m}}{2} a_0^2 + \frac{r^{-1} b_0 p_{0,m}}{2 w_m} d_m^2 - s_{0,m}\sqrt{\frac{b_0 c_m}{w_m}}\, a_0 d_m \right]
   + \sum_{n \ge 1} \left[ \frac{r^{-1} b_n p_{0,n}}{2} d_0^2 + \frac{r c_0 p_{0,n}}{2 w_n} a_n^2 - s_{0,n}\sqrt{\frac{c_0 b_n}{w_n}}\, d_0 a_n \right].
   \]

3. **Stage-1 first-tail-to-higher-tail sectors ($n = 1, m = k \in \{2, \dots, M\}$ and $m = 1, n = k \in \{2, \dots, M\}$):**
   Instead of discarding the ordinary-transfer difference square at $(1, k)$ and bounding its reflected transfer $k_{k+2}$ by Schur's test, we retain the full $(1, k)$ interaction
   \[
   p_{1,k} = k_{k-1} - k_{k+2} = \frac{3 w_k}{\bigl(t + (k-1)^2\bigr)\bigl(t + (k+2)^2\bigr)},
   \qquad
   s_{1,k} = k_{k-1} + k_{k+2},
   \]
   giving
   \[
   Q_{1,{\rm high}} = \sum_{k=2}^{M} \left[ \frac{r c_k p_{1,k}}{6} a_1^2 + \frac{r^{-1} b_1 p_{1,k}}{2 w_k} d_k^2 - s_{1,k}\sqrt{\frac{b_1 c_k}{3 w_k}}\, a_1 d_k \right]
   + \sum_{k=2}^{M} \left[ \frac{r^{-1} b_k p_{1,k}}{6} d_1^2 + \frac{r c_1 p_{1,k}}{2 w_k} a_k^2 - s_{1,k}\sqrt{\frac{c_1 b_k}{3 w_k}}\, d_1 a_k \right].
   \]

4. **Remaining tail-tail sector (all $(n, m) \in \mathbb{N}_{\ge 1}^2$ except $\{1\} \times \{2, \dots, M\}$ and $\{2, \dots, M\} \times \{1\}$):**
   For every remaining pair $(n, m)$, the exact algebraic identity (`reflected_tail_decomposition` in `Finite.lean`)
   \[
   (a - b)(X^2 + Y^2) - 2(a + b) X Y = (a-b)(X - Y)^2 - 4 b X Y
   \]
   with $a = k_{|n-m|}(t)$ and $b = k_{n+m+1}(t)$ splits the $(n, m)$ contribution to $\frac{1}{2}I_{xy}(t)$ into $\frac{1}{2}b_n c_m p_{nm}(X_n - Y_m)^2 - 2 b_n c_m k_{n+m+1} X_n Y_m$, where $p_{nm}=k_{|n-m|}-k_{n+m+1}\ge0$. Dropping precisely this nonnegative difference square and rewriting the reflected cross term in unit-bare-norm coordinates $a_n = r^{-1/2}\sqrt{w_n b_n}\,X_n$ and $d_m = r^{1/2}\sqrt{w_m c_m}\,Y_m$ via weighted AM--GM gives
   \[
   -2 b_n c_m k_{n+m+1} X_n Y_m
   = -\frac{2\sqrt{b_n c_m}}{\sqrt{w_n w_m}}\, k_{n+m+1} a_n d_m
   \ge -\frac{c_m}{w_m} k_{n+m+1} a_n^2 - \frac{b_n}{w_n} k_{n+m+1} d_m^2
   \ge -\frac{a_n^2}{w_m(n+m+1)^2} - \frac{d_m^2}{w_n(n+m+1)^2},
   \]
   since $b_n, c_m \in [0, 1]$ and $k_{n+m+1}(t) \le (n+m+1)^{-2}$. Summing the coefficients of $a_n^2$ and $d_n^2$ over the remaining pairs $(n, m)$ yields the exact rational reduced Hankel Schur row sums:
   - For $1 \le n \le M$:
     \[
     q_{{\rm red},n} := \sum_{\substack{m=1 \\ (n,m) \notin {\rm Stage\,1}}}^{M} \frac{1}{(2m+1)(n+m+1)^2} + \frac{1}{43(n+21)},
     \]
     since for $m \ge 21$, $2m + 1 \ge 43$ and $\sum_{m=21}^{\infty} (n+m+1)^{-2} \le \int_{20}^{\infty} (n+x+1)^{-2}\,dx = (n+21)^{-1}$.
   - For $n > M$: no Stage-1 entries are removed from row $n$, so the full Hankel Schur bound
     \[
     q_{{\rm full},n} \le q_{{\rm full},1} \le \sum_{m=1}^{16} \frac{1}{(2m+1)(m+2)^2} + \frac{1}{4\cdot 16^2} < \frac{33}{500} =: q_{\rm cap}
     \]
     applies.

After reserving $d_* \in (0, 1)$ of the bare norm of every mode $n \ge 0$ and subtracting the reduced Hankel Schur cost from each tail mode $n \ge 1$, the remaining available bare-norm share on mode $n \in \{1, \dots, M\}$ is at least
\[
{\rm rem}_n = 1 - d_* - q_{{\rm red},n} > 0.
\]
We split ${\rm rem}_n$ between Stage 0 and Stage 1 with the exact rational fraction $\sigma = 7/10$:
\[
\alpha_{0,1} = {\rm rem}_1, \quad \beta_{1,1} = 0;
\qquad
\alpha_{0,k} = \frac{7}{10}\,{\rm rem}_k, \quad \beta_{1,k} = \frac{3}{10}\,{\rm rem}_k \quad (2 \le k \le M).
\]
For $n > M$, the entire remaining bare norm $1 - d_* - q_{\rm cap} =: \alpha_{\rm ext}$ goes to Stage 0 ($\beta_{1,n} = 0$).

### Lemma 1 (Stage-1 completion and cancellation of $b_1, c_1$)
**PROVED (`stage1_transfer_identity` and `stage1_completion_numerator` in `Finite.lean`).**
For each $k \in \{2, \dots, M\}$, combine the Stage-1 bare-norm share $\beta_{1,k} d_k^2$ with the $(a_1, d_k)$ block in $Q_{1,{\rm high}}$:
\[
\frac{r c_k p_{1,k}}{6} a_1^2 + \left(\beta_{1,k} + \frac{r^{-1} b_1 p_{1,k}}{2 w_k}\right) d_k^2 - s_{1,k}\sqrt{\frac{b_1 c_k}{3 w_k}}\, a_1 d_k.
\]
Completing the square in $d_k$ and using $b_1 \in (0, 1] \implies \frac{\beta_{1,k}}{b_1} \ge \beta_{1,k}$ yields a lower bound $c_k D_{1,k}^{(X)}(r, t)\, a_1^2$, where
\begin{equation}
D_{1,k}^{(X)}(r, t)
:= \frac{r p_{1,k}}{6} - \frac{s_{1,k}^2}{12 w_k \bigl(\beta_{1,k} + r^{-1} p_{1,k}/(2 w_k)\bigr)}
= \frac{2 p_{1,k}\bigl(r w_k \beta_{1,k} - \frac{2}{3 w_k}\bigr)}{12 w_k \beta_{1,k} + 6 r^{-1} p_{1,k}}.
\label{eq:D1}
\end{equation}
*Proof.* The difference of squares satisfies
\[
s_{1,k}^2 - p_{1,k}^2 = 4 k_{k-1} k_{k+2} = \frac{4}{\bigl(t + (k-1)^2\bigr)\bigl(t + (k+2)^2\bigr)} = \frac{4 p_{1,k}}{3 w_k},
\]
so $3 w_k (s_{1,k}^2 - p_{1,k}^2) = 4 p_{1,k}$. Putting the two terms of \eqref{eq:D1} over the common denominator $12 w_k \beta_{1,k} + 6 r^{-1} p_{1,k}$ gives the numerator
\[
2 r w_k \beta_{1,k} p_{1,k} + p_{1,k}^2 - s_{1,k}^2 = 2 p_{1,k}\left(r w_k \beta_{1,k} - \frac{2}{3 w_k}\right).
\]
Notice that $b_1 \in (0, 1]$ cancels from the leading term $\frac{r c_k p_{1,k}}{6} a_1^2$ because $a_1^2 = r^{-1} w_1 b_1 X_1^2$ already carries the factor $b_1$, while in the negative square-completion term $\frac{s_{1,k}^2 b_1}{12 w_k (\beta_{1,k} + r^{-1} b_1 p_{1,k}/(2 w_k))}$ the increasing function $b_1 \mapsto \frac{b_1}{\beta_{1,k} + C b_1}$ on $[0, 1]$ attains its maximum at $b_1 = 1$. By symmetry ($r \leftrightarrow r^{-1}$, $b \leftrightarrow c$), completing the square in $a_k$ adds $\sum_{k=2}^{M} b_k D_{1,k}^{(Y)}(r, t)\, d_1^2 = \sum_{k=2}^{M} b_k D_{1,k}^{(X)}(r^{-1}, t)\, d_1^2$ to the Stage-0 bare norm of $d_1^2$. $\square$

For $r \ge 1/3$ and $k \ge 2$, $w_k \ge 5$, so the threshold in \eqref{eq:D1} is $\beta_{1,k} > \frac{2}{3 r w_k^2} \le \frac{2}{25} = 0.08$. In our allocations $\beta_{1,k} = \frac{3}{10}{\rm rem}_k > 0.19 > 0.08$, so every $D_{1,k}^{(X)}(r, t)$ and $D_{1,k}^{(Y)}(r, t)$ ($2 \le k \le 20$) is strictly positive on $[1/3, 3] \times [0, 1/2]$, strictly increasing in $r$ (resp. $r^{-1}$) and in $p_{1,k}$, and therefore strictly decreasing in $t$.

Consequently, before Stage 0 begins, the effective Stage-0 normalized bare-mass coefficients of $a_1^2$ and $d_1^2$ are boosted to
\begin{equation}
\alpha_1^{{\rm eff},(X)}(r, t; c) = \alpha_{0,1} + \sum_{k=2}^{M} c_k D_{1,k}^{(X)}(r, t),
\qquad
\alpha_1^{{\rm eff},(Y)}(r, t; b) = \alpha_{0,1} + \sum_{k=2}^{M} b_k D_{1,k}^{(Y)}(r, t).
\label{eq:alpha1_eff}
\end{equation}

### Lemma 2 (Stage-0 completion with retained own-head weight $b_0, c_0$)
**PROVED (`stage0_transfer_identity` and `stage0_completion_numerator` in `Finite.lean`).**
For each $m \ge 1$, let $\alpha_m^{(Y)} = \alpha_1^{{\rm eff},(Y)}(r, t; b)$ for $m = 1$, $\alpha_m^{(Y)} = \alpha_{0,m}$ for $2 \le m \le M$, and $\alpha_m^{(Y)} = \alpha_{\rm ext}$ for $m > M$. Combining $\alpha_m^{(Y)} d_m^2$ with the $(a_0, d_m)$ block in $Q_{0,{\rm tail}}$ and completing the square in $d_m$ yields a contribution $c_m \hat{D}_{0,m}^{(X)}(r, t; b_0)\, a_0^2$, where for any positive lower bound $\tilde{\alpha}_m \le \alpha_m^{(Y)} / b_0$:
1. **Pointwise lower bound:**
   \begin{equation}
   \hat{D}_{0,m}^{(X)}(r, t; b_0)
   \ge \frac{r p_{0,m}}{2} - \frac{s_{0,m}^2}{4 w_m \bigl(\tilde{\alpha}_m + r^{-1} p_{0,m}/(2 w_m)\bigr)}
   = \frac{2 p_{0,m}\bigl(r w_m \tilde{\alpha}_m - \frac{2}{w_m}\bigr)}{4 w_m \tilde{\alpha}_m + 2 r^{-1} p_{0,m}}.
   \label{eq:D0}
   \end{equation}
2. **Sign-split rectangle bounds on $(r, t) \in [r_{\rm lo}, r_{\rm hi}] \times [t_{\rm lo}, t_{\rm hi}]$:**
   - If $r_{\rm lo} w_m \tilde{\alpha}_m - 2/w_m \ge 0$, then \eqref{eq:D0} is nonnegative and increasing in $r$ and $p_{0,m}$, so
     \begin{equation}
     \hat{D}_{0,m}^{(X)}(r, t; b_0) \ge \frac{2 p_{0,m}(t_{\rm hi})\bigl(r_{\rm lo} w_m \tilde{\alpha}_m - \frac{2}{w_m}\bigr)}{4 w_m \tilde{\alpha}_m + 2 r_{\rm lo}^{-1} p_{0,m}(t_{\rm hi})} \ge 0,
     \label{eq:D0_pos}
     \end{equation}
     which multiplies any nonnegative lower envelope for $c_m / c_0$.
   - If $r_{\rm lo} w_m \tilde{\alpha}_m - 2/w_m < 0$, dropping $2 r^{-1} p_{0,m} > 0$ from the denominator of \eqref{eq:D0} and multiplying by the upper Harnack ratio $c_m / c_0 \le w_m$ from \eqref{eq:harnack} gives the outward negative lower bound
     \begin{equation}
     \frac{c_m}{c_0} \hat{D}_{0,m}^{(X)}(r, t; b_0)
     \ge w_m \frac{2 p_{0,m}(t_{\rm lo})\bigl(r_{\rm lo} w_m \tilde{\alpha}_m - 2/w_m\bigr)}{4 w_m \tilde{\alpha}_m} < 0.
     \label{eq:D0_neg}
     \end{equation}
     (Notice why the upper Harnack law $c_m / c_0 \le w_m$ in \eqref{eq:harnack} is essential here: for a negative completion $\hat{D}_{0,m}^{(X)} < 0$, bounding $\frac{c_m}{c_0}\hat{D}_{0,m}^{(X)}$ from below requires an upper bound on $c_m / c_0$, not a lower bound; without $c_m \le w_m c_0$, a vector in the synthetic lower-only cone with $c_0 = 1/10, c_1 = 1, r = 1/3, t = 1/10$ and $\tilde{\alpha}_1 = 3/5$ gives $(c_1/c_0)\hat{D}_{0,1}^{(X)}(1/3, 1/10) = -500/6309$, strictly below the expression $-50/1353$ evaluated at $c_1/c_0 = w_1 = 3$.)
*Proof.* Clearing denominators in $s_{0,m}^2 - p_{0,m}^2 = 4 k_m k_{m+1}$ via the polynomial identity `stage0_transfer_identity` in `Finite.lean` gives $w_m(s_{0,m}^2 - p_{0,m}^2) = 4 p_{0,m}$, so the exact square completion in $d_m$ leaves the coefficient of $a_0^2$ equal to
\[
\frac{r c_m p_{0,m}}{2} - \frac{s_{0,m}^2 b_0 c_m}{4 w_m \bigl(\alpha_m^{(Y)} + r^{-1} b_0 p_{0,m}/(2 w_m)\bigr)}
= c_m \left[ \frac{r p_{0,m}}{2} - \frac{s_{0,m}^2}{4 w_m \bigl(\alpha_m^{(Y)}/b_0 + r^{-1} p_{0,m}/(2 w_m)\bigr)} \right].
\]
Since the bracketed expression is strictly increasing with respect to $\alpha_m^{(Y)}/b_0 > 0$, replacing $\alpha_m^{(Y)}/b_0$ by $\tilde{\alpha}_m \le \alpha_m^{(Y)}/b_0$ and applying `stage0_completion_numerator` gives \eqref{eq:D0}. Both rectangle branches \eqref{eq:D0_pos} and \eqref{eq:D0_neg} follow immediately from the sign of the numerator and the monotonicity of $p_{0,m}(t)$ on $[t_{\rm lo}, t_{\rm hi}]$.
For $m > M = 20$ on $r \in [1/3, 3]$, $w_m \ge 43$ and $\tilde{\alpha}_m \ge \alpha_{\rm ext} > 3/5$, so $r w_m \tilde{\alpha}_m \ge \frac{1}{3}\cdot 43 \cdot \frac{3}{5} = \frac{43}{5} > \frac{2}{43} = \frac{2}{w_m}$, and all unretained $m > 20$ Stage-0 completions are strictly positive. $\square$

---

## 3. Joint $(b_0, c_0)$ sub-box decomposition and decoupled $2 \times 2$ head criterion

**PROVED (`head_square` and `decoupled_head_criterion` in `Finite.lean`).**
After extracting the Stage-1 and Stage-0 completions and reserving $d_* (a_0^2 + d_0^2)$ from the heads, the remaining head quadratic form in $(a_0, d_0)$ is
\[
Q_{\rm head}(a_0, d_0)
= \left(1 - d_* + \frac{r c_0 p_{00}}{2} + \sum_{m=1}^{M} c_m \hat{D}_{0,m}^{(X)}(r, t; b_0)\right) a_0^2
+ \left(1 - d_* + \frac{r^{-1} b_0 p_{00}}{2} + \sum_{m=1}^{M} b_m \hat{D}_{0,m}^{(Y)}(r, t; c_0)\right) d_0^2
- \sqrt{b_0 c_0}\, S_{00} a_0 d_0.
\]
Define the normalized per-unit-opposite-head reserves
\begin{equation}
\Phi_X(r, t; b, c) := \frac{1 - d_*}{c_0} + \sum_{m=1}^{M} \frac{c_m}{c_0} \hat{D}_{0,m}^{(X)}(r, t; b_0),
\qquad
\Phi_Y(r, t; b, c) := \frac{1 - d_*}{b_0} + \sum_{m=1}^{M} \frac{b_m}{b_0} \hat{D}_{0,m}^{(Y)}(r, t; c_0).
\label{eq:Phi_XY}
\end{equation}
In the scaled head coordinates $\tilde{a}_0 = \sqrt{c_0}\,a_0$ and $\tilde{d}_0 = \sqrt{b_0}\,d_0$ (if $b_0 = 0$ or $c_0 = 0$, the cross term vanishes and $\Phi_X, \Phi_Y > 0$ follows immediately from $1 - d_* > 0$), the $2 \times 2$ matrix of $Q_{\rm head}$ is
\[
M_{\rm head} = \begin{pmatrix}
\Phi_X + \beta r & -(\beta + u) \\
-(\beta + u) & \Phi_Y + \beta r^{-1}
\end{pmatrix},
\qquad u = \frac{1}{1 + t}, \quad \beta = \frac{1}{2t(1 + t)}.
\]

### Lemma 3 (Decoupled $2 \times 2$ head positivity and linearized $r$-convexity)
**PROVED (`decoupled_head_criterion` in `Finite.lean`).** Suppose $\Phi_X > 0$ and $\Phi_Y > 0$. Setting $A = 2t(1+t) > 0$, clearing the positive factors $A^2 > 0$ and $r > 0$ via the polynomial identity `decoupled_head_criterion` in `Finite.lean`,
\[
(A\Phi_X + r)(A r \Phi_Y + 1) - r(1 + 2t)^2
= A^2 r \Phi_X \Phi_Y + A(\Phi_X + r^2\Phi_Y) - 4 r t(1+t)
= 2 r t(1+t)\bigl[2t(1+t)\Phi_X\Phi_Y + (r^{-1}\Phi_X + r\Phi_Y) - 2\bigr],
\]
shows that $M_{\rm head} \succeq 0$ if and only if $\det M_{\rm head} \ge 0$, because its diagonals are positive. Dividing the displayed polynomial by $A r > 0$ gives the equivalent condition
\begin{equation}
2t(1 + t)\Phi_X \Phi_Y + (r^{-1}\Phi_X + r\Phi_Y) - 2 \ge 0.
\label{eq:det_head}
\end{equation}
Moreover, on any rectangle $(r, t) \in [r_{\rm lo}, r_{\rm hi}] \times [t_{\rm lo}, t_{\rm hi}] \subset [1/3, 3] \times [0, 1/2]$, if $\Phi_X^{\rm lo} > 0$ and $\Phi_Y^{\rm lo} > 0$ are uniform lower bounds for $\Phi_X$ and $\Phi_Y$ on the rectangle, then:
1. For the midpoint $r_0 = \frac{1}{2}(r_{\rm lo} + r_{\rm hi})$, the first-order convexity inequality $r^{-1} \ge 2 r_0^{-1} - r_0^{-2} r$ gives
   \[
   r^{-1} \Phi_X^{\rm lo} + r \Phi_Y^{\rm lo}
   \ge \frac{2 \Phi_X^{\rm lo}}{r_0} + r\left(\Phi_Y^{\rm lo} - \frac{\Phi_X^{\rm lo}}{r_0^2}\right)
   \ge g_{\rm lo} := \frac{2\Phi_X^{\rm lo}}{r_0} + r_*\left(\Phi_Y^{\rm lo}-\frac{\Phi_X^{\rm lo}}{r_0^2}\right),
   \]
   where $r_*=r_{\rm lo}$ when $\Phi_Y^{\rm lo}-\Phi_X^{\rm lo}/r_0^2\ge0$, and $r_*=r_{\rm hi}$ otherwise. This is the minimum of the affine supporting tangent over the entire interval, exactly as evaluated in both checkers.
2. Since $\Phi_X^{\rm lo} > 0$, $\Phi_Y^{\rm lo} > 0$, and $t \mapsto 2t(1+t)$ is nondecreasing on $t \ge t_{\rm lo} \ge 0$, we have $2t(1+t)\Phi_X\Phi_Y \ge 2 t_{\rm lo}(1 + t_{\rm lo})\Phi_X^{\rm lo}\Phi_Y^{\rm lo} \ge 0$. Therefore \eqref{eq:det_head} holds across the entire rectangle whenever the single unconditional sufficient condition
   \[
   2 t_{\rm lo}(1 + t_{\rm lo})\Phi_X^{\rm lo}\Phi_Y^{\rm lo} + g_{\rm lo} \ge 2
   \]
   holds (`replay.py` and `independent.py` verify this exact monotone sufficient condition on every sub-box of every accepted leaf). $\square$
### Partition of $(b_0, c_0) \in [0, 1]^2$ into $9 \times 9$ sub-boxes
We partition $[0, 1]$ into $S = 9$ contiguous rational sub-intervals $I_a = [\beta_{a,\rm lo}, \beta_{a,\rm hi}]$:
\[
\left[0, \frac{2}{3}\right],\ \left[\frac{2}{3}, \frac{4}{5}\right],\ \left[\frac{4}{5}, \frac{17}{20}\right],\ \left[\frac{17}{20}, \frac{22}{25}\right],\ \left[\frac{22}{25}, \frac{9}{10}\right],\ \left[\frac{9}{10}, \frac{23}{25}\right],\ \left[\frac{23}{25}, \frac{47}{50}\right],\ \left[\frac{47}{50}, \frac{24}{25}\right],\ \left[\frac{24}{25}, 1\right].
\]
On each sub-interval $I_a$:
- Since $b_0 \le \beta_{a,\rm hi}$, we have $\alpha_{0,m}/b_0 \ge \alpha_{0,m}/\beta_{a,\rm hi}$ for every $2 \le m \le M$.
- Whenever all coefficients $D_{1,k}$ and $\hat{D}_{0,m}$ are nonnegative, each reserve sum of the form
  \[
  \frac{A}{b_0} + \sum_{m=1}^{M} \frac{b_m}{b_0} C_m
  \]
  is bounded below by substituting a piecewiese-linear lower envelope $b_m \ge \max_\ell(\lambda_{\ell,m} b_0 + c_{\ell,m})$ (where $c_{\ell,m} \le 0$):
  \[
  \frac{A}{b_0} + \sum_{m=1}^{M} C_m \max_\ell\left(\lambda_{\ell,m} + \frac{c_{\ell,m}}{b_0}\right).
  \]
  In the reciprocal variable $u_0 = 1/b_0 \in [\beta_{a,\rm hi}^{-1}, \beta_{a,\rm lo}^{-1}]$ (with $\beta_{0,\rm lo}^{-1} = +\infty$, where the positive slope $A - \sum_m (w_m - 1) C_m$ or $A > 0$ dominates as $u_0 \to \infty$), this lower bound is a maximum of affine functions of $u_0$, hence is **convex in $u_0$** and takes its exact minimum on $I_a$ at the finitely many breakpoints of $\max_\ell(\lambda_{\ell,m} + c_{\ell,m} u_0)$ restricted to $I_a$!

---

## 4. Centered Poisson slope-cone envelope

**PROVED.** On the centered even-DOS class, the self-consistent widths obey the two-sided width cone $w_m^{-1} \le s_{m,x}/s_{0,x} \le w_m$ (`../../CENTERED-GAP/RIGOR.md`, `../../SHARP-Q2/RIGOR.md`). Rescaling energy by $\pi\tau s_{0,x}$, the normalized even DOS profile $\tilde{N}_x(\xi) \in [0, 1]$ has first Matsubara average and $m$-th Matsubara average
\[
b_{0,x} = \frac{2}{\pi} \int_0^\infty \frac{\tilde{N}_x(\xi)}{1 + \xi^2}\,d\xi,
\qquad
b_{m,x} = \frac{2}{\pi} \int_0^\infty \frac{\tilde{N}_x(\xi)\, \sigma_m}{\sigma_m^2 + \xi^2}\,d\xi,
\qquad \sigma_m := \frac{s_{m,x}}{s_{0,x}} \in [w_m^{-1}, w_m].
\]

### Proposition 1 (Centered Poisson lower envelope and tangent inequalities)
For every $\sigma_m \in [w_m^{-1}, w_m]$ and every measurable $\tilde{N}_x : [0, \infty) \to [0, 1]$ with first average $b_0 \in [0, 1]$,
\begin{equation}
b_{m,x} \ge \phi_{w_m}(b_0) := \frac{2}{\pi} \arctan\!\left(\frac{1}{w_m} \tan\frac{\pi b_0}{2}\right).
\label{eq:phi_wm}
\end{equation}
Furthermore, $\phi_w : [0, 1] \to [0, 1]$ is strictly increasing and strictly convex on $[0, 1]$ for every $w > 1$. For any reference point $\xi_0 > 0$, setting $b_0(\xi_0) = \frac{2}{\pi}\arctan(\xi_0)$ and
\begin{equation}
\lambda_m(\xi_0) := \phi_{w_m}'\bigl(b_0(\xi_0)\bigr) = \frac{w_m(1 + \xi_0^2)}{w_m^2 + \xi_0^2},
\qquad
c_m(\xi_0) := \frac{2}{\pi}\arctan\!\left(\frac{\xi_0}{w_m}\right) - \lambda_m(\xi_0)\,\frac{2}{\pi}\arctan(\xi_0),
\label{eq:tangent_line}
\end{equation}
every centered weight sequence satisfies the linear lower bound
\begin{equation}
b_{m,x} \ge \lambda_m(\xi_0)\, b_{0,x} + c_m(\xi_0),
\qquad c_m(\xi_0) < 0.
\label{eq:tangent_ineq}
\end{equation}
*Proof.* Let $W = \tan(\pi b_0 / 2) \in [0, \infty]$ and consider the kernel ratio $R_\sigma(\xi) = \frac{\sigma/(\sigma^2 + \xi^2)}{1/(1 + \xi^2)} = \frac{\sigma(1 + \xi^2)}{\sigma^2 + \xi^2}$ on $\xi \in [0, \infty)$ at fixed $b_0$ and fixed $\sigma = \sigma_m \in [w_m^{-1}, w_m]$:
- **Case 1 ($\sigma \in [1, w_m]$):** $R_\sigma(\xi)$ is nondecreasing in $\xi \in [0, \infty)$. By the bathtub / Neyman--Pearson rearrangement lemma (`../../OPTIMUM2/RIGOR.md` Proposition 1), among all measurable $\tilde{N}_x \in [0, 1]$ with $\frac{2}{\pi}\int_0^\infty \frac{\tilde{N}_x(\xi)}{1 + \xi^2}\,d\xi = b_0$, the integral $b_{m,x}$ is minimized by the central box indicator $\tilde{N}_x(\xi) = \mathbf{1}_{[0, W]}(\xi)$, giving $b_{m,x} \ge \frac{2}{\pi}\arctan(W/\sigma) \ge \frac{2}{\pi}\arctan(W/w_m) = \phi_{w_m}(b_0)$ since $\sigma \le w_m$.
- **Case 2 ($\sigma \in [w_m^{-1}, 1)$):** $R_\sigma(\xi)$ is strictly decreasing in $\xi \in [0, \infty)$. By the bathtub lemma, $b_{m,x}$ at fixed $b_0$ is minimized by the outer ray indicator $\tilde{N}_x(\xi) = \mathbf{1}_{[a, \infty)}(\xi)$ with $\frac{2}{\pi}(\frac{\pi}{2} - \arctan a) = b_0$, so $a = \cot(\pi b_0 / 2) = W^{-1}$. Evaluating the integral on $[a, \infty)$ gives
  \[
  b_{m,x} \ge \frac{2}{\pi}\int_a^\infty \frac{\sigma}{\sigma^2 + \xi^2}\,d\xi
  = \frac{2}{\pi}\left(\frac{\pi}{2} - \arctan\frac{a}{\sigma}\right)
  = \frac{2}{\pi}\arctan\frac{\sigma}{a}
  = \frac{2}{\pi}\arctan(\sigma W)
  \ge \frac{2}{\pi}\arctan(W / w_m)
  = \phi_{w_m}(b_0),
  \]
  since $\sigma \ge w_m^{-1}$ and $\arctan$ is strictly increasing. Thus \eqref{eq:phi_wm} holds across the entire two-sided width cone $\sigma_m \in [w_m^{-1}, w_m]$.

Differentiating $\phi_w(b_0)$ via $\xi_0 = \tan(\pi b_0 / 2)$ ($d\xi_0 / db_0 = \frac{\pi}{2}(1 + \xi_0^2)$) gives
\[
\phi_w'(b_0) = \frac{2}{\pi} \frac{1/w}{1 + (\xi_0/w)^2} \cdot \frac{\pi}{2}(1 + \xi_0^2) = \frac{w(1 + \xi_0^2)}{w^2 + \xi_0^2},
\qquad
\phi_w''(b_0) = \frac{\pi w(w^2 - 1)\xi_0(1 + \xi_0^2)}{(w^2 + \xi_0^2)^2} > 0
\quad (\xi_0 > 0,\ w > 1).
\]
By strict convexity of $\phi_{w_m}$ on $[0, 1]$, its supporting tangent line at $b_0(\xi_0)$ lies globally below $\phi_{w_m}(b_0)$ on $[0, 1]$, and its intercept at $b_0 = 0$ is $c_m(\xi_0) = \phi_{w_m}(b_0(\xi_0)) - \phi_{w_m}'(b_0(\xi_0))b_0(\xi_0) < \phi_{w_m}(0) = 0$. Replacing $c_m(\xi_0)$ by any rational lower bound $c_{m,\rm lo}(\xi_0) \le c_m(\xi_0)$ preserves both the lower bound \eqref{eq:tangent_ineq} and $c_{m,\rm lo}(\xi_0) \le 0$. $\square$

In both `replay.py` and `independent.py`, for the centered class we combine the two universal Harnack lines ($\lambda = 1/w_m, c = 0$ and $\lambda = w_m, c = -(w_m - 1)$) with three certified tangent lines at $\xi_0 \in \{5, 8, 12\}$, whose rational intercepts $c_{m,\rm lo}(\xi_0)$ are enclosed from below using certified rational alternating-series bounds for $\pi$ and $\arctan$.

---

## 5. Exterior regions ($t \ge 1/2$ and $r \notin [1/3, 3]$)

**PROVED (`exact_budgets_and_ceilings` in `Finite.lean`, Section 4 of `replay.py`).**
Outside the compact residual domain $(r, t) \in [1/3, 3] \times [0, 1/2]$, even the unassisted Stage-0 bound with binary masks $(b_0, c_0) \in \{0, 1\}^2$ and universal Harnack inequalities exceeds both $d_{\rm cent,*} = 317/1000$ and $d_{\rm caus,*} = 119/400$:

1. **Unallocated exterior Stage-0 tail share:**
   When no Stage-1 cascade is used, the full Hankel Schur bound $q_{\rm cap} = 33/500$ leaves on every tail mode $m \ge 1$ the Stage-0 bare-mass share
   \[
   \alpha_{\rm ext,cent} = 1 - \frac{317}{1000} - \frac{33}{500} = \frac{617}{1000} > \frac{3}{5},
   \qquad
   \alpha_{\rm ext,caus} = 1 - \frac{119}{400} - \frac{33}{500} = \frac{1273}{2000} > \frac{617}{1000}.
   \]

2. **Far imbalance $r \in (0, 1/4]$ (and $r \ge 4$ by $r \leftrightarrow r^{-1}$ symmetry):**
   On $r \le 1/4$, $R = r + r^{-1} \ge 17/4$, so the two-head reserve is at least $e_2(r, t) \ge 1 - 2/R \ge 9/17$. For every $m \ge 1$ and $t \ge 0$, the unassisted Stage-0 completion satisfies
   \[
   D_{0,m}(r, t) = \frac{2 p_{0,m}(t)\bigl(r w_m \alpha_{\rm ext} - 2/w_m\bigr)}{4 w_m \alpha_{\rm ext} + 2 r^{-1} p_{0,m}(t)}
   \ge \frac{2 p_{0,m}(t)(-2/w_m)}{4 w_m \alpha_{\rm ext}}
   = -\frac{p_{0,m}(t)}{w_m^2 \alpha_{\rm ext}}
   \ge -\frac{p_{0,m}(0)}{w_m^2 \alpha_{\rm ext}}
   = -\frac{1}{\alpha_{\rm ext}(2m+1)m^2(m+1)^2}.
   \]
   Since $(2m+1)m^2(m+1)^2 > 2m^5$, summing the first $20$ terms in exact rational arithmetic and bounding the $m \ge 21$ tail by $\frac{1}{2}\sum_{m=21}^\infty m^{-5} \le \frac{1}{2}\int_{20}^\infty x^{-5}\,dx = \frac{1}{8\cdot 20^4} = \frac{1}{1280000}$ gives (`replay.py` lines 214--221)
   \[
   \sum_{m=1}^{\infty} \frac{1}{(2m+1)m^2(m+1)^2}
   \le \sum_{m=1}^{20} \frac{1}{(2m+1)m^2(m+1)^2} + \frac{1}{8\cdot 20^4}
   = \frac{221683723759377431023}{2453474125083421440000}
   \approx 0.0903550 < \frac{7}{75} = 0.09\overline{3}.
   \]
   Therefore the total Stage-0 loss is strictly bounded by $L < \frac{7/75}{\alpha_{\rm ext}}$, leaving
   \[
   \frac{9}{17} - \frac{7/75}{617/1000} = \frac{11899}{31467} \approx 0.378142 > \frac{317}{1000} > \frac{119}{400}.
   \]

3. **Intermediate imbalance $r \in [1/4, 1/3]$ (and $r \in [3, 4]$ by symmetry):**
   On $r \in [1/4, 1/3]$, $R \ge 10/3$, so $e_2(r, t) \ge 1 - 2/R \ge 2/5 = 0.4$. Since $\alpha_{\rm ext} \ge 617/1000$ and $r \ge 1/4$, for every $m \ge 2$ one has $r \alpha_{\rm ext}/2 \ge 617/8000 = 0.077125 > 1/25 \ge w_m^{-2}$, so all $m \ge 2$ Stage-0 completions are strictly positive. Only $m = 1$ has a negative completion, bounded below by
   \[
   D_1(r, t) = \frac{p_{0,1}(r\alpha_{\rm ext}/2 - 1/9)}{\alpha_{\rm ext} + p_{0,1}/(6r)}
   \ge \frac{\frac{3}{4}(\alpha_{\rm ext}/8 - 1/9)}{\alpha_{\rm ext}},
   \]
   giving total head reserve at least
   \[
   \frac{2}{5} + \frac{3}{4}\cdot\frac{\alpha_{\rm ext,cent}/8 - 1/9}{\alpha_{\rm ext,cent}} = \frac{106229}{296160} \approx 0.358688 > \frac{317}{1000} > \frac{119}{400}.
   \]

4. **Large spectral energy $t \ge 1/2$ on $r \in [1/3, 3]$:**
   For $t \ge 1/2$, dropping the positive rank-one $\beta$ square gives the head reserve $e_0(t) = t/(1 + t) \ge 1/3$. For $r \ge 1/3$, every $m \ge 2$ completion is positive, while for $m = 1$, $p_{0,1}(t) = \frac{3}{(t+1)(t+4)} \le \frac{3}{(3/2)(9/2)} = \frac{4}{9}$ and $r\alpha_{\rm ext}/2 - 1/9 \ge \alpha_{\rm ext}/6 - 1/9 < 0$, so
   \[
   e_0(t) + D_1(r, t) \ge \frac{1}{3} + \frac{4}{9}\cdot\frac{\alpha_{\rm ext,cent}/6 - 1/9}{\alpha_{\rm ext,cent}} = \frac{16361}{49977} \approx 0.327371 > \frac{317}{1000} > \frac{119}{400}.
   \]

---

## 6. Certified coverings on $[1/3, 3] \times [0, 1/2]$ and strict onset ceilings

**CERTIFIED (`replay.py`, `independent.py`, `verify.py`, `VERIFICATION.txt`).**
On the compact residual domain $(r, t) \in [1/3, 3] \times [0, 1/2]$:
- **Primary replay (`replay.py`):** uses symmetry $r \leftrightarrow r^{-1}$ (which swaps channels $X \leftrightarrow Y$) to reduce $[1/3, 3] \times [0, 1/2]$ to $[1, 3] \times [0, 1/2]$ (area $1$). On every rectangle $[r_{\rm lo}, r_{\rm hi}] \times [t_{\rm lo}, t_{\rm hi}]$, all Stage-1 and Stage-0 completions are evaluated outward using scaled integer arithmetic ($Q = 10^{12}$ with exact floor/ceiling rounding from `fractions.Fraction`). Binary subdivision certifies:
  - `causal` ($d_{\rm caus,*} = 119/400$): **99 leaves**, maximum depth **17**, exact area $1$.
  - `centered` ($d_{\rm cent,*} = 317/1000$): **107 leaves**, maximum depth **17**, exact area $1$.
- **Independent replay (`independent.py`):** covers the **entire** $[1/3, 3] \times [0, 1/2]$ domain (area $4/3$) directly by 4-way quadtree subdivision without invoking $r \leftrightarrow r^{-1}$ symmetry, including the negative-completion outward bound \eqref{eq:D0_neg} for $Y$-head Stage-0 completions when $r < 1$. It also verifies 72 supplementary $6 \times 6$ rational trial forms (frequencies $0,1,2$ on each side, with 21 opposing normal weights retained in each side's diagonal entries) via exact symmetric Schur elimination, together with a separate indefinite $6\times6$ rejection guard, and verifies that near-tight full-head Poisson admissibility forces $b_1/b_0$ toward $1$. These trial forms are supplementary samples, not the infinite-frequency reserve certificate. Quadtree subdivision certifies:
  - `causal` ($d_{\rm caus,*} = 119/400$): **245 leaves**, maximum depth **9**, exact area $4/3$.
  - `centered` ($d_{\rm cent,*} = 317/1000$): **287 leaves**, maximum depth **9**, exact area $4/3$.
- **Tightest certified leaf slacks (primary binary covering on $[1, 3] \times [0, 1/2]$):** across all leaves and all $9 \times 9 = 81$ head sub-boxes $I_a \times I_b$, the minimum certified determinant slack $2 t_{\rm lo}(1 + t_{\rm lo})\Phi_X^{\rm lo}\Phi_Y^{\rm lo} + g_{\rm lo} - 2$ occurs on the balanced small-$t$ corner $r \approx 1, t \approx 0$ at 0-based head sub-box $(a, b) = (6, 5)$ ($I_6 = [23/25, 47/50]$, $I_5 = [9/10, 23/25]$):

  | Class and reserve | Tightest leaf $r$ interval | Tightest leaf $t$ interval | Worst head box $(a, b)$ | Certified lower slack $2 t_{\rm lo}(1 + t_{\rm lo})\Phi_X^{\rm lo}\Phi_Y^{\rm lo} + g_{\rm lo} - 2$ |
  |---|---|---|---|---|
  | Centered, $d_{\rm cent,*} = 317/1000$ | $[271/256,\ 17/16]$ | $[0,\ 1/256]$ | $(6, 5)$ | $0.000233258273 > 0$ |
  | Centered, $d_{\rm cent} = 19/60$ | $[139/128,\ 9/8]$ | $[0,\ 1/128]$ | $(6, 5)$ | $0.000201053907 > 0$ |
  | Causal, $d_{\rm caus,*} = 119/400$ | $[265/256,\ 67/64]$ | $[0,\ 1/128]$ | $(6, 5)$ | $0.000283795384 > 0$ |
  | Causal, $d_{\rm caus} = 19/64$ | $[133/128,\ 35/32]$ | $[0,\ 1/64]$ | $(6, 5)$ | $0.000180823397 > 0$ |
To convert the certified reserves into strict onset inequalities without requiring a uniform-in-temperature row bound, we report the strictly smaller reserves (`exact_budgets_and_ceilings` in `Finite.lean`):
\[
d_{\rm cent} = \frac{19}{60} < d_{\rm cent,*} = \frac{317}{1000},
\qquad
2(1 - d_{\rm cent}) = \frac{41}{30},
\qquad
1 - \frac{1 - d_{\rm cent,*}}{1 - d_{\rm cent}} = \frac{1}{2050} > 0,
\]
\[
d_{\rm caus} = \frac{19}{64} < d_{\rm caus,*} = \frac{119}{400},
\qquad
2(1 - d_{\rm caus}) = \frac{45}{32},
\qquad
1 - \frac{1 - d_{\rm caus,*}}{1 - d_{\rm caus}} = \frac{1}{1125} > 0.
\]
Thus at every temperature $\tau \ge C_{\rm cent}\sqrt{\rho}$ on the centered class (resp. $\tau \ge C_{\rm caus}\sqrt{\rho}$ on the causal class), the bare stability form satisfies $Q(v) \ge \frac{1}{2050}\sum_x \eta_x \|g_x\|^2$ (resp. $\frac{1}{1125}\sum_x \eta_x \|g_x\|^2$). Since $g_{n,x}=\sqrt{w_n b_{n,x}}v_{n,x}$ and $s_{n,x}\le w_n+R(\tau)\le(1+R(\tau))w_n$, the normal-coordinate norm obeys $\sum_x\eta_x\|(s_x b_x)^{1/2}v_x\|^2\le(1+R(\tau))\sum_x\eta_x\|g_x\|^2$. This converts the bare gap into a positive symmetrized operator gap $\ge [2050(1 + R(\tau))]^{-1} > 0$ (resp. $[1125(1 + R(\tau))]^{-1} > 0$), proving $\tau_{\rm lin} < C_{\rm cent}\sqrt{\rho}$ and $\tau_{\rm lin} < C_{\rm caus}\sqrt{\rho}$.

---

## 7. Initial NUMERICAL scan and what remains OPEN

**NUMERICAL (`search_finite_energy.py`, `SEARCH.json`).**
Before committing to the proof route above, we executed a fast numerical scan of 50 single-channel and two-channel finite-energy Einstein and Lorentzian glue instances at $\tau = C_\infty\sqrt{\rho}$, solving the exact nonlinear centered normal equation to relative residual $< 10^{-13}$ (`n_norm = 128`, `n_pair = 64`) across flat, centered box, symmetric shell, and notched DOS profiles:
- Single-channel rows (20 rows, already covered analytically by Zenodo 24364879): maximum pairing eigenvalue $0.999491612731$ (soft Einstein $t = 0.02$, flat DOS).
- Two-channel rows (30 rows across off-diagonal, asymmetric rank-one, and mixed two-mode Einstein+Lorentzian glue): maximum pairing eigenvalue $0.999491612731$ (off-diagonal soft Einstein $t = 0.02$, flat+flat DOS); maximum non-flat finite-energy gain $0.999450370390$ (off-diagonal soft Einstein $t = 0.02$, wide-box+narrow-box).
- Every tested finite-energy row has pairing eigenvalue strictly below $1$ at $\tau = C_\infty\sqrt{\rho}$.

**OPEN questions:**
1. Whether the sharp centered finite-channel constant $S = \sup k_B T_{\rm lin}/\sqrt{\rho}$ at finite glue energies equals $C_\infty \approx 0.18272625$ or lies strictly inside $[C_\infty, C_{\rm cent}] \approx [0.18272625, 0.18605932]$ remains **OPEN**.
2. Closing the remaining $1.824\%$ gap between $C_\infty$ and $C_{\rm cent}$ (if $S = C_\infty$) requires either exploiting the full multi-frequency self-consistency of the nonlinear normal equation $s_{n,x} = w_n + \sum_m K_{nm}^{\rm norm} b_{m,x}$ beyond the pairwise width-cone relation $w_m^{-1} \le s_{m,x}/s_{0,x} \le w_m$, or a global operator rearrangement across channels with unequal DOS profiles.

---

## 8. Verification boundary and References to verify

**Verification boundary:**
- `Finite.lean` (verified by Lean 4.35.0, exit code 0) checks the exact polynomial identities after clearing denominators (`head_square`, `decoupled_head_criterion`, `reflected_tail_decomposition`, `cascade01_width_diff`, `stage0_transfer_identity`, `stage1_transfer_identity`, `stage0_completion_numerator`, `stage1_completion_numerator`) and the exact rational budget and ceiling inequalities (`exact_budgets_and_ceilings`). It does not formalize infinite-dimensional operators, the branch-and-bound covering loop, or the analytic Poisson convexity.
- `replay.py` and `independent.py` (verified by `verify.py`, exit code 0, with 21 guarded mutants—13 in `replay.py` and 8 in `independent.py`—all exiting with code 1 and `REJECT:`) certify the finite rational coverings, closure-pin integrity, decimal constant enclosures, reflected-tail cross-term bound, first-tail Stage-1 threshold, upper-Harnack negative-completion sign split, two-sided width-cone bathtub lower bound, and indefinite 6x6 Schur-complement rejection in exact rational/integer arithmetic.
**References to verify:**
1. D. Rens, *Sharp Eliashberg pairing bounds for an energy-dependent density of states*, Zenodo (2026), DOI `10.5281/zenodo.24364879`, supplement DOI `10.5281/zenodo.24364881` (local source: `PUBLISH/finite-repo/README.md` and `POLISH2/inputs-mb/finite-proofs.tex`). Cited for the single-channel scalar finite-energy theorem $k_B T_{\rm lin} < C_\infty\sqrt{\rho}$.
2. D. Rens, *What a glue moment bounds in linear Eliashberg theory*, Zenodo (2026), DOI `10.5281/zenodo.24372930`, supplement DOI `10.5281/zenodo.24372929` (local source: `PUBLISH/univ-glue/repo/main.tex` and `RIGOR.md`). Cited for the reduced multichannel Eliashberg stability form, the soft-limit constant $C_\infty = 0.18272624777232952$, and the initial centered upper bound $U = \sqrt{6466219/3576050}/(2\pi)$.
3. D. Rens, *A $\sqrt{377/250}/(2\pi)$ moment-norm ceiling for causal multichannel Eliashberg branches*, Zenodo (2026), DOI `10.5281/zenodo.24374032`, supplement DOI `10.5281/zenodo.24374034` (local source: `PUBLISH/optimum/source/main.tex` and `UNIV-GLUE/OPTIMUM/RIGOR.md`). Cited for the single-stage causal multichannel ceiling $C_{21} = \sqrt{377/250}/(2\pi)$ and the global edge-allocation step.
4. Local predecessor artifacts pinned in `closure.json`: `UNIV-GLUE/SHARP-Q2/RIGOR.md`, `UNIV-GLUE/CENTERED-GAP/RIGOR.md`, `UNIV-GLUE/GENERAL/RIGOR.md`, `UNIV-GLUE/OPTIMUM/RIGOR.md`, `UNIV-GLUE/OPTIMUM2/RIGOR.md`, and `UNIV-GLUE/note-centered/PHYS-Q2/RIGOR.md`. The separately cited, under-review `UNIV-GLUE/OPTIMUM2/SOFT-LIMIT/RIGOR.md` specializes the predecessor angular-locking result to a fixed connected graph, finite inverse moment, and centered boxes $W_x=d_x\epsilon^{-p}$ with $p>1$. With unit channel heights its lower endpoint is $C_\infty\sqrt{(\mathbf{1}^T J\mathbf{1})/(d\rho)}$; for fixed weighted heights the numerator/denominator become the constant-vector Rayleigh quotient in $\ell^2(\eta)$. This is a common-soft lower construction, not a universal finite-energy upper theorem.
5. **[INFERENCE] Prior-work boundary:** Consistent with the predecessor notes in `UNIV-GLUE/OPTIMUM/` and `UNIV-GLUE/OPTIMUM2/`, no new external third-party full texts were fetched or used; the Cauchy--Schwarz, Schur test, Harnack, and Poisson rearrangement steps are classical mathematical tools applied to the local SuperMatics reduced Eliashberg model.
