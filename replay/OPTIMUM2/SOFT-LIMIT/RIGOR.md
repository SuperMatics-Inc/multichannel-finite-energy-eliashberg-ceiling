# Sharp common-soft limits on noncommon finite-band Eliashberg branches

## Verdict and scope

**PROVED, analytic proof below:** a fixed connected reciprocal multichannel glue, softened at fixed total moment, has an exact averaged-moment linear onset limit when its finite centered DOS boxes grow faster than the normal widths. The result gives a matching lower construction for C_infinity with genuinely different channel spectra, DOS and leading normal widths. It also proves why a connected unbalanced moment matrix cannot match the moment-norm endpoint by this particular soft construction.

**OPEN:** this does not prove the C_infinity upper bound for all finite-energy coupled energy-dependent DOS branches. A matching lower sequence establishes that a universal coefficient cannot be smaller. It does not exclude a larger optimal coefficient.

**Explicit corollary of the predecessor:** `../../note-centered/PHYS-Q2/RIGOR.md` §3, **“The stability limit and the collapsed sites”**, supplies angular locking and the limiting component moment. Its §4, **“Exact limiting coefficient for scaled centered boxes”**, supplies the matching finite-support upper tests. Here the scaled boxes grow to the flat profile beta=1, so the component formula becomes m_bar. The growing-width construction and the examples below are an explicit corollary of those labeled results, not a rediscovery of angular averaging. A complete specialized proof is included to make the lower construction independently reviewable. No novelty or priority claim is made. All temperatures are energies tau=k_B T; none is a measurement or a Hamiltonian thermodynamic transition.

## 1. Hypotheses and units

There are d>=2 angular channels with counting measure. The extension to fixed positive channel measures replaces every sum by its weighted sum, and the average below by the corresponding weighted Rayleigh quotient. We prove the counting-measure version used by the explicit examples.

Let dM_xy(e) be finite reciprocal entrywise-nonnegative moment measures on positive base boson energies e. Assume

- J_xy = integral dM_xy(e) and L_xy = integral e^(-2) dM_xy(e) are finite;
- the undirected graph of positive off-diagonal J_xy is connected;
- the same retarded interaction enters the normal and the pairing equations;
- there is no additional instantaneous pairing interaction or normal self-energy;
- the radial DOS in channel x is r_x^epsilon(xi)=1_{|xi|<W_x^epsilon}, with W_x^epsilon=d_x epsilon^(-p), fixed d_x>0 in energy units, and fixed p>1;
- the real normal shift is zero, as on the centered, reflected, even-DOS branch.

The positive measures imply that the off-diagonal graphs of J and L coincide. Zero-transfer mass is finite before softening. The measures and the graph do not vary with epsilon. Replace each base energy by sqrt(epsilon)*e without changing its moment weight. At test energy tau>0 put a=(2*pi*tau)^2 and

K_j^epsilon = integral dM(e)/(epsilon*e^2+a*j^2).

Thus K_0=L/epsilon, K_j -> J/(a*j^2) for j>=1, and 0<=K_j<=J/(a*j^2). J has units of energy squared and its angular norm rho=||J|| is the height-weighted total glue moment norm. K is dimensionless. Write w_n=2n+1 and s_{n,x}=u_{n,x}/(pi*tau). The exact box Poisson factor is

b_{n,x}=(2/pi) atan[W_x^epsilon/(pi*tau*s_{n,x})].

No second DOS-height factor multiplies J. Define P_nm=K_|n-m|-K_(n+m+1)>=0 and H_nm=K_(n+m+1). The centered normal equation is

s_{n,x}=w_n+sum_(y,m>=0) P_nm(x,y)b_{m,y}.                                      (1)

All sums below are full half-line sums, not a Matsubara truncation.

## 2. Existence, uniqueness and the wide-box limit

The frequency row telescopes to K_0+2*sum_(j=1)^n K_j. If R_x^epsilon is the row sum of L/epsilon+2*sum_(j>=1)K_j and R=max_x R_x, then the normal map preserves

w_n <= s_{n,x} <= w_n+R,
R_x <= (L*1)_x/epsilon + (4/a)(J*1)_x,                                       (2)

where sum_(j>=1)j^(-2)<=2 follows by a decreasing-function integral comparison. These bounds are finite at every positive epsilon and tau.

For any nonnegative integrable DOS, its Poisson average B(s) obeys |s B'(s)|<=B(s), by differentiating s/(s^2+xi^2). In logarithmic coordinates the normal map has absolute derivative row sum at most R/(1+R): bound the derivative numerator by the normal contribution and divide by w_n plus that contribution. The invariant log-width box, measured relative to log w_n, is complete in the supremum norm. Banach contraction gives one full positive centered solution. This argument is on the infinite frequency space. It does not infer a full solution from a finite root.

On any fixed compact tau interval inside (0,infinity), (2) implies pi*tau*s_{n,x}/W_x^epsilon ->0 uniformly for each fixed n,x, because p>1. Therefore

b_{n,x} -> 1                                                               (3)

coordinatewise and uniformly on such tau intervals. The scalar inequality 1-(2/pi)atan(W/u)=(2/pi)atan(u/W)<=2u/(pi W) gives an explicit bound if needed. Multiplying (1) by epsilon and isolating its same-frequency term gives, for each fixed n,

epsilon*s_n = L b_n + O_n(epsilon).

The omitted term is bounded by nonzero transfer rows and epsilon*w_n. Consequently

epsilon*s_{n,x} -> (L*1)_x.                                                (4)

This limit depends on the inverse moment L, not merely on J. Unequal limiting row sums mean that the normal branches remain noncommon at leading order, even though their Poisson factors tend to one.

## 3. Exact stability form and the scalar reference operator

For finite-support real v define g_{n,x}=sqrt(w_n b_{n,x})v_{n,x}. Expanding the ordered difference squares and using (1) gives the exact form

Q_epsilon(g)=||g||^2 + (1/2) sum_(x,y,n,m) P_nm(x,y)b_{n,x}b_{m,y}(v_{n,x}-v_{m,y})^2
              -2 sum_(x,y,n,m) H_nm(x,y)b_{n,x}b_{m,y}v_{n,x}v_{m,y}.         (5)

The factor 1/2 belongs to ordered pairs; the reflection factor is 2. Algebraically this is sum s_n b_n v_n^2 - sum (P+2H)_nm b_n b_m v_n v_m. The original even linearized pairing form is positively congruent to (5), by the bounded invertible multiplier sqrt(s_n/w_n) at each fixed epsilon. Thus a negative Q is an actual even linear instability, not a sign test on a surrogate normal-width operator. Nonnegative reflected pairing entries permit folding absolute values, so the leading full-frequency pairing spectral top occurs in this even sector.

For fixed epsilon the difference form is bounded: its diagonal normal cost is at most R/w_n and its cross cost is bounded by the same positive rows. Finite-support tests are dense. Write a_epsilon(tau)=inf_(||g||=1)Q_epsilon(g).

Let A be the scalar flat soft operator on l^2(N_0). Its off-diagonal entries are

A_nm=[|n-m|^(-2)+(n+m+1)^(-2)]/sqrt(w_n w_m), n!=m,

and its diagonal is

A_nn=[w_n^(-2)-2*sum_(j=1)^n j^(-2)]/w_n.

A is bounded self-adjoint and compact. The off-diagonal difference part is Hilbert–Schmidt because sum_(n>=0)1/(w_n w_(n+j))<=sum_n w_n^(-2) and sum_j j^(-4)<infinity. The Hankel part is Hilbert–Schmidt since sum_(n,m)(n+m+1)^(-4)<infinity; dropping w factors already suffices. The diagonal tends to zero. Let kappa=sup sigma(A)>0 and C_infinity=sqrt(kappa)/(2*pi). The established scalar companion computes C_infinity=0.18272624777232952 to its certified precision. This proof uses the exact spectral definition, not the decimal. The finite test v=(1,1/5,1/16), u_n=sqrt(w_n)v_n, gives <u,(A-13/10)u>=247/48000>0; replay verifies this rational identity. In particular kappa>0 and finite-support tests approximate it.

## 4. Exact limit of the full stability bottom

Put m_bar=(1/d) sum_(x,y)J_xy. We prove

lim_(epsilon->0) a_epsilon(tau)=1-kappa*m_bar/a
                               =1-C_infinity^2*m_bar/tau^2.                (6)

This is the moment average, not rho in general.

### Lower inequality

Take a sequence epsilon->0 and finite-support unit near-minimizers g^epsilon. The Hankel operators in (5), after the g transformation, converge in Hilbert–Schmidt norm to the operator with entries J_xy/[a(n+m+1)^2 sqrt(w_n w_m)]. Indeed b->1 pointwise and the squared entries are dominated by a fixed summable kernel. The same holds on a compact tau interval, or along tau_epsilon->tau_*>0.

The Hankel terms are uniformly bounded. Near-minimizers have uniformly bounded total difference cost, since Q is bounded above by any fixed finite test. Extract a weak subsequential limit g. At each fixed n and off-diagonal angular edge, the same-frequency coefficient is

P_nn=L/epsilon-K_(2n+1),

and diverges positively as epsilon->0. Its endpoint weights tend to one. Hence bounded difference cost forces v_{n,x}^epsilon-v_{n,y}^epsilon->0. Connectivity gives one limiting angular value v_n at each frequency; g_{n,x}=sqrt(w_n)v_n for all x.

Fatou's lemma retains the nonnegative off-frequency limiting squares. Compactness of the Hankel term gives convergence of its quadratic value even for weakly convergent g^epsilon. Therefore the liminf is at least 1-||g||^2 plus the limiting form on g. Set u_n=sqrt(d*w_n)v_n. Expanding the remaining squares, or telescoping their scalar normal rows, makes that form exactly

||u||^2 - (m_bar/a)<u,A u>.

Since ||u||=||g||<=1, this yields a lower bound 1-(kappa*m_bar/a)||g||^2 >=1-kappa*m_bar/a. No mass lost to high Matsubara frequencies gives a hidden negative contribution: the lost norm contributes free positive bare cost and the negative Hankel operator is compact. This proves the lower inequality in (6).

### Upper inequality

Take any scalar finite-support u and set v_{n,x}=u_n/sqrt(d*w_n) on every angular channel. The same-frequency angular squares vanish identically at every finite epsilon, despite the unequal b and s. The remaining squares converge by dominated convergence. For the normal rows connecting the finite support to its complement, only nonzero frequency transfers occur and their rows are summable. Both endpoints outside the support have zero v. The bare norm converges to ||u||^2 and the Hankel term also converges. Normalize g^epsilon by its actual finite-epsilon norm. Its Rayleigh form tends to 1-(m_bar/a)<u,A u>/||u||^2. Approximate kappa by finite scalar support to obtain the reverse inequality in (6).

The preceding lower argument works along every tau_epsilon->tau_*>0; the finite upper tests work uniformly on compact positive tau intervals. A subsequence contradiction therefore proves locally uniform convergence in tau in (6).

## 5. Largest linear onset, not merely a fixed-temperature test

Define tau_lin^epsilon=sup{tau>0:a_epsilon(tau)<0}. This is the highest normal-branch instability energy. It is not a nonlinear free-energy transition. The set is nonempty for small epsilon by (6) at any fixed tau<C_infinity*sqrt(m_bar).

A uniform high-temperature bound follows directly from (5). The normalized Hankel frequency kernel 1/[(n+m+1)^2 sqrt(w_n w_m)] has Schur norm at most 2. Apply the angular norm rho after taking absolute values and use b<=1. Then Q>=1-4rho/a=1-rho/(pi^2 tau^2) on unit vectors. In particular all tau>=sqrt(rho) are strictly stable. This bound is only a compactness aid, not the claimed sharp ceiling.
Applied to the scalar soft form itself, the same bound gives kappa<=4: its nonnegative difference cost may be dropped, and twice its Hankel norm is at most 4. Thus C_infinity<=1/pi<1, and tau_*<sqrt(rho) follows from m_bar<=rho without relying on the decimal constant.


Fix delta>0 small enough that tau_*=C_infinity*sqrt(m_bar)>delta. Equation (6) makes a_epsilon(tau_*-delta)<0 for sufficiently small epsilon. Its locally uniform convergence makes a_epsilon positive on [tau_*+delta,sqrt(rho)], and the preceding high-temperature bound excludes higher onsets. Hence

lim_(epsilon->0) tau_lin^epsilon=C_infinity*sqrt(m_bar).                    (7)

No monotonicity or uniqueness of the onset in temperature is assumed. This sup argument covers possible additional instability intervals.

Since J is symmetric with nonnegative entries, m_bar=<1,J1>/||1||^2<=rho. Equality holds precisely when the constant angular vector is a top eigenvector, equivalently when all row sums equal rho, in this connected setting. A reducible graph has one locked constant per component and (7) uses the maximum component average instead; that is why connectivity was stated.

## 6. A matching passive-glue construction with noncommon branches

Use two channels, p=2, W_1=epsilon^(-2) and W_2=2*epsilon^(-2) in a fixed energy unit. There are two base boson energies 1 and 2, with moment atoms

M_1=[[1/2,0],[0,0]],
M_2=[[1/4,1/4],[1/4,3/4]].

Both matrices are PSD and entrywise nonnegative: the second quadratic is (x+y)^2/4+y^2/2. Their total and inverse moments are

J=M_1+M_2=[[3/4,1/4],[1/4,3/4]],
L=M_1+M_2/4=[[9/16,1/16],[1/16,3/16]].

J has eigenvalues 1 and 1/2, so rho=m_bar=1; its graph is connected. L has row sums 5/8 and 1/4. Consequently the exact positive normal branch has, at every fixed Matsubara index and positive test temperature,

epsilon*s_(n,1)->5/8, epsilon*s_(n,2)->1/4.

The leading normal widths are genuinely unequal. The channel spectra are not proportional: only the first channel has the first boson. The DOS boxes differ and are finite at every epsilon. Nevertheless (7) gives

lim tau_lin^epsilon/sqrt(rho)=C_infinity.                                  (8)

This is a matching lower construction in the class with actual energy-dependent DOS and noncommon self-consistent widths, not only a flat common-branch reduction.

Passivity is explicit in DOS-height-weighted interaction units. Put

W^R(z)=M_1/(z^2-epsilon)+M_2/(z^2-4epsilon).

It is analytic in the upper half-plane with the retarded boundary prescription. On positive real frequency its negative imaginary spectral measure is pi*M_1/(2sqrt(epsilon)) delta_(sqrt(epsilon)) + pi*M_2/(4sqrt(epsilon)) delta_(2sqrt(epsilon)), PSD because the atoms are PSD. The Matsubara attraction is -W^R(i nu)=M_1/(nu^2+epsilon)+M_2/(nu^2+4epsilon), exactly the K used above. The instantaneous limit is zero. For B=-Im W^R/pi, the total moment 2*integral omega B(domega)=M_1+M_2=J. Thus the sign, factors of two and spectral normalization are fixed. The same K is used in (1) and in pairing. This is a prescribed linear Eliashberg model, not a claim of a stable underlying lattice or an exact beyond-Migdal many-body realization.
### Consistency with the released theorems

**PROVED scope check:** at each positive epsilon this family satisfies the released general centered-branch hypotheses: finitely many positive-height channels; finite, even integrable DOS fixed across test temperatures for that member; energy-independent radial vertices; reciprocal nonnegative finite positive-energy moment atoms; finite zero transfer and bounded normal rows on every compact positive-temperature interval; the full matched positive imaginary normal solution; and zero reflected real shift. The atoms are additionally PSD, so this is genuinely passive glue rather than only a positive Cooper-channel envelope. There is no extra normal term or static pairing subtraction.

It violates the *common normal branch* assumption of the earlier balanced common-branch sharp reduction: the leading widths tend to different constants by (4). It also violates an *angle-independent retarded transfer within the coupled sector*: M_1 is present only on one channel. It does not violate the hypotheses of the general UNIV-GLUE or OPTIMUM bounds, which do not require common spectra, DOS or widths. The energy DOS is not literally flat at finite epsilon; those flat-DOS sharp statements therefore do not directly describe a finite member. The sufficiently broad-box theorem at fixed glue is likewise not a uniform bandwidth estimate over this epsilon-dependent glue sequence.

The released general upper coefficients are above C_infinity, so a sequence approaching C_infinity is consistent with their strict finite-member inequalities. It is also consistent with PHYS-Q2's common-soft asymptotic upper bound. Equation (8) asserts a limit only, not equality, excess, or a sharp upper comparison at any finite epsilon. No released theorem is contradicted, and no causal shifted radial-propagator representation is inferred merely from passivity: we use the established centered class.


## 7. A failed lower construction, with an exact reason

Keep the same growing finite boxes but take a single base energy 1 and

J=L=U=[[2,1],[1,1]].

U is positive definite because its quadratic is x^2+(x+y)^2. Its angular norm is rho=(3+sqrt(5))/2, its row sums are 3 and 2, and m_bar=5/2. By (7),

lim tau_lin^epsilon/sqrt(rho)=C_infinity*sqrt[5/(3+sqrt(5))] < C_infinity.     (9)
**NUMERICAL orientation, rounded:** the exact multiplier sqrt[5/(3+sqrt(5))] is approximately 0.977. The proof of strictness uses rationals, not this rounding.

### General soft ratio, proved in this class

**PROVED, not conjectured:** under sections 1–5,

lim tau_lin^epsilon/sqrt(||J||) = C_infinity*sqrt[(1^T J 1)/(d*||J||)].       (10)

For fixed positive channel measures eta_x, use the reciprocal operator (J_eta f)_x=sum_y eta_y J_xy f_y on l^2(eta), with rho=||J_eta||. Every normal angular sum gains eta_y and every form angular sum gains eta_x eta_y. The limiting constant-vector coordinate is u_n=sqrt[(sum_x eta_x)w_n]v_n. The same proof therefore gives the exact ratio C_infinity*sqrt[(sum_xy eta_x eta_y J_xy)/(rho*sum_x eta_x)]. This is the Rayleigh quotient of the constant angular vector divided by rho. It reaches C_infinity exactly when the weighted row sums are all rho on the connected component. With disconnected components, replace the numerator-to-measure average by the largest component average; isolated vertices with no self-coupling contribute zero.

For finite scaled-box profiles beta rather than beta=1, PHYS-Q2 §4, **“Exact limiting coefficient for scaled centered boxes”**, instead gives j_D=(sum_(x,y in D) eta_x eta_y J_xy beta_x beta_y)/(sum_(x in D) eta_x beta_x) on each active component. The general ratio there is C_infinity*sqrt[max_D j_D/rho]. Neither formula is asserted outside the fixed connected common-soft hypotheses, for example for a graph degenerating with epsilon, external shifts, or a finite-energy optimum.


The ratio squared is less than 191/200: this is equivalent to sqrt(5)>427/191, whose squared gap is 5*191^2-427^2=76>0. Exact replay also brackets sqrt(5) between 2236/1000 and 2237/1000.

**KILLED as a matching construction, not as a physical model:** keeping this connected unbalanced moment matrix fixed and merely making every box wider while commonly softening cannot reach C_infinity*sqrt(rho). The zero-transfer angular Laplacian penalizes its nonconstant Perron direction with a diverging cost. The allowed limiting angular direction is constant, and its Rayleigh quotient is 5/2, strictly below rho. The Perron moment bound is an upper comparison here, not the realized soft coefficient. The exact same distinction explains why some rank-one or diagonal-heavy coupled numerical families have not supplied a matching lower point. It does not classify all finite-energy parameter choices in those families.

The prior fixed-scaled-box family W_x=D_x/epsilon also fails to match the endpoint for fixed finite D_x on any positive connected component: its positive limiting weights are beta_x<1 and its exact effective moment is (sum_xy J_xy beta_x beta_y)/(sum_x beta_x)<rho. This last assertion is the already proved PHYS-Q2 section 4 result, not a new inference from our numerical searches.

## 8. Replay, formalization and scientific boundary

`replay.py` and `independent.py` check exact matrix identities, the ordered-pair/reflection cancellation on rational finite forms, the scalar rational Ritz lower test, moment and inverse-moment rows, rational spectral brackets and the strict failed-construction gap. Guarded sign and coefficient mutants must be rejected. `Finite.lean` checks the finite polynomial and rational arithmetic core. Neither script nor Lean certifies an infinite cutoff by sampling; sections 2–5 supply that analytic argument. The numerical DOS searches elsewhere remain NUMERICAL and are not used to prove (6)–(9).

The only imported scalar value is the established spectral identity C_infinity=sqrt(kappa)/(2*pi). The finite companion proof and the earlier PHYS-Q2 proof are recorded as external analytic inputs, not silently re-certified. A separate independent review is required before any release of this note. No OPTIMUM byte is modified.

## 9. Prior work and items to verify

**[INFERENCE], no novelty claim:** the internal PHYS-Q2 proof is the direct mathematical predecessor and already contains angular locking and the exact scaled-box averaged coefficient. The validated paperpipe extraction for Scheurer, DOI 10.1103/PhysRevB.93.174509, concerns the positive angle-resolved linearized instability; it does not validate novelty of this soft limit. The cached Allen-Dynes extraction, DOI 10.1103/PhysRevB.12.905, is abstract-only and cannot by itself verify the zero-Coulomb coefficient or a multiband average theorem. Parent-provided caches for DOI 10.1103/PhysRevB.108.214508 and DOI 10.3367/UFNe.2021.05.039007 have parsed text but no validated results.json in the inspected directories. No claim about their detailed contents is used here.

Prior work to verify before any public manuscript: multiband strong-coupling angular averaging and its hypotheses; the full 2023 three-band paper; the 2021 bounds review; GSI Theorems 7–10; and literature from the last 60 days. No Golubov or Dolgov statement is cited or attributed. Main owns the new paper-reader intake and novelty review. This note is a local proof artifact, not a public priority assertion.

**[INFERENCE], targeted recent metadata search, not a full-text novelty clearance:** on 19 November 2026, the queries `"Eliashberg" "multiband" "soft" after:2026-09-20 before:2026-11-20` and `"Eliashberg" "strong coupling" "average" after:2026-09-20 before:2026-11-20` returned, among other metadata, the [18 November multiband uniqueness record](https://zenodo.org/records/24326099), [arXiv:2610.09656](https://arxiv.org/abs/2610.09656), and the already cited [arXiv:2611.20496](https://arxiv.org/abs/2611.20496). These hits were not read as new paper content or used in the proof. Search date filters also returned an older generic chapter, so coverage is demonstrably imperfect. Main's multiband strong-coupling paper-reader intake is pending. No theorem attribution or novelty conclusion follows from the search.
