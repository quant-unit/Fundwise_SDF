# JBES Referee Reports: Response & Decision Record (Rounds 1–3)

**Manuscript:** *Cash-Flow Pricing of Non-Traded Assets: Finite-Pool Identification, Horizon Selection, and Inference*  
**Date:** 2026-11-28 (covering Round-1 report after `bd6be4b`, Round-1 commit `a4d444e`, Round-2 commit `971680e`, and Round-3 completed final polish)  
**Archived reports in this folder:**
- Round 1: `jbes_methodological_referee_report.pdf` (8 pp.; verdict: *"Major revision before submission"*)
- Round 2: `jbes_second_round_referee_report.pdf` (4 pp.; verdict: *"A focused revision before submission"*; closed all major Round-1 objections)
- Round 3: `jbes_third_round_referee_report.pdf` (3 pp., GPT 6.2 Astra; closes main objection to $U_{v,H}$, retains VC application as written, accepts unified T1–T8 structure)

---

## 1. Executive Summary

Across two AI-generated methodological referee rounds, almost all mathematical and expositional catches were **genuinely valid** and directly improved the paper's rigor (zero re-estimation or new Monte Carlo runs were needed in either round). However, not every suggestion was accepted uncritically:

- **One major Round-1 suspicion was factually incorrect** (suspected FW double-sizing in R) and was refuted by auditing the estimation code.
- **One substantive Round-2 framing (Item 5) was half-wrong** ($(1 - h_v)^{-1}$ vintage reweighting): we proved algebraically that normalizing by $(1 - h_v)^{-1}$ *preserves* equal unit weight on the population signal and shared shock relative to L2, while *declining* the suggested unnormalized off-diagonal comparison and sharpening the text to distinguish population signal weight from finite-sample cross-fund noise scaling.
- **One Round-2 editorial preference (p. 1) was declined** (moving Theorem 1 parts T2–T6 to the Online Supplement), keeping T1–T8 unified in §2.2 so the main narrative remains self-contained.
- **Three technically true but mild/conservative Round-2 qualifications** (Items 1, 3, and 4) were adopted losslessly because they harden general theorem statements without weakening any benchmark or empirical conclusion.

---

## 2. Points We Pushed Back On, Declined, or Qualified — and Why

### 2.1 Round-2 Item 5: The $(1 - h_v)^{-1}$ Reweighting Claim *(Half-Wrong — Qualified & Declined Unnormalized Baseline)*

**Referee claim (Round 2, p. 4):** Moving from $L_v = PV_v^2$ to $U_v = (L_v - \sum_i \omega_i^2 PV_i^2)/(1 - h_v)$ (with $h_v = \sum_i \omega_i^2$) drops own-fund products *and* multiplies vintage $v$ by $(1 - h_v)^{-1}$ ($2$ for $n_v = 2$ vs. $10/9$ for $n_v = 10$ under EW, and analogously under FW). Without comparing against the unnormalized off-diagonal numerator $N_v := L_v - \sum_i \omega_i^2 PV_i^2$, the referee argued that diagonal removal and vintage reweighting cannot be separated.

**Why the framing is half-wrong:** Under the paper's decomposition $PV_i = \mu_v + \eta_v + \nu_i$ (common conditional mean $\mu_v$, shared vintage shock $\eta_v$, idiosyncratic error $\nu_i$, with $\sum_i \omega_i = 1$), exact algebraic expansion gives:

$$
\begin{aligned}
L_v &= \mathbf{1}\cdot(\mu_v + \eta_v)^2 + 2(\mu_v + \eta_v)\sum_i \omega_i\nu_i + \sum_i \omega_i^2\nu_i^2 + \sum_{i \neq j}\omega_i\omega_j \nu_i\nu_j, \\
N_v &= \mathbf{(1 - h_v)}\cdot(\mu_v + \eta_v)^2 + 2(\mu_v + \eta_v)\sum_i \omega_i(1-\omega_i)\nu_i + \sum_{i \neq j}\omega_i\omega_j \nu_i\nu_j, \\
U_v &= \mathbf{1}\cdot(\mu_v + \eta_v)^2 + 2(\mu_v + \eta_v)\sum_i \frac{\omega_i(1-\omega_i)}{1 - h_v}\nu_i + \frac{1}{1 - h_v}\sum_{i \neq j}\omega_i\omega_j \nu_i\nu_j.
\end{aligned}
$$

Taking conditional expectations given $\mathcal{G}_v$:
- $\mathbb{E}[L_v \mid \mathcal{G}_v] = \mathbf{1}\cdot\mu_v^2 + \mathbf{1}\cdot c_v + h_v\sigma_{\nu,v}^2$,
- $\mathbb{E}[N_v \mid \mathcal{G}_v] = \mathbf{(1 - h_v)}\cdot\mu_v^2 + (1 - h_v)c_v$,
- $\mathbb{E}[U_v \mid \mathcal{G}_v] = \mathbf{1}\cdot\mu_v^2 + \mathbf{1}\cdot c_v$.

**Decision and rationale:**
1. **Push back on signal reweighting & decline an unnormalized-$N_v$ estimation run:** Dividing by $(1 - h_v)$ does **not** alter relative population vintage signal weights between $L_v$ and $U_v$—both put exact weight **$1$** on $(\mu_v^2 + c_v)$ in every vintage. The unnormalized numerator $N_v$ would actually be a *worse* counterfactual because it distorts population signal weights by $(1 - h_v)$ (giving two-fund vintages half the signal weight of large vintages).
2. **Adopt the valid finite-sample qualification precisely:** In a finite realized sample (and under heterogeneous fund conditional means $\mu_{v,i} \neq \mu_{v,j}$), $(1 - h_v)^{-1}$ *does* rescale the surviving realized off-diagonal cross-fund noise $\sum_{i \neq j}\omega_i\omega_j \nu_i\nu_j$ (and heterogeneous-fund mean products), giving small/concentrated cohorts larger realized noise leverage.
3. **Manuscript edit:** Right after eq. (11) (`sdf_paper.tex` §2.2), we state explicitly that $(1 - h_v)^{-1}$ *preserves identical unit weight on the common vintage signal $\mu_v^2 + c_v$ relative to L2 while scaling realized cross-fund noise and within-vintage heterogeneity in small or concentrated cohorts*, and in §§5.2–5.3 we refer factually to $+0.435$/$+0.551$ as the *"matched 23-vintage difference between pooled L2 and the normalized pair estimator $U_{v,0}$."*

---

### 2.2 Round-2 p. 1 "Editorial Priority": Move Theorem 1 Parts T2–T6 to the Supplement *(Declined)*

**Referee suggestion (Round 2, p. 1):** *"T1 and T7–T8 can now carry the main narrative. Moving the intermediate decompositions to the supplement would give the identification assumptions more room."*

**Decision and rationale:** **Declined.** Keeping T1–T8 unified in Theorem 1 (§2.2) preserves the complete, self-contained logical bridge explaining *why* L2 fails at inception (T3–T4), *why* horizon averaging merely rescales and shifts rather than removes noise bias (T1–T2, T5–T6), and *why* $U_{v,H}$ solves the bounded-$n_v$ problem under $c_v = 0$ (T7) while retaining a shared-shock score term under $c_v > 0$ (T8). Moving T2–T6 out would sever direct in-text links used by Table 1, Figure 2, §4.2, and the Fix-$A$ decomposition in §5.2, and the full unified Theorem 1 already fits cleanly inside the 40-page JBES cap.

---

### 2.3 Round-1 Suspected FW Double-Sizing Bug in R Code *(Incorrect Suspicion — Refuted by Code Audit)*

**Referee concern (Round 1, p. 3):** Because Appendix A Table A.1 previously contained contradictory notation (writing $CF_{v,i,t}$ in dollars in one row while multiplying by $s_{v,i}$ in another), the Round-1 referee flagged a potential double-sizing bug in our FW estimation code that would require re-estimating all FW results.

**Decision and rationale:** **Refuted; zero re-estimation.** Direct audit of `r_project/empirical/prep_preqin.R` (L294–305), `r_project/helper/pool_vintage.R` (L37–66), and `r_project/helper/objective.R` (L137–146) verified that all fund cash flows are normalized per unit of commitment (`Size = 10` standardization in Preqin) and weighted once by $\omega_{v,i} = s_{v,i}/S_v$ in FW. We fixed only the LaTeX notation in Table A.1, §2.2, §3.4, and table notes—and in Round 2 additionally restricted the exact DLP (2012) inverse-FoF-size dollar-portfolio equivalence strictly to FW (noting EW averages per-dollar fund NPVs equally).

---

### 2.4 Technically Valid but Mild / Conservative Round-2 Points *(Adopted Losslessly)*

Three Round-2 points were mathematically true in full generality, even though they did not threaten our main benchmarks. Because adopting them cleanly cost zero substantive concessions and closed off any pedantic referee objections, we incorporated all three:

1. **Item 3 — Stylized $12.65$~yr duration root vs. additive-$\beta$ MC calibration ($-0.008$ bias):** In the constant-return positive-drift corollary (Prop. 2(i)), $\bar F/G$ cancels and $\alpha$ and $\beta$ share the exact same duration equation $\bar\tau(H_0^\ast; r) = T_{\mathrm{eff}}$—explaining why $12.65$~yr and the simulated additive-$\beta$ zero-crossing ($\approx 12.5$~yr) nearly coincide. Under stochastic returns, however, $S_{v,t}^{(\beta)} = \sum_u F_u/G_u$ is not proportional to $t$, so distinguishing the analytical $12.65$~yr as *motivation* and the VC shared-shock Monte Carlo (`u_vc_shared`, $-0.008$ bias) as *calibration evidence* is strictly more transparent.
2. **Item 4 — Local vs. global optimum in eq. (17) + uniform differentiability:** In our constant-return benchmark (Prop. 2(i)), $\bar\tau(H; r)$ is strictly increasing so $H_0^\ast$ is unique; for an arbitrary smooth $B_0(H)$ with multiple simple zeros (e.g. $B_0(H) = (H-1)(H-2)$), Taylor expansion around one zero yields only a *local* MSE optimum. Stating eq. (17) as a local optimum near $H_0^\ast$ (and global when $H_0^\ast$ is the unique globally separated zero on a compact horizon domain) plus stating uniform differentiability of the $O(V^{-1})$ expansion hardens Prop. 2(iii) at zero cost.
3. **Item 1 — Calling diffuse infill "sufficient" rather than "necessary" (Intro):** In our positive-duration, positive-drift deal model (T6), the $H=0$ variance gradient is strictly nonzero, so diffuse infill is genuinely needed at $H=0$ there; in full generality, however, a variance function could have zero derivative at $\theta_0$ or cancel at $H = H_0^\ast$. Describing diffuse infill ($\sum_i \omega_{v,i}^2 \to 0$) as a *sufficient* route to L2 consistency absent shared shocks is strictly accurate and lossless.

---

### 2.5 Double-Blind Attribution Policy (`\citet{TP24}`) *(Followed with Strict Third-Person Neutrality)*

Round 1 rightly asked for explicit attribution to Tausch & Pietz (2024, eq. 4) where eq. (8) defines horizon averaging. We added standard third-person bibliographic citation (`\citet[eq.~4]{TP24}`) with zero self-identifying phrasing, preserving full JBES double-blind compliance while giving proper scholarly credit (signed off in Round 2, p. 1).

---

## 3. Summary of Points Followed (Genuine Errors, Slips, and Framing Improvements)

### 3.1 Round 1 (Committed in `a4d444e` — All Signed Off in Round 2 p. 1)

| Area | Genuine Issue Identified | Fix Applied |
| :--- | :--- | :--- |
| **Prop. 2(iii) efficiency** | False Cauchy–Schwarz $H=0$ variance-minimality proof dropped the random product-rule term $2A_{v,H}(\nabla_\theta A_{v,H})U_{v,0}$ and omitted $2B_0'b_1$ in actual finite-$V$ MSE. | Replaced with exact sample-score identity, verified two-fund counterexample ($V_0 = \sigma^2/2 + \sigma^4 > V_{\mathrm{later}} = \sigma^2/2 + \sigma^4/4$), actual-MSE numerator $\mathcal{V}' + 2B_0'b_1$, and empirical core-MC RMSE framing. |
| **Prop. 2(i)–(ii) duration & calendar shocks** | Needed symmetric ordered-pair $\Sigma_v$, $T_{\mathrm{eff}}>0$ sufficiency (signed-duration counterexample), cap attainability, score-root vs. global-min distinction, stylized multiplicative-shift scoping for Brownian shrinkage, and exact exit-overlap covariance condition (with late-short-deal counterexample). | Incorporated all conditions, counterexamples, and scoping across §3.1 and Appendix C.11(i)–(ii). |
| **Theorem 1(T7), Lemma C.1, A2, Merton bound** | T7 needed sign-aligned proportionality (opposing-slope counterexample), local spanning, and global separation; Lemma C.1 needed single-origin-caller $U_{v,0}\to 0$ and $2\tau^\ast > s_v^{\min}$ coercivity; A2 needed $n_v=1$ allowed for linear moments; Merton $\sigma/\sqrt{T}$ was not a formal cash-flow $\alpha$ lower bound. | Fixed T7, Lemma C.1, A2, Condition (P)/(X), and removed Merton lower-bound claim from main text (kept as informal heuristic in Table F.8). |
| **Empirical interpretation, attribution & structure** | Needed "$\alpha$ fixed at zero" everywhere, VC sensitivity pair instead of "band"/bracketing, matched 23-vintage L2–$U$ contrast alongside DLP 2.73 caveat, 27 vs. 23 VC vintage distinction, AR vs. Wald clarity, T8 centering & fixed-$b$ finite-$V$ approximation caveats, two-fund lead-in, and Table 1 comparison. | Implemented across Abstract (184 words), main text (40 pp.), Supplement (67 pp.), and R table-note scripts (zero numerical changes). |

### 3.2 Round 2 (Focused Revision After `a4d444e`)

| Round-2 Item | Genuine Slip / Issue Identified | Fix Applied |
| :--- | :--- | :--- |
| **Item 1 (Table 1, Lemma C.2, Prop. 2(iv), App. B)** | Table 1 overstated all 3 targets (L2 "varies with $\theta$" contradicted Prop. 2(i); $U_{v,H}$ "$c_v>0$" ignored derivative and omitted T7 global separation; AR row overstated realized coverage and guaranteed scalar non-uniqueness); Lemma C.2 lacked Herfindahl diffusion $\sum \omega_i^2 \to 0$; Prop. 2(iv) ignored $A(\beta)$ moving/removing $\beta^\ast$ at $H>0$; App. B lacked T7 reference. | Updated Table 1 (all 3 rows), T3, Intro, Conclusion, Lemma C.2, Table B.1, Prop. 2(iv), Table C.1, and Appendix B (`sdf_paper_appendix_sdf.tex`). |
| **Item 2 (§2.2 example & §C.4 T3 proof)** | Return unit ambiguity ($R^\dagger = R_0(1+\sigma^2/n)^{1/T}$ per-period vs. holding-period); T3 proof still cited old $n=1$ coefficient $2\sigma_N^2$ instead of $\sigma_N^2$; conditioning on $\mathcal{F}_T^{\mathrm{pub}}$ failed when $K=2T+1 > T$. | Stated both return conventions, conditioned on $\mathcal{F}^{\mathrm{pub}}$ ($\mathcal{F}_{\max(T,K-1)}^{\mathrm{pub}}$), and corrected T3 proof to $\sigma_N^2$ (general $2\sigma_N^2/n$). |
| **Item 3 (Figure 1B, Tables C.2–C.3, Prop. 2(ii))** | Stale *"Zero-bias horizon"* axis in Figure 1B and captions of Tables C.2–C.3; $12.65$~yr motivation vs. $-0.008$ additive-$\beta$ MC calibration conflated; Prop. 2(ii) needed i.i.d.-$e_u$/equal-outlay scoping for raw overlap length, nonnegative vs. strictly positive common-entry covariance, and linear-component scoping in $\Lambda_H^{\mathrm{cross}}$. | Rebuilt Figure 1 (`figure_optimal_H_theory.tex` L126 -> `Zero-score horizon`), updated Tables C.2–C.3 captions & §D.4, separated motivation from calibration evidence, and tightened Prop. 2(ii) (main & §C.11). |
| **Item 4 (Prop. 2(iii) & §C.11)** | Missing $q_v$ in $2q_v A_{v,H}(\nabla_\theta A_{v,H})U_{v,0}$ (Prop. 2(iii) & eq. (C.3)); random $(\sum q_v)^{-1}$ pulled outside unconditional $\mathbb{E}[\cdot]$ in §C.11(ii)(a); local vs. global optimum & uniform differentiability in eq. (17). | Restored $q_v$, conditioned on $\mathcal{W}$ (or kept $\tilde q_v$ inside $\mathbb{E}[\cdot]$), and stated eq. (17) as a local optimum with uniform differentiability. |
| **Item 5 (Matched 23-vintage VC L2–$U$)** | Moving from $L_v$ to $U_v$ combines diagonal removal with $(1-h_v)^{-1}$ normalization. | Qualified eq. (11) (accurately noting preservation of unit weight on $\mu_v^2+c_v$ alongside $(1-h_v)^{-1}$ noise/heterogeneity scaling) and §§5.2–5.3 phrasing; declined unnormalized-$N_v$ run (see §2.1 above). |
| **Item 6 (Theorem alignment)** | T4 wrote $\mu_v = \mathbb{E}[PV_i\mid\mathcal{G}_v]$ without assuming equal fund means; (X) said "requires" instead of "sufficient"; (P) needed selected-cohort note for future-reporting filter; §5.3 overstated coercivity & Remark C.1(iii) conflated parameter-dependent $A_{v,H}(\alpha,\beta)$ with eq. (18)'s fixed $z_{v,H}$; Lemma C.3 follow-up lacked finite-dependence condition (AR(1) never hits zero); §3.4 & table notes didn't restrict DLP inverse-size identity strictly to FW. | Defined pool mean $\mu_v = \mathbb{E}[PV_v\mid\mathcal{G}_v]$ in T4/§C.5; fixed (X), (P), §5.3, Conclusion, Fig. C.1, Fig. F.2 (`step5_diagnostics.R`), Remark C.1(iii), §3.2, post-Lemma C.3, §2.2, §3.4, Table A.1, and `main_exhibits_u.R`. |

---

## 4. Final Build Verification & Archive State

- **Main blind PDF:** `2025 sdf paper/tex/jbes_main_blind.pdf` — **40 pages** (double-spaced, $\le 40$ JBES limit), **3 figures + 5 tables**, ~183-word abstract ($\le 200$), **0 LaTeX/BibTeX warnings**, **0 overfull hboxes >10pt**, clean blind scan.
- **Supplement blind PDF:** `2025 sdf paper/tex/jbes_supplement_blind.pdf` — **67 pages**, **0 overfull hboxes >10pt**, clean blind scan.
- **Zero numerical drift:** `r_project/results/figures_bovc_us_g1/main_u/exhibits_data.csv` is bit-identical before and after regenerating table notes; `2025 sdf paper/jbes/jbes_replication.zip` rebuilt with clean leak scan.
- **Archived files (`2025 sdf paper/jbes/archive/`):**
  - `jbes_methodological_referee_report.pdf` (Round 1)
  - `jbes_second_round_referee_report.pdf` (Round 2)
  - `jbes_referee_response_decisions.md` (this record)

---

## 5. Round-3 Independent Verification & Implemented Decisions (2026-11-28, after `971680e`)

**Verdict on Round-3 Referee Report (GPT 6.2 Astra):** Every single mathematical claim, counterexample, and bookkeeping catch in the 3-page Round-3 report was independently checked from first principles against `sdf_paper.tex` and `sdf_paper_online_appendices.tex`. **Astra is 100% mathematically right on all points.** Together with the author's editorial decision on fixed-$b$ positioning, all Round-3 edits were implemented cleanly (zero re-estimation, zero R/MC changes; completed checklist in `../jbes_referee_revision_plan.md` §6).

### 5.1 Author Editorial Decision — Fixed-$b$ Inference Positioning (Abstract, Keywords, Cover Letter)

- **Context:** Fixed-$b$ Bartlett inference (Kiefer & Vogelsang 2005; Lazarus, Lewis, Stock, & Watson 2018) is a standard off-the-shelf finite-sample HAC correction (analogous to Student-$t$ vs. Normal critical values when bandwidth ratio $b = S/V \approx 0.26\text{--}0.30$ is non-negligible at $V = 23\text{--}27$).
- **Implemented decision:**
  1. **Keep intact in the body (§3.3, §4, §5, Appendix D):** Essential methodological craft that preempts the main JBES/finance referee objection about small-$V$ overlapping vintages (preventing 13.3%–14.1% naive-Normal over-rejection).
  2. **Drop `Fixed-$b$ inference` from Keywords (`sdf_paper.tex` L82):** Keywords route the paper to referees. Advertising `Fixed-$b$ inference` invites routing to HAR-asymptotics theorists (who would note Corollary 1 applies KV05 rather than inventing new fixed-$b$ theory, and may fixate on trending $n_v$ / A5 homogeneity) rather than financial econometrics, asset-pricing, or panel-GMM referees. Replaced with `Horizon selection` (matching the title and Proposition 2).
  3. **Revise Abstract (`sdf_paper.tex` L77) and Cover Letter (`jbes_cover_letter.tex/.md`):** Replaced *"and develop one-dimensional fixed-$b$ and identification-robust inference"* with *"and establish valid small-sample (fixed-$b$) and identification-robust inference"* (preserving the JBES signal without overselling fixed-$b$ as our invention). Updated the cover letter's editorial-track sentence to focus on nonlinear financial econometrics, asset-pricing identification, horizon selection, and identification-robust inference.

### 5.2 Independent Algebraic Verification & Applied Fixes (Astra's Three Substantive Findings + Two Bookkeeping Items)

| Round-3 Item | Independent Verification | Applied Fix |
| :--- | :--- | :--- |
| **1. Table 1 L2 variance conditioning (p. 2)** | **100% verified.** Counterexample ($R_f=1$, $F\in\{\pm 1\}$ w.p. $1/2$, $Y\mid(F=+1)\in\{3/8,11/8\}$ w.p. $1/2$, $Y\mid(F=-1)=9/8$, $PV(\beta)=Y/(1+\beta F)-1$) satisfies (P) at $\beta=0$ ($\mathbb{E}[Y]=1$) while violating (X) ($\mu_+=-1/8, \mu_-=+1/8$). Exact derivatives at $\beta=0$: $\partial_\beta\mathbb{E}[\operatorname{Var}(PV\mid F)]|_0 = -1/4 \neq 0$, $\partial_\beta\mathbb{E}[\mu_F^2]|_0 = +1/4$, total L2 score $= 0$, curvature $L''(0) = 91/32 > 0$ (unique local recovery at $\beta=0$ despite nonzero $\mathcal{G}_v$-variance slope!). Meanwhile T3's $\partial_\beta\operatorname{Var}(PV\mid\mathcal{F}_0)|_0 = 0$ holds. Table 1 L331 had conflated $\mathcal{F}_{t_v^0}$ (T3) and $\mathcal{G}_v$ ($\sigma_{\mathrm{eff},v}^2$, T4–T6) and omitted vintage averaging. | **Adopted.** Split Table 1 L2 population-target condition cleanly: at $H=0$ under (P), vintage-avg. $\nabla_\theta\mathbb{E}[\operatorname{Var}(PV_v\mid\mathcal{F}_{t_v^0})]|_{\theta_0} \neq 0$ (T3); at any $H \ge 0$ under (X) & T4–T6, vintage-avg. $\nabla_\theta\mathbb{E}[A_{v,H}^2\sigma_{\mathrm{eff},v}^2]|_{\theta_0} \neq 0$. Added explicit $\mathcal{F}_{t_v^0}$ vs. $\mathcal{G}_v$ distinction in Supplement §C.4 (L162) and vintage-average phrasing in T3 (L290). Changed Abstract (L75), Conclusion (L657), and Cover Letter to *"can shift the population target."* |
| **2. Prop. 2(ii)(b) $\gamma_e \ge 0$ condition (p. 2)** | **100% verified.** Counterexample (flat returns, common entry $a=0$, equal outlays, $b\in\{1,2\}$ w.p. $1/2$, AR(1) calendar shock with $\gamma_e(0)=1, \rho=-0.8$) gives $\Gamma_e(1,1)=1, \Gamma_e(1,2)=0.2, \Gamma_e(2,2)=0.4 \implies \bar{O}(1)=0.6 > \bar{O}(2)=0.3 \implies T_{\mathrm{eff}}^{\mathrm{cal}} = (0.6 + 2\cdot 0.3)/(0.6 + 0.3) = 4/3 < 3/2 = T_{\mathrm{eff}}^{\mathrm{static}}$ (exact covariance check: $-0.075/0.45 = -1/6$). Supplement §C.11(ii)(b) L433 already states `(with $\gamma_e \ge 0$)`; Main L408 dropped it. | **Adopted.** Inserted *"with nonnegative calendar autocovariances"* into the common-entry parenthesis in Main Prop. 2(ii)(b) (L408). Retained calibration results unchanged. |
| **3. Finite-vintage $O(1/V)$ bias centering & root ties (p. 3)** | **100% verified.** Counterexample ($H=0$, $n_v=2$, identical payoff $1+\eta_v$, $R_0=1$, 1-dependent $\eta_v = \sqrt{c_v/2}(e_v+e_{v-1})$ with $c_v = c + d/\sqrt{v}$) satisfies ALL of A1–A5 (including A5 partial-sum homogeneity since $V^{-1}\sum_{v=1}^{\lfloor\tau V\rfloor} v^{-1/2} \to 0$). Yet finite-$V$ target $R_V^\dagger = 1 + V^{-1}\sum c_v$ drifts from $R^\dagger = 1 + c$ at rate $R_V^\dagger - R^\dagger \sim 2d/\sqrt{V} \gg O(1/V)$, while sampling bias $\mathbb{E}[\hat R_V - R_V^\dagger] = O(1/V)$. Moreover, even if $V^{-1}\sum \mathbb{E}[s_{v,H}(\theta_H^\dagger)] = \gamma_1(H)/V + o(V^{-1})$, the deterministic $-\mathcal{H}_H^{-1}\gamma_1(H)/V$ term enters the $1/V$ bias coefficient. Verified also that the signs in Appendix C.11(ii)(c) L453 ($+\mathcal{H}^{-2}C_{hs} - \frac{1}{2}\mathcal{H}^{-3}\mathcal{D}\Lambda$) are exact when centered on $\theta_{V,H}^\dagger$, and verified citation: **Bao & Ullah (2007), *Journal of Econometrics*, 140(2), 650–669** (time-series second-order bias/MSE, distinct from spatial `bao2007finite` in 137(2)). | **Adopted.** Added `bao2007timeseries` to `sdf_paper_refs.bib`. Centered the stochastic $O(1/V)$ expansion in §C.11(ii)(c) on $\theta_{V,H}^\dagger$ plus $(\theta_{V,H}^\dagger - \theta_H^\dagger)$, stated the explicit matching rate condition $V^{-1}\sum_v \mathbb{E}[s_{v,H}(\theta_H^\dagger)] = \gamma_1(H)/V + o(V^{-1})$ (with $-\mathcal{H}_H^{-1}\gamma_1(H)$ in $b_1(H)$, vanishing under $o(V^{-1})$) plus higher-order regularity across Prop. 2(ii)(c), Prop. 2(iii), and §C.11(ii)(c)–(iii). Added *"ties unresolved at $O(V^{-1})$"* to competing-root ranking in Prop. 2(iii) and §C.11(iii). |
| **4a. Normalization bookkeeping (p. 3)** | **100% verified.** Eq. (10) (L259) divides $\widehat{U}_H$ by $\sum_{v:n_v\ge 2} q_v$, while §3.1 (L371) and §3.3 (L463) divide by $V$ (if $q_v \equiv 2$, $Q_H$ is twice the limit of $\widehat{U}_H$, contradicting §C.11(ii)(a) L419 $\widehat{U}_H \xrightarrow{p} Q_H$). Also §C.11(ii)(a) L417 places random eligibility $\sum_{v:n_v\ge 2}$ outside an unconditional expectation, and L463 omits $q_v$ from L2 while L371 includes it. | **Adopted.** Harmonized L371/L463 and §C.11(ii)(a), denoting the limiting mean weight by $\mu_q := \lim_{V\to\infty} V^{-1}\sum_v\mathbb{E}[q_v] > 0$ (avoiding symbol clash with A2's upper bound $\bar{q}$) so $\widehat{U}_H \xrightarrow{p} Q_H/\mu_q$ (same minimizer; $\mu_q=1$ under $q_v\equiv 1$), aligned $q_v$ in L463, and moved $E_v = \mathbf{1}_{\{n_v\ge 2\}}$ and the denominator inside the unconditional expectation in §C.11(ii)(a). |
| **4b. Fund-size weighting & "can fail" bookkeeping (p. 3)** | **100% verified.** Main L109, L319, L407 (and Supp Table C.1) still wrote `"$O(1/n_v)$ idiosyncratic noise"` or `"$\sigma_{\nu,v}^2/n_v$"` rather than general $h_v\sigma_{\nu,v}^2$ (which equals $\sigma_{\nu,v}^2/n_v$ under EW). Main L484 said per-vintage centering `"fails"` for L2 when $n_v$ varies, whereas varying $n_v$ permits—rather than forces—nonzero vintage score means. | **Adopted.** Used `"idiosyncratic diagonal contribution"` ($h_v\sigma_{\nu,v}^2$, equal to $\sigma_{\nu,v}^2/n_v$ under EW) at Main L109, L319, L407 and Table C.1; changed `"fails"` to `"can fail"` at Main L484 and softened Remark D.1 L745 to `"are generally nonzero and need only average to zero"`. |

### 5.3 Round-3 Final Build Verification & Archive State

- **Main blind PDF:** `2025 sdf paper/tex/jbes_main_blind.pdf` — **40 pages** (double-spaced, $\le 40$ JBES cap; 32 lines on p. 40), **3 figures + 5 tables**, **186-word abstract** ($\le 200$), **0 LaTeX/BibTeX warnings**, **0 overfull hboxes >10pt**, **0 `??`**, clean blind scan (only third-person `\citet{TP24}` bibliography/citations).
- **Supplement blind PDF:** `2025 sdf paper/tex/jbes_supplement_blind.pdf` — **68 pages**, **0 overfull hboxes >10pt**, **0 `??`**, clean blind scan; Bao & Ullah (2007, *JoE* 140(2), 650–669) verified in bibliography.
- **Cover letter:** `2025 sdf paper/jbes/jbes_cover_letter.pdf` — **2 pages**, synced to 40 pp. / 186-word abstract / 68 pp. supplement.
- **Zero R/numerical drift:** No R scripts, CSVs, or replication code modified.
- **Archived Round-3 report:** `2025 sdf paper/jbes/archive/jbes_third_round_referee_report.pdf`.
