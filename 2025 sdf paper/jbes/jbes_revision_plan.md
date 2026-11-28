# Detailed Revision Plan for JBES Submission

*Date:* November 26, 2026  
*Target Venue:* *Journal of Business & Economic Statistics* (JBES — ASA / Taylor & Francis)  
*Editorial Alignment:* Joint Editors Yingying Fan, Michal Kolesár, **Dacheng Xiu** (Chicago Booth — financial econometrics & empirical asset pricing); Finance AE **Andrew Patton** (Duke).  
*Companion Analysis:* `jbes_fit_analysis.md` (same folder)  
*Source Files:* `2025 sdf paper/tex/sdf_paper.tex`, `sdf_paper_online_appendices.tex`, `sdf_paper_appendix_sdf.tex`, `sdf_paper_refs.bib`, `r_project/`

> **Status (updated 2026-11-27):** Sections 0–7 are the original work plan of 2026-11-26. Section 8 is the running log of what has been done since, what was decided and why, and what is still open. Read Section 8 first when resuming.

---

## 0. Strategic Targets & Hard Constraints

| Metric | Current State (`sdf_paper.tex`) | JBES Target |
| :--- | :--- | :--- |
| **Double-spaced main manuscript** | ~58–62 pp. (37 pp. text + 4 pp. refs single-spaced) | **34–38 pp.** inclusive of tables, figures, references (strictly **≤40 pp.**) |
| **Main-text floats** | 14 (7 tables + 7 figures) | **6–7** (3–4 tables + 3 figures) |
| **Online Supplementary Material** | 39 pp. standalone (`sdf_paper_appendix.pdf`) | **~45–50 pp.** standalone PDF (all proofs, primitive verification lemmas, moved floats, extra MC & empirics) |
| **Abstract** | 303 words, contains `\citet{DLP12}` & math | **160–195 words**, plain text, zero citations, minimal notation |
| **Section sequence** | Intro → Theory → Data/Imputation/NAVs + MC → Empirics → Conclusion | **Intro (+ Notation) → Theory (ID) → Theory (Optimal $H$ & Inference) → Monte Carlo → Empirical Application → Conclusion** |
| **Title convention** | Contains "SDF" acronym | **Descriptive methodological title, NO "SDF"** |
| **Blind review & replication** | Contains author GitHub URL (`quant-unit`, line 812) | **100% anonymized** main + supplement PDFs; public R replication archive + Editor Preqin exemption letter |

---

## 1. Title Strategy (Explicitly Avoiding "SDF")

### Why Dropping "SDF" from the Title Makes Strong Sense for JBES
1. **Eliminates the KN16 vs. DLP12 semantic trap at the door:** Putting "SDF" in the title invites finance referees to ask immediately whether the kernel prices public equities and Treasury bills (the Korteweg–Nagel 2016 paradigm). Leading with **cash-flow pricing / discounting models** avoids that conflation while Section 2.1 and Appendix B still prove cleanly when a reciprocal benchmark deflator constitutes an asset-class-specific SDF on the private-market payoff span.
2. **Broadens econometric appeal across JBES editors and referees:** "SDF" reads as narrow asset-pricing jargon to general econometricians (e.g., M-estimation, panel, or weak-ID referees). Framing around **finite-pool identification, horizon selection, and inference for pooled non-traded cash flows** speaks directly to both Dacheng Xiu's financial-econometrics pool and JBES's broader methodological audience.
3. **Matches JBES house style:** None of the 10 JBES benchmark papers uses domain-specific acronyms in the title (even Massacci 2025 uses *"Instability of Factor Strength in Asset Returns"*).
4. **Preserve discoverability in metadata:** Retain `Stochastic discount factor` and `Benchmark deflator` in the **Keywords** and Section 1–2 text so finance readers still find and cite the paper cleanly.

### Candidate JBES Titles (No "SDF")
1. *(Top pick — cleanest general econometric framing):*  
   **Identification, Pooling, and Inference in Cash-Flow Pricing Models for Non-Traded Assets**
2. *(Hafner–Herwartz–Wang 2026 style — method + explicit empirical anchor):*  
   **Identification and Inference for Pooled Non-Traded Cash-Flow Pricing Models, with an Application to Private Equity**
3. *(Highlights all three theoretical pillars — finite pools, horizon theory, inference):*  
   **Finite-Pool Identification, Horizon Selection, and Inference in Nonlinear Cash-Flow Discounting Models**
4. *(Tightest variation on DLP12 / Massacci style):*  
   **Finite-Pool Identification and Bias Correction in Cash-Flow Asset Pricing Estimators**

*(Final title selection deferred until WP1–WP3 edits are complete, so the title mirrors the exact finalized structure.)*

---

## 2. Work Package 1 (WP1): Asymptotic Theory Hardening & $V = 400$ Diagnosis *(Highest Priority — Gap 1)*

**Objective:** Close the three specific technical vulnerabilities a tough JBES econometrician would flag in Section 2.4 (`sdf_paper.tex:585-617`) and Appendix D (`sdf_paper_online_appendices.tex:532-570`).

### WP1.1 — Replace High-Level "Black-Box" `Assumption 2` with Primitive DGP Conditions + Verification Lemma
- **Current issue:** `Assumption 2` (`ass:asymptotics_H`) assumes directly that the nonlinear criterion $\ell_{v,H}(\theta)$, score $\nabla_\theta \ell_{v,H}(\theta)$, and Hessian $\nabla_{\theta\theta}^2 \ell_{v,H}(\theta)$ satisfy Jenish–Prucha (2009, 2012) NED/mixing and uniform $(2+\delta)$-moment envelopes, and proves Proposition 3 in 10 lines.
- **Concrete deliverables in `sdf_paper_online_appendices.tex` (Appendix D) and `sdf_paper.tex` (Sec. 3.2):**
  1. State **primitive DGP assumptions** cleanly numbered as `Assumptions A1–A5`:
     - **A1 (Compact positive-return domain & global identification):** $\Theta \subset \mathbb{R}^p$ compact with $\inf_{u,\theta\in\Theta} G_u(\theta) \ge \underline{G} > 0$, and unique population minimizer $\theta_H^\dagger \in \mathrm{int}(\Theta)$.
     - **A2 (Bounded cash-flow support & pool structure):** Deterministic or bounded fund lifetimes $t - t_v^0 \le L_{\max} < \infty$, bounded deal counts/capital weights, bounded fund counts $2 \le n_v \le \bar{n}$ (for $U_{v,H}$), bounded vintage weights $0 < \underline{q} \le q_v \le \bar{q}$, and non-vanishing minimum off-diagonal denominator $1 - \sum_i \omega_{v,i}^2 \ge \underline{d} > 0$.
     - **A3 (Weak dependence of underlying processes):** Public factor vector $F_u$ and common calendar/cohort shocks $(\eta_v, e_u)$ form a 1D $\alpha$-mixing (or $L_2$-NED on an $\alpha$-mixing base) time-series process along calendar time with mixing coefficients decaying at rate $\alpha(k) = O(k^{-r/(r-2)})$ (or finite $m$-dependence when shocks have bounded memory), while fund-idiosyncratic shocks $\nu_i$ are independent across funds and vintages given $(\mathcal{F}^{\mathrm{pub}}, \mathcal{W})$.
     - **A4 (Primitive moment & smoothness bounds):** Finite $(8 + 4\delta)$-th (or explicit $r$-th, $r > 4$) uniform moments on single-vintage compounded gross return factors $\sup_{\theta \in \Theta} \max_{0 \le k \le L_{\max} \vee 12H} (R^b_{t_v^0, t_v^0+k}(\theta) + D_{v,t_v^0+k}(\theta))$ and cash-flow innovations $e_{i,t}$, ensuring finite $(2+\delta)$ moments of the 4th-degree nonlinear products $\ell_{v,H}(\theta) \propto A_{v,H}(\theta)^2 PV_i(\theta) PV_j(\theta)$ and their first two derivatives.
     - **A5 (Nonsingular curvature & long-run variance):** Positive definite $\mathcal{H}_H = \lim_{V\to\infty} V^{-1}\sum_{v=1}^V \mathbb{E}[\nabla_{\theta\theta}^2 \ell_{v,H}(\theta_H^\dagger)]$ and $\Lambda_H = \lim_{V\to\infty} \operatorname{Var}(V^{-1/2}\sum_{v=1}^V s_{v,H}(\theta_H^\dagger))$.
  2. Prove a dedicated **Lemma D.1 (Inheritance of mixing, Lipschitz envelopes, and uniform moments)** in Appendix D:
     - Because each vintage's cash flows and horizon window span at most $L_H := \max(L_{\max}, 12H)$ months from $t_v^0$, $\ell_{v,H}(\theta)$, $s_{v,H}(\theta)$, and $\nabla_{\theta\theta}^2 \ell_{v,H}(\theta)$ are measurable functions of a **finite block** of underlying processes of fixed length $L_H$.
     - Therefore, mixing/NED of $(F_u, e_u, \eta_v)$ passes *directly* to $(\ell_{v,H}, s_{v,H}, \nabla_{\theta\theta}^2 \ell_{v,H})$ with shift at most $\lceil L_H/12 \rceil$ vintage years (and exact $\lceil L_H/12 \rceil$-dependence when factors and shocks are independent across months!).
     - Prove explicit Lipschitz dominance bounds $\|\nabla_\theta \ell_{v,H}(\theta)\| \le B_v$ and $\|\nabla_{\theta\theta}^2 \ell_{v,H}(\theta)\| \le B_v^{(2)}$ via Hölder's inequality from A1–A4, verifying Andrews (1992) / Pötscher–Prucha (1997) stochastic equicontinuity for the ULLN cleanly without hand-waving.

### WP1.2 — Reconcile 1D Vintage Calendar vs. Spatial HAC & Formalize Fixed-$b$ Inference
- **Current issue:** Mixing 2D spatial random-field citations (`JP09`, `JP12`, `KS11`) with 1D regular-lattice `KV05` fixed-$b$ critical values without a formal fixed-$b$ proposition.
- **Concrete deliverables:**
  1. Frame the vintage asymptotics cleanly on the **1D calendar-time vintage lattice** $y_v \in \mathbb{Z}$ (where overlap across vintages arises purely from overlapping active fund lifetimes along calendar time—which is 1D local time-series dependence, *not* a 2D spatial field!). Cite Kim & Sun (2011) / Conley (1999) cleanly only as the calendar-distance kernel notation $K(|y_v - y_w|/S)$ allowing occasional empty vintage years (irregularly spaced 1D calendar grid).
  2. Add **Proposition 3(iii) / Corollary 1 (Fixed-$b$ Asymptotic Distribution & AR Moment Inversion):**
     - Under A1–A5 plus a functional CLT (FCLT) for partial sums $V^{-1/2}\sum_{v=1}^{\lfloor rV \rfloor} s_{v,H}(\theta_H^\dagger) \Rightarrow \Lambda_H^{1/2} W_p(r)$ on $D[0,1]^p$ (Davidson 1994 / Wooldridge & White 1988) and bandwidth $S = \lfloor b V \rfloor$ with fixed $b \in (0,1]$, prove (in Appendix D, following Kiefer & Vogelsang 2005, Theorem 1–2) that:
       - The Bartlett long-run variance estimator $\hat\Lambda_H(b)$ using centered estimated scores $\tilde{s}_{v,H} = s_{v,H}(\hat\theta_H) - \bar{s}_H(\hat\theta_H)$ satisfies $\hat\Lambda_H(b) \Rightarrow \Lambda_H^{1/2} \Xi_p(b) \Lambda_H^{1/2}$ (with plug-in estimation error $O_p(V^{-1/2})$ negligible under the integrated-Hessian expansion);
       - Hence any linear restriction $t$-statistic converges to the pivotal Kiefer–Vogelsang fixed-$b$ distribution $t_{\mathrm{KV}}(b)$ (using the Lazarus–Lewis–Stock–Watson 2018 rule $S = \lceil 1.3\sqrt{V} \rceil$);
       - Similarly, at the true parameter $\theta_0$, the centered inception moment statistic in Equation (14) converges to $t_{\mathrm{KV}}(b)$ under Condition (P) alone, validating the identification-robust Anderson–Rubin confidence set inversion without requiring non-zero slope or unique roots (Stock & Wright 2000).
     - Explicitly note in a brief **Remark** how occasional missing vintage years (e.g., 1989 in BO, 1983–84 in VC) affect the calendar-distance kernel vs. index-distance kernel (verify empirically whether index-lag Bartlett $K(|v-w|/S)$ vs. calendar-lag Bartlett $K(|y_v-y_w|/S)$ changes any headline SE/AR interval by more than 0.01!).

### WP1.3 — Diagnose and Resolve the $V = 400$ Studentization / Heavy-Tail Puzzle (`sdf_paper.tex:803`)
- **Current issue:** Our MC currently admits that even at $V = 400$, equal-vintage $U_{v,H}$ rejects in 6.8%–10.8% of draws with $N_{\mathrm{eff}} = 13\text{--}27$ and a left-tail excess, with the untested guess *"heavy-tailed compounded weights (cf. Martin 2012) are a candidate mechanism that we have not tested."* At JBES, leaving an untested failure of nominal size at $V = 400$ is a major red flag.
- **Concrete diagnostic script in `r_project/` (estimated 1–2 hours compute):**
  1. Inspect the existing growing-$V$ simulation code and output (`r_project/simulation/` and `r_project/results/mc_step4/summary.csv` / `inference_checks.csv`) to verify **exact root cause** of the $V = 400$ low $N_{\mathrm{eff}}$ and 6.8%–10.8% rejection:
     - **Check A (Tail exponent / 4th-moment check):** Compute empirical score kurtosis, Hill tail index, and max vintage influence share $\max_v |\mathrm{IF}_v| / \sum_w |\mathrm{IF}_w|$ across $V \in \{25, 100, 400, 1600\}$ under both (i) resampled empirical market paths (which include extreme drawdown clusters like 2000–2002 / 2008–2009 where $G_u(\theta)$ drops steeply) and (ii) Gaussian/bounded-support simulated factor paths.
     - **Check B (Leverage / skewness of studentization):** Check $\mathrm{Corr}(\hat\beta - \beta_0, \widehat{\mathrm{SE}})$—for nonlinear discounting, when $\hat\beta$ undershoots $\beta_0$, $D_{v,t}(\hat\beta)$ increases, or vice versa, creating strong negative correlation between estimation error and sandwich SE (classic asymmetric studentization distortion in nonlinear M-estimators with skewed influence functions).
     - **Check C (Fix / Practical Diagnostic):** Test whether (i) $N_{\mathrm{eff}} \to \infty$ slowly as $V \to 1600$ (confirming slow-rate convergence from heavy-tailed compounding kurtosis rather than inconsistency of the CLT) OR (ii) whether a simple score-level / studentization check (e.g., reporting $N_{\mathrm{eff}}$ alongside fixed-$b$ SEs, or testing whether $H=0$ vs. $H>0$ differs because $A_{v,H}^2$ adds extra lognormal compounding tail mass!) sharply separates $H=0$ from $H=15, 30$.
  2. **Convert the finding into a clean econometric contribution:** Report the exact mechanism quantitatively (replacing *"that we have not tested"* with an actual verified mechanism + diagnostic rule: report $N_{\mathrm{eff}}$ and pair fixed-$b$ Wald inference with the self-normalized/AR check whenever $N_{\mathrm{eff}}$ is small).

---

## 3. Work Package 2 (WP2): Econometric Reframing & Literature Positioning *(Gap 2)*

**Objective:** Transform the framing from *"auditing DLP12 in private equity"* to *"general identification, finite-pool bias correction, horizon selection, and weak-ID inference for nonlinear cash-flow pricing models of pooled non-traded assets."*

### WP2.1 — Reframe Introduction & Add Notation Block (`sdf_paper.tex`, Sec. 1)
1. **Opening paragraph:** Frame the broad econometric problem first: across private equity, venture capital, private credit, infrastructure, real estate, loan/CLO portfolios, and project finance, researchers do not observe periodic traded returns—only **irregular, pooled cash-flow streams** grouped into finite vintage/cohort pools ($n_v < \infty$). Estimating nonlinear present-value discounting parameters $\theta$ by minimizing squared pooled pricing errors creates an econometric problem with no analogue in traded-return factor regressions.
2. **Sharp 4-point contribution list** (mirroring Mehrabani & Parsaeian 2025, Liu 2025, Massacci 2025):
   - *(i) Finite-pool identification failure & horizon rescaling* (Theorem 1, T1–T6): squared pooled NPV fails consistency at finite $n_v$ even at inception; horizon averaging rescales by $A_{v,H}(\theta)$ and adds a Scale-Jacobian pull.
   - *(ii) Cross-fund pair-product estimator $U_{v,H}$* (Theorem 1, T7–T8): leave-own-fund-out U-statistic criterion restores exact identification at $H=0$ for bounded $n_v \ge 2$ when $c_v = 0$.
   - *(iii) Exact optimal-horizon theory & impossibility of cross-vintage LOO* (Proposition 2): multi-deal duration matching $H_0^\ast$, calendar-overlap duration lengthening, MSE-optimal horizon $H_{\mathrm{MSE}}^\ast(V)$, and proof that cross-vintage leave-band-out pair products collapse to rank 1 and become finite-sample indefinite.
   - *(iv) Primitive increasing-domain asymptotics, fixed-$b$ inference & identification-robust AR sets* (Proposition 3 + Corollary 1), validated via 72-design Monte Carlo and illustrated on North American Buyout and VC cash flows.
3. **Add Econometric Literature Link (new paragraph in Related Literature):**
   - Explicitly connect $U_{v,H}$ to **leave-one-out / jackknife quadratic-form debiasing** and **U-statistic M-estimation**:
     - Angrist, Imbens, & Krueger (1999, *Econometrica*) — JIVE;
     - Chao, Swanson, Hausman, Newey, & Woutersen (2012, *Econometrica*) & Cattaneo, Jansson, & Newey (2018, *JASA*) — many-instrument / variance-component bias correction;
     - Kline, Saggio, & Sølvsten (2020, *Econometrica*) — leave-out estimation of variance components in two-way models;
     - Mikusheva & Sun (2022, *Review of Economic Studies*) — jackknife Anderson–Rubin weak-ID inference;
     - Honoré & Powell (1994, *Journal of Econometrics*) — pairwise-difference / U-statistic extremum estimators.
   - Articulate concisely **why pooled cash-flow discounting is fundamentally harder than linear JIVE / variance-component LOO**: (a) nonlinear parameter-dependent discounting and forward compounding $A_{v,H}(\theta)^2$, (b) non-coercive boundary collapse ($\alpha \to +\infty$) at $H=0$, (c) exact duration-matching trade-off under shared/calendar-time shocks (Prop. 2(i)–(iii)), and (d) rank-1 collapse and indefiniteness when extending LOO across vintages (Prop. 2(iv)).
4. **Add a concise `Notation.` paragraph** at the end of Section 1 (4–6 lines, matching 6 of 10 JBES benchmarks: $\|\cdot\|$, $\nabla_\theta$, $\nabla_{\theta\theta}^2$, $\xrightarrow{p}$, $\xrightarrow{d}$, $\Rightarrow$, $\mathcal{F}_t$, $\mathcal{G}_v$, $PV_v$, $\bar\epsilon_{v,H}$), pointing to Appendix A for the full glossary.

### WP2.2 — Add BibTeX Entries (`sdf_paper_refs.bib`)
Add verified entries for:
- `AIK99` (Angrist, Imbens, & Krueger 1999, *Econometrica*)
- `KSS20` (Kline, Saggio, & Sølvsten 2020, *Econometrica*)
- `MS22` (Mikusheva & Sun 2022, *ReStud*)
- `HP94` (Honoré & Powell 1994, *JoE*)
- `Andrews92` (Andrews 1992, *Econometric Theory* — generic ULLN) or `Davidson94` (Davidson 1994, *Stochastic Limit Theory*)
- `WooldridgeWhite88` (if used for mixing FCLT)

---

## 4. Work Package 3 (WP3): Reordering, Float Pruning & ≤38-Page Double-Spaced Budget *(Gap 3)*

**Objective:** Match the exact structural rhythm of the 10 JBES benchmarks and bring the double-spaced manuscript cleanly under JBES's 40-page ceiling (**target: 34–38 double-spaced pages** inclusive of references, tables, and figures).

### WP3.1 — Target Section Architecture & Page Budget (12pt Double-Spaced, 1-in. Margins)

| New Section | Content | Main-Text Floats | Target Double-Spaced Pages |
| :--- | :--- | :--- | :---: |
| **Title Page** | Separate file (JBES blind review): Title, Authors, Affiliations, Acknowledgments, Funding, Disclosure | — | *(excluded from blind MS)* |
| **Abstract + Keywords** | 160–195 words plain-text abstract, 5–6 Keywords, JEL codes | — | **1.0 pp.** |
| **1. Introduction** | General non-traded cash-flow problem, 4 numbered contributions, Econometric + PE literature (incl. LOO/JIVE & KN16/DLP12 boundary), **Notation** block, Roadmap | — | **4.5 pp.** |
| **2. Model and Finite-Pool Identification** | 2.1 Pooled cash-flow pricing (`Ass. 1`, `Lemma 1`, `Prop. 1` — **move both proofs to Online Supp**); 2.2 1-deal motivating algebra, `Theorem 1 (T1–T8)`, cross-fund $U_{v,H}$, admissible domain $\Theta$ | **Fig. 1** (`fig:exploding_alpha` *or* merge with Prop. 2 — see below) | **7.0 pp.** |
| **3. Optimal Horizon Theory and Increasing-Domain Inference** | 3.1 `Proposition 2 (i–iv)` (duration matching, cross-vintage dependence, MSE expansion, rank-1 collapse) + economic intuition; 3.2 Identification-robust AR moment diagnostic; 3.3 Primitive `Assumptions A1–A5`, `Proposition 3` (Consistency, CLT, Fixed-$b$ Corollary); 3.4 Brief positioning vs. DLP12 / KN16 / KN24 (prose; full taxonomy table in Supp B) | **Fig. 1** (`fig:optimal_H_theory`), **Tab. 1** (`tab:regime_summary`) | **7.5 pp.** |
| **4. Monte Carlo Study** *(moved BEFORE Empirics!)* | 4.1 Calibration & DGP design (concise); 4.2 Single-factor core & bias decomposition; 4.3 Joint $(\alpha,\beta)$, shared vintage shocks & $H>0$ coercivity; 4.4 Hardening summary (FW sizes, path-correlated shocks, $V=400$ studentization diagnosis, two-factor summary) | **Tab. 2** (`tab:mc_core`), **Fig. 2** (`fig:mc_decomposition`), **Tab. 3** (`tab:mc_alpha_beta`) | **6.5 pp.** |
| **5. Empirical Application: Private Equity Cash Flows** | 5.1 Data & construction (1 concise paragraph: Preqin BO/VC $V \le 2010$, constant-composition, within-vintage median FW imputation, terminal NAV haircut summary — full Sec. 3.1–3.3 details moved to Online Supp); 5.2 Single-factor Buyout & VC ($U_{v,H}$, L2, Fix-$A$, AR sets & matched-vintage gaps); 5.3 Weak identification of $(\alpha,\beta)$ (incl. VC 1993–97) and two-factor $q^5$ summary | **Fig. 3** (`fig:u_beta_profiles`), **Tab. 4** (`tab:u_l2_beta` + compact AR/matched-gap bottom panel), **Tab. 5** (`tab:u_l2_alpha_beta` — *or* keep 4 tables total by folding into compact format) | **6.5 pp.** |
| **6. Conclusion & End Matter** | Concise 2-paragraph conclusion + **Supplementary Materials** + **Data Availability** + **Disclosure** | — | **1.5 pp.** |
| **References** | Cleaned ASA/JBES bibliography (~45–55 key refs) | — | **3.5 pp.** |
| **TOTAL** | | **3 Figures + 4–5 Tables (7–8 floats max; target 6–7)** | **~37.0–38.0 pp.** *(≤40 pp. ceiling!)* |

### WP3.2 — Exact Float Disposition Plan (14 Current Floats → 6–7 Main-Text Floats)

| Current Float | Label | Action | Rationale |
| :--- | :--- | :--- | :--- |
| Fig. 1 (NPV timeline schematic) | `fig:npvs` | **Move to Online Supp** (App. A/C) | Purely pedagogical; equations (5)–(6) are crystal clear to JBES readers |
| Fig. 2 (Exploding-$\alpha$ surface) | `fig:exploding_alpha` | **Move to Online Supp** (App. C.9 alongside `lem:coercivity`) | `fig:optimal_H_theory` already carries the main theoretical visual load |
| Fig. 3 (Optimal-$H$ 3-panel figure) | `fig:optimal_H_theory` | **KEEP → Main Fig. 1** | Visualizes our crown-jewel Prop. 2 (geometry + $H_0^\ast$ curve + Regimes 1–3b RMSE) |
| Tab. 1 (Dependence regimes) | `tab:regime_summary` | **KEEP → Main Tab. 1** | High-density synthesis of Thm 1 + Prop 2 across Regimes 1–3 |
| Tab. 2 (Paradigm taxonomy) | `tab:comparison` | **Move to Online Supp** (App. B) | Keep crisp 1-paragraph distinction in Sec. 3.4; full comparison table belongs with the SDF/KN16/KN24 analysis in App. B |
| Fig. 4 (Terminal NAV haircuts) | `fig:nav_disc_BOVC_MKT` | **Move to Online Supp** (App. F) | PE-specific robustness check; summarize 0%–100% range in 2 sentences in Sec. 5.1 |
| Tab. 3 (Core MC results) | `tab:mc_core` | **KEEP → Main Tab. 2** | Core single-factor MC proof of concept |
| Fig. 5 (Bias decomposition) | `fig:mc_decomposition` | **KEEP → Main Fig. 2** | Directly illustrates Thm 1(T3–T6) target shift vs. finite-sample bias |
| Fig. 6 (MC margins M1–M7) | `fig:mc_margins` | **Move to Online Supp** (App. D/E) | Summarize one-margin-at-a-time robustness in 1 crisp paragraph in Sec. 4.2 |
| Tab. 4 (Joint $(\alpha,\beta)$ MC) | `tab:mc_alpha_beta` | **KEEP → Main Tab. 3** | Essential proof of exploding-$\alpha$ + Prop. 2 $H>0$ rescue under shared shocks |
| Fig. 7 (Empirical beta profiles) | `fig:u_beta_profiles` | **KEEP → Main Fig. 3** | Headline empirical illustration for Buyout & VC |
| Tab. 5 (Empirical single-factor) | `tab:u_l2_beta` | **KEEP → Main Tab. 4** | Headline Buyout/VC estimates + fold compact 4-row AR confidence-set panel underneath |
| Tab. 6 (IV & matched-vintage gaps) | `tab:iv_comparison` | **Merge key AR sets into Tab. 4; full table → Online Supp** (App. F) | Eliminates table duplication while preserving exact AR sets & root-fragility stats in prose/panel |
| Tab. 7 (Empirical $(\alpha,\beta)$) | `tab:u_l2_alpha_beta` | **KEEP → Main Tab. 5** *(compact)* | Anchors the weak-ID empirical takeaway (Buyout ridge + VC $+2\%$ bound) |

**Result:** **3 Figures + 5 Tables (or 4 Tables if `tab:regime_summary` is compressed) = 7–8 floats**, cutting 6–7 full pages of float space immediately while keeping every core scientific result visible in the main paper!

### WP3.3 — Prose De-Duplication Targets (~25% Word Reduction)
- **Introduction:** Cut the 4 empirical bullets (`sdf_paper.tex:132-154`, currently ~750 words repeating every decimal of Sections 3–4!) down to **1 tight empirical-illustration paragraph (~150 words)**—just like Massacci (2025), Liu (2025), and Hafner et al. (2026).
- **Section 3.1–3.3 (Data/Imputation/NAVs/CV):** Compress 4 single-spaced pages down to **~0.75 double-spaced page** in Section 5.1, moving the detailed fund-size imputation breakdown (lines 679–699), terminal-NAV share numbers (lines 713–732), and $hv$-block CV setup (lines 734–737) to Appendix F with clear forward pointers.
- **Conclusion:** Tighten from 1.5 single-spaced pages to **~0.75 double-spaced page** (2–3 crisp paragraphs summarizing methodological contributions, empirical lessons, and future directions).

---

## 5. Work Package 4 (WP4): JBES Formatting, Double-Blind Anonymization & Replication Package *(Gap 4)*

### WP4.1 — Rewrite Abstract (≤200 Words, Plain Text)
- Remove `\citet{DLP12}`, drop excessive empirical decimals and heavy math notation, lead with the general econometric problem (finite-pool bias in nonlinear cash-flow pricing models, cross-fund pair-product $U_{v,H}$, optimal-horizon duration matching & cross-vintage rank-1 collapse, increasing-domain fixed-$b$ & AR inference), and close with 1–2 crisp sentences on the Preqin Buyout/VC application.

### WP4.2 — JBES LaTeX Template & Blind-Review Hygiene
- Replace JFEC-specific `\ifSubmissionMode` (`endnotes` + `endfloat` tables-at-end) with a clean **`\ifJBESMode`** switch (12pt, 1-in. margins, `\doublespacing`, inline floats or ASA standard layout, standard footnotes, `natbib` ASA/JASA-style author-year citations via `chicago` or `apalike`).
- Maintain two completely clean, anonymized build targets:
  1. `jbes_main_blind.pdf` (no author names, no affiliations, no acknowledgments, no GitHub links, no PDF author metadata).
  2. `jbes_supplement_blind.pdf` (standalone Online Supplementary Material with `xr-hyper` cross-references synchronized to `jbes_main_blind.aux` and zero identifying info).
  3. `jbes_titlepage.pdf` (separate non-blind title page with authors, affiliations, corresponding author email, acknowledgments, funding, and conflict-of-interest disclosure).
- **Anonymization grep checklist before submission:** Search all `.tex` and `.pdf` files for `Tausch`, `Bohnert`, `AssetMetrix`, `Munich`, `quant-unit`, `github.com`, `jfinecon`, `JFEC`, `revision-2026`.

### WP4.3 — Mandatory JBES End Matter & ASA Replication Archive
1. Add standard JBES closing sections to the main manuscript:
   - **Supplementary Materials** (concise bullet list of Online Supplement Appendices A–F + Replication Code Archive).
   - **Data Availability Statement** (state clearly that U.S. $q^5$ factors are publicly downloadable from `global-q.org`, while fund-level PE cash flows are proprietary Preqin data subject to commercial licensing; provide exact Preqin export criteria and note that synthetic Preqin-calibrated cash-flow panels are provided in the replication archive to verify all code end-to-end).
   - **Disclosure Statement** (standard ASA no-conflict statement; keep blinded in blind copy).
2. **Replication Package (`r_project/` → clean `jbes_replication/` zip):**
   - Include clean, self-contained R scripts reproducing: (i) all Monte Carlo tables/figures (including the WP1.3 $V=400$ tail diagnostic), (ii) the Prop. 2 lookup tables ($H_0^\ast$ and $T_{\mathrm{eff}}$), and (iii) the full empirical estimation pipeline on a synthetic Preqin-calibrated dataset (plus exact Preqin ingestion script for licensed users).
   - Include a top-level `README.md` with R version, package versions (`renv` or `sessionInfo`), runtime table, and script-to-Table/Figure mapping.
   - *(Optional bonus, inspired by Virolainen 2025 `gmvarkit`):* Wrap the core $U_{v,H}$ / L2 / fixed-$b$ SHAC / AR-set estimator functions into a clean, documented standalone R module (or lightweight R package skeleton) included in the supplement.

---

## 6. Work Package 5 (WP5): Cover Letter & Final Pre-Submission Verification

1. **Pre-submission QA checks:**
   - Main blind PDF strictly **≤40 double-spaced pages** (target 35–38 pp.), zero LaTeX/BibTeX warnings, zero overfull hboxes in tables, zero `??` cross-references across Main ↔ Online Supplement.
   - All numerical claims in Abstract, Intro, MC, Empirics, Conclusion, and Online Supplement cross-checked against `r_project/results/` CSVs.
   - Confirm whether calendar-lag vs. index-lag Bartlett SHAC produces identical substantive conclusions (WP1.2 check).
2. **JBES Cover Letter (`jbes_cover_letter.md` / `.pdf`):**
   - Explicitly position the paper as general financial/panel econometrics for irregular non-traded cash-flow data.
   - Highlight the four novel methodological contributions (especially Proposition 2's duration-matching & cross-vintage impossibility results and the primitive increasing-domain fixed-$b$ framework).
   - Respectfully note alignment with the financial econometrics editorial track (e.g., Dacheng Xiu / Andrew Patton, where permitted by the submission portal).
   - Include a formal **ASA Data-Exemption Request** for proprietary Preqin fund-level cash flows, documenting our full public simulation & synthetic-data replication archive.

---

## 7. Recommended Execution Order & Effort Estimate

| Phase | Work Package | Est. Effort | Output / Gate |
| :---: | :--- | :---: | :--- |
| **Phase 1** | **WP1.3** ($V=400$ studentization & calendar-vs-index SHAC diagnostic in R) | 0.5–1 day | Resolved root cause + exact diagnostic numbers (no untested speculation!) |
| **Phase 2** | **WP1.1 + WP1.2** (Primitive Assumptions A1–A5, Lemma D.1, 1D Fixed-$b$ Corollary & proof in Appendix D) | 2–3 days | Referee-proof asymptotic & fixed-$b$ theory in `sdf_paper_online_appendices.tex` |
| **Phase 3** | **WP2 + WP3** (Reframe Intro/Abstract/Lit + LOO/JIVE citations + Reorder Secs 1–6 + Prune to 6–7 floats & ≤38 double-spaced pp.) | 3–4 days | Clean JBES main manuscript + synchronized Online Supplement |
| **Phase 4** | **WP4** (JBES blind template, anonymize all files, finalize non-SDF Title, assemble Replication Archive + README) | 1–2 days | `jbes_main_blind.pdf`, `jbes_supplement_blind.pdf`, `jbes_titlepage.pdf`, `jbes_replication.zip` |
| **Phase 5** | **WP5** (End-to-end numerical audit, visual PDF inspection, Cover Letter + Preqin exemption request) | 0.5–1 day | Submission-ready JBES package |

**Total estimated timeline:** **~7–11 focused working days** (with every single upgrade also transferring 1-to-1 to JFEC / JEF / *Econometric Reviews* if ever needed as a fallback).

---

## 8. Pre-Submission Audit, Fixes, and Decisions Log (from 2026-11-27)

Repository `~/github/Fundwise_SDF`, branch `revision-2026`. Nothing in this section has been pushed unless stated otherwise. Rule followed throughout: every number in the paper is recomputed from a file or produced by a script before it is stated.

### 8.1 Commits

| Commit | Content |
| :--- | :--- |
| `5db48f6` | JBES refinement 1: the package as audited. |
| `8585b0d` | Critical audit fixes (8.2). |
| `718cf64` | Clean empirical rerun from `8585b0d`; all exhibits regenerated (8.3). |
| `cc1d53c` | Moderate audit items (8.4). |
| `97c3c7f` | KV05 cubic-approximation footnote (8.5). |
| `5510a05` | Horizon-invariance test, two-moment test, VC beta band, and this Section 8 log (8.7–8.8). |
| `733dbc1` | The test outputs regenerated on the clean commit `5510a05`. Only the provenance columns changed; the table body and the PDFs are unchanged. |

Nothing above is pushed; `revision-2026` is 7 commits ahead of `origin` after the commit that adds this table row.

### 8.2 Referee-style audit of `5db48f6` and the critical fixes (`8585b0d`)

The audit checked Appendix D (A1–A5, Lemma D.1, Prop 3, Cor 1 / Prop D.2), the code, every headline number against the CSVs, the fixed-*b* critical values, and blind/packaging hygiene.

1. **Text claims contradicted by the CSVs**, corrected:
   - *Calendar- vs. index-lag Bartlett weighting.* "No test decision changes" was wrong: 5 of 80 p-values cross a conventional threshold. No H = 0 decision at 5% changes.
   - *Score test.* It was called calibrated, but it rejects 7.8% at V = 800, H = 0.
   - *SD scaling.* "Exact V^{-1/2}" became "close to the V^{-1/2} rate".
   - *Introduction.* "AR tests restore nominal size" was removed: the AR replay rejects 9–22%.
2. **Unsourced numbers.**
   - The kurtosis of β̂−1 at V = 400 (2.97–3.47) is now written by `simulation/mc_v400_error_moments.R` to `results/mc_step4/v400_error_moments.csv`.
   - The Table AS1 note is computed in `mc_paper.R` instead of being hard-coded.
3. **A1/A4 logical gap.**
   - A fixed Θ with an almost-sure floor on gross returns is impossible with unbounded (e.g. Gaussian) factors, so A1 now assumes bounded factor support.
   - The A4 remark "automatically satisfied with lognormal shocks" was wrong; it is replaced by explicit envelope bounds.
4. **Replication archive.**
   - The archive could not rebuild what the README claimed, so the missing configs, scripts and small precomputed inputs were added.
   - The README runtimes were corrected: the Monte Carlo study takes about 9.6 h and the empirical estimation about 5 h.
   - `jbes/build_jbes_replication.sh` refreshes the archive from `r_project/`, leak-scans it and zips it.
5. **PDFs.** Recompiled with 0 `??`. The main text is 39 pp and the cover letter says 39.

Verified clean in the audit:
- The Hölder (2,4,4) step needs exactly 16+8δ moments; 8+4δ would not suffice.
- The proofs of Lemma D.1(a), the KV functional and summation by parts, the sparse-gap embedding, the H₀* formula and bounds, and the M_m eigenvalues all check out.
- The fixed-*b* critical values 2.790 / 2.867 / 2.722 were reproduced independently.
- All headline and Section 4 numbers match the files.

### 8.3 Clean empirical rerun and provenance (`718cf64`)

- **Why.** Every empirical number had come from caches written with uncommitted changes (`git.dirty` TRUE).
- **What was rerun.** The 26 exhibit scenarios were rerun from `8585b0d` in a clean scratch worktree. All 806 cache files (37,128 rows) are bit-identical to the earlier dirty caches. The superseded caches were moved to `r_project/results/data_out_2026-emp-superseded-731ea2f/` (gitignored).
- **Pipeline fact (non-obvious).** The exhibit scripts read `0_*_summary.csv`, written by `results/analyze_result.R`, and nothing ran that step for the current result folders. The new `results/analyze_empirical_results.R` must run between `run_empirical_study.R` and `export_g1_g2_tables.R`; archive README §4.3 lists it.
- **L1_Ridge tables.** They are not exhibits of the paper, and their scenarios were not re-estimated after the Nov 22 data change. `export_g1_g2_tables.R` now skips them with a message, and the tracked L1 tables are stale (Nov 18).
- **Generation stamps.** Each exhibit's generation stamp reads "(uncommitted changes)" because earlier scripts in the chain had already rewritten tracked outputs. Only the estimation-provenance line matters, and it reads `8585b0d`, uncommitted FALSE.
- **Known check quirk.** Check 3 of `check_main_u.R` reports "same sequence FALSE" for the E1 tables because `U.p1` is recorded first. The tables themselves are correct.

### 8.4 Moderate audit items (`cc1d53c`)

1. **Per-vintage score centering.**
   - *Change.* Prop 3, Cor 1, Prop D.1 and D.2 now require E[s_v(θ†)] = 0 *for every vintage*. New Remark `rem:per_vintage_centering`.
   - *Why.* The definition of θ† only zeroes the cross-vintage average. Centering holds for U under T7 for any pool sizes, because each vintage's expected criterion is stationary at θ₀. It fails for L2 when n_v varies, because the T3–T6 variance tilt moves each vintage's own minimizer. Hence significance markers go on U only.
2. **A5 homogeneity versus trending n_v.**
   - *Change.* "Exact KV reference distribution used throughout" is dropped.
   - *Why.* n_v rises from 2 to 41 (Buyout) and 56 (VC U), so the KV values are an approximation. The justification is size by design: the core Monte Carlo designs use the observed calendars and fund counts, and U_{v,0} rejects 4.6% / 6.4% there.
3. **Table AS6 rejection rates.**
   - *Change.* The table now has Rej.% rows.
   - *Why.* Horizon choice fixes bias but not inference. At H = 12.5–15 under persistent cohort or calendar shocks, tests reject a true α in 18.4–20.8% of draws at V = 25 (19.0–23.8% around the pseudo-true value) and in up to 17.2% at V = 100. The main text states this too.
4. **The −5.12%/yr at H = 30, V = 100** is a finite-sample effect, not a contradiction of T7. The target shift is 0.00 and there are no bound hits, but the SD jumps from 0.30 to 7.74. The Table AS6 note says so, and the caption was shortened into a note so the table fits its page.
5. **β at H = 15 under shared shocks.** It is biased by −0.074 (Buyout 2c) and −0.191 (VC), entirely because the target itself shifts (−0.083 / −0.199); tests of β = 1 reject 9.6%. This is stated in the main text and the Table AS3 note.
6. **Symmetric deal-shock ablation.** It makes 22% of payoffs negative (37% for deals held ≥ 11 years; the light version 2.3%), which is outside Assumption 1. The shares are computed exactly in `mc_v400_table.R` and the ablation is labelled as such.
7. **"Lowest RMSE"** now reads "among estimators with near-nominal tests". The text discloses that q_v = n_v has RMSE 0.056 / 0.086, and that a jackknife still rejects 10.3%.
8. **Exposition.**
   - *Prop 2(i).* The closed-form H₀* is now marked as continuous-time. `simulation/optimal_H_theory/discrete_vs_continuous_H0.R` shows the monthly-grid roots are at most 0.42 months shorter.
   - *Prop 2(iv).* "Driving minimizers to ∂Θ" is replaced by the Monte Carlo evidence (42–100%); indefiniteness alone does not prove it.
   - *Bandwidth ratio.* b = S/V is the conservative choice relative to S/V_cal (critical values 2.753 at 7/26 and 2.665 at 7/29). `helper/fixedb_q1_omega.rds` gained one cache entry (M = 241); the existing entries are unchanged.
9. **Page counts.** The main text is still 39 pp; the supplement went from 59 to 60 pp.

### 8.5 KV05 footnote (`97c3c7f`)

- **Source.** Kiefer and Vogelsang (2005), Table 1, p. 1146, Bartlett 97.5% quantile: cv(b) = 1.9600 + 2.9694b + 0.4160b² − 0.5324b³. The copy is in `jfinecon_referees/papers/`.
- **Size of the gap.** The cubic gives 2.812 / 2.887 / 2.749, i.e. 0.020–0.027 above ours.
- **Nature of the gap.** It is a level difference between the two sets of simulated values, not a property of the cubic form. A cubic refitted to our grid matches it within 0.023, whereas KV's cubic lies above our whole grid (mean 0.027). The cause is not traced.
- **Effect on results.** No main-sample Buyout/VC U test changes its 5% decision. Five robustness-sample tests would, and none of them is printed as a test.
- **Script.** `results/fixedb_kv05_comparison.R`.

### 8.6 Decisions (with reasons)

| Date | Decision | Why |
| :--- | :--- | :--- |
| 2026-11-19 | q_v = n_v U columns: point estimates only, no SE or stars. The jackknife is kept in reserve for a referee. | q_v = n_v over-rejects about 20% under Bartlett; the jackknife brings this down to 10.3%, still twice nominal. |
| 2026-11-27 | **No PitchBook data before submission.** Keep it as an independent replication sample for an R&R. | It does not resolve the VC identification question, and the complete-history rule adds only about 5 vintages (V ≈ 30, N_eff stays single-digit). The cost is weeks: license, a new prep pipeline, recalibrated MC designs (size by design rests on matching the sample's calendar and counts), regenerated exhibits, and a new audit. First check PitchBook's coverage of complete cash-flow histories per vintage for North American BO/VC. |
| 2026-11-27 | **"Which VC beta": do options 1 and 3 now, option 2 as a short appendix result** (8.7). | Option 1 turns the referee's question into a test the paper runs itself. |
| 2026-11-27 | **No AssetMetrix affiliation disclosure** in the disclosure statement. | The author's decision. |
| 2026-11-27 | **Opening sentence of the introduction narrowed** to "private capital funds (private equity including buyout and venture capital, real estate, infrastructure, natural resources) and project finance". Private credit, commercial real estate in general, and CLOs were dropped from the list. | The old claim ("researchers rarely observe periodic traded returns" across seven areas) was false as written for areas where traded prices exist: commercial real estate (REITs; transaction and appraisal indices), CLOs (traded debt tranches; J.P. Morgan CLOIE, a daily total-return index on observable prices since 2014), and parts of private credit (listed BDCs, syndicated-loan indices). The new sentence says what holds for funds and project finance: sparse secondary trades, appraisal or book values, cash flows as the reliable record. The author set the scope. The abstract's list ("private equity, venture capital, private credit, real estate, and infrastructure") was left unchanged: it says these assets generate pooled cash-flow streams, which is true for private credit funds. |
| 2026-11-28 | **Reposition fixed-$b$ inference out of Keywords and soften in Abstract / Cover Letter** (keep intact in §3.3, §4, §5, and Appendix D). | Fixed-$b$ Bartlett inference (KV05, LLSW18) is an off-the-shelf finite-$V$ HAC correction rather than our methodological invention. Keeping it in the body is essential for $V=23\text{--}27$ overlapping vintages, but listing `Fixed-$b$ inference` in Keywords risks misrouting to HAR-asymptotics theorists instead of financial econometrics / asset-pricing / panel-GMM referees, and saying "develop one-dimensional fixed-$b$ ... inference" in the Abstract both adds opaque jargon and slightly oversells. Drop from Keywords (replace with a core-contribution keyword) and rephrase in Abstract/Cover Letter as small-sample overlap- and identification-robust inference. |
| 2026-11-28 | **Accept all Round-3 referee catches (GPT 6.2 Astra) after independent verification**; plan Round-3 polish pass (zero re-estimation). | All three counterexamples (Table 1 L2 conditioning under (P) vs. (X); Prop. 2(ii)(b) $\gamma_e \ge 0$ for common-entry calendar shocks; finite-vintage $O(1/V)$ bias centering under heterogeneous $c_v = c + d/\sqrt{v}$) and both bookkeeping items (normalization $Q_H/\bar{q}$ + random eligibility indicator; $h_v\sigma_{\nu,v}^2$ diagonal contribution + "can fail" at L484) verified exact from first principles. Full record in `archive/jbes_referee_response_decisions.md` §5 and checklist in `jbes_referee_revision_plan.md` §6. |
| open | Push `revision-2026`. | The author's decision. |

### 8.7 "Which VC beta" weakness: plan

- **The problem.** VC U β is 1.493 / 1.742 (EW/FW) at H = 0, 1.382 / 1.583 at H = 12.5 and 1.250 / 1.423 at H = 15 (`upair_VC_US_data.csv`). The headline at H = 0 is valid only without shared vintage shocks (T7). Yet the paper uses shared shocks to explain the VC alpha ceiling and calibrates its VC Monte Carlo with σ_W = 0.50. The paper has no test of which case holds.
- **Option 1: horizon-invariance test.**
  - *Logic.* Under (X) and c_v = 0, U targets the same β at every H, so a significant β̂_0 − β̂_H rejects the assumptions behind the headline.
  - *Method.* The SE comes from one Bartlett sandwich over the stacked per-vintage scores of both horizons, with the fixed-*b* critical value. Buyout serves as a check.
  - *Validation.* Size and power are checked by Monte Carlo on the core designs (c_v = 0) and the shared-shock designs.
- **Option 2: joint two-moment test at H = 0** (short appendix result). At β₀ the mean pricing moment is zero either way, but the pair-product moment is zero only if c_v = 0. Inverting a joint fixed-*b* q = 2 test over β, an empty set means shared shocks. Low power is expected.
- **Option 3: framing.** Present VC β as a band from H₀* to H = 0, with the bias direction known: under c_v > 0, U_{v,0} is biased upward. Conclusions to check at both ends: β > 0 rejected, β = 1 not rejected, point estimates below the squared-loss 2.73 of DLP12. Flag the 18–21% over-rejection at H = 12.5–15 under persistent shocks.

### 8.8 "Which VC beta": results (2026-11-27, after `97c3c7f`)

**Code.** All of it is in `r_project/` and copied to the replication archive; README §4.1 step 6 and §4.3 step 5 list it.
- `helper/horizon_invariance.R` implements both tests.
  - `u.horizon.contrast()`: the stacked M-estimator (β̂_0, β̂_H) of the separable criterion Q_0 + Q_H on the same vintage units. The SE is the q = 1 contrast of the stacked Bartlett sandwich: influence differences d_v = g_{v,0}/h_0 − g_{v,H}/h_H, se² = d′Kd/V², calendar-distance kernel, S = ⌈1.3√V⌉, b = S/V, fixed-*b* p-value.
  - `u.two.moment.test()`: per vintage (n_v ≥ 2) at H = 0, the moments are g1 = Σe/Σs (vintage mean error) and g2 = U_v. A centered fixed-*b* Wald statistic with q = 2 is inverted over a 401-point β grid on the admissible domain.
- `results/horizon_invariance_test.R` runs both tests on the Preqin panels (~30 s) and writes `results/figures_bovc_us_g1/horizon_invariance_test.csv` and `two_moment_test.csv`. Checks:
  - its influence code reproduces the production SEs to 1e-8;
  - all 40 β̂ match `upair_*_US_data.csv` within 5e-8.
- `simulation/mc_horizon_invariance.R` runs the Monte Carlo (8.6 min on 10 workers) on nine Step 4 designs with common random numbers.
  - Replication r is `mc.draw(design, r)`, and replication 1 reproduces the production estimates exactly. The archive reproduces it from a fresh unzip.
  - Output: `results/mc_step4/horizon_invariance_mc.csv`. The replication rows go to `simulation/data_out_mc/horizon_invariance/` (gitignored).
- `results/horizon_invariance_tables.R` writes Supplement Table D.1 (`table_OA_horizon_invariance.tex`, label `tab:horizon_invariance`) and the printed numbers to `horizon_invariance_exhibit_data.csv`.

**Preqin results (Supplement Table D.1, Panels A–B).**

| Panel | H₁ | β̂_0 | β̂_H₁ | Diff | SE | Corr | t | p |
| :--- | :---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| VC EW | 15 | 1.493 | 1.250 | 0.243 | 0.142 | 0.91 | 1.71 | 0.200 |
| VC FW | 15 | 1.742 | 1.423 | 0.319 | 0.147 | 0.93 | 2.18 | 0.115 |
| VC EW | 12.5 | 1.493 | 1.382 | 0.111 | 0.121 | 0.95 | 0.92 | 0.468 |
| VC FW | 12.5 | 1.742 | 1.583 | 0.159 | 0.130 | 0.95 | 1.22 | 0.343 |
| Buyout (all) | 5–20 | | | | | | | ≥ 0.513 |

- **Other horizons.** Over the full profile (H₁ = 5, 10, 12.5, 15, 20), the largest |t| is 2.58 (VC FW, H₁ = 5, p = 0.071), where the VC profile peaks. This is disclosed in the supplement text.
- **Two-moment test.** Its largest p-value is 0.52–0.86, and its sets cover 27–65% of the domain: no rejection, not informative.

**Monte Carlo (Supplement Table D.1, Panel C; 5% rejection rates).**
- **Size (c_v = 0).** EW cores: 3.8% (BO) and 2.0% (VC) at H₁ = 15, 4.8% and 3.4% at 12.5. FW with fund sizes: 3.2% and 1.4% at 15, 7.0% and 3.6% at 12.5. SE/SD is 0.64–0.72, so the SEs are too small, but the fixed-*b* critical value compensates.
- **Power at H₁ = 15.**

  | Design | Power at 15 | Power at 12.5 |
  | :--- | ---: | ---: |
  | m4_shared (c) | 8.4% | 6.8% |
  | m4p_shared2 (2c, persistent) | 15.6% | 12.4% |
  | m4_shared2 (2c) | 20.8% | 15.2% |
  | u_vc_shared (VC scale, σ_W = 0.50) | **33.2%** (49.6% at 10%) | 18.4% |

  Power is lower at 12.5 because the target overshoots beyond H₀*.
- **Biases in u_vc_shared.** U_{v,0} +0.249, U_{v,12.5} −0.008, U_{v,15} −0.135. So **the zero-bias horizon for β is about 12.5 yr, not 15**; the same holds in m4_shared2 (+0.143 / −0.011 / −0.056). Mean diff at H₁ = 15: −0.007 under c_v = 0 versus +0.384 at VC scale.
- **Path-correlated design A′ ((X) fails).** Both targets move about equally (+0.32 at H = 0, +0.35 at H = 15), so the test has no power there. Its 16–21% rejections come from SE/SD = 0.30, and this is disclosed.
- **Two-moment test.** 0% rejections in every design, so it has no power at V ≈ 25. It is kept in the supplement only as a documented negative result.

**Decisions taken from these results (2026-11-27).**
1. **Primary test horizon H₁ = 15 yr**, fixed before the power results; 12.5 is secondary. Why: the test is valid at any H₁ under the null, and the overshoot beyond H₀* makes 15 the more powerful choice.
2. **The shared-shock end of the VC band is H = 12.5, not 15**: VC β is 1.382 / 1.583 at 12.5 and 1.493 / 1.742 at inception. Why: under the VC-scale shared shock, U_{v,12.5} is unbiased (−0.008) while U_{v,15} overshoots (−0.135). This replaces the "1.25–1.42" I quoted earlier in the conversation.
3. **The band is reported in the abstract, introduction, empirical section and conclusion.**
   - Abstract: "Venture Capital betas of 1.49–1.74 at inception, or 1.38–1.58 if vintages share common shocks"; it now has 198 words (limit 200), and the cover letter says 198.
   - Both ends reject β = 0 at 5% (H = 12.5: p₀ = 0.017 / 0.014) and neither rejects β = 1 (p₁ = 0.418 / 0.269 at 12.5; ≥ 0.116 at 0).
4. **Non-rejection is reported as weak evidence.** The VC gaps (0.24 / 0.32) lie between the c_v = 0 expectation (≈ −0.01) and the VC-scale expectation (0.38), and the power against VC-scale shocks is only 33%.
5. **Placement.**
   - Main text, Section 3.3: a short paragraph "Horizon-invariance test", with no displayed equation to save space.
   - Supplement Appendix D: new subsection "Specification tests of the U_{v,H} identifying conditions", with Remark D.2 (validity and consistency: stacked separable criterion, block-diagonal Hessian, per-vintage centering under T7, Cor 1(i) applied to a = (1, −1)), Remark D.3 (two-moment test) and Table D.1.
   - This renumbers the later D floats: the V = 400 table is now D.2, the mc_paper D tables are D.3–D.9, and the calendar-vs-index remark is D.4. Archive README updated. The README's "Remark D.1" had been stale since `cc1d53c`.
6. **Page counts.** The main text is 40 pp (plan limit ≤ 40; it was 39), the supplement 62 pp; the cover letter says 40.

### 8.9 Still open (referee-level, not errors)

1. Inference with V ≈ 25 and N_eff of 4.5–8.9; the AR sets reach the domain ceiling.
   - The horizon-invariance test has only modest power (≤ 33%) against VC-scale shared shocks.
   - To raise power, a referee could ask for more vintages. PitchBook data could later extend V to about 30.
2. Mixed messages and scope: "the horizon displays the bias" versus "a unique zero-bias horizon"; the paper reads like two papers; U's novelty over leave-one-out quadratic forms rests on the compounding argument.

### 8.10 Methodological referee-report revision pass (2026-11-27, after `bd6be4b`)

Full item-by-item checklist and verification log: `2025 sdf paper/jbes/jbes_referee_revision_plan.md` (addressing `jbes_methodological_referee_report.pdf`). All P0–P3 items completed with **zero re-estimation** (verified in `prep_preqin.R`, `pool_vintage.R`, and `objective.R` that FW fund-size weighting has no double-sizing bug; the issue was purely a notation contradiction in Appendix A Table A.1).

- **P0 mathematical & formal fixes:** Deleted the false Cauchy–Schwarz $H=0$ efficiency theorem in Prop. 2(iii) / Appendix C.11(iii) and added the exact sample-score identity $s_{v,H} = A_{v,H}^2 s_{v,0} + 2A_{v,H}(\nabla A_{v,H})U_{v,0}$ plus the two-fund counterexample and actual-MSE $2B_0'b_1$ term; restricted Prop. 2(i) to symmetric pair-averaged $\Sigma_v$, $T_{\mathrm{eff}} > 0$ conditions, attainable-range caps, score-root vs. global-min distinction, and moved Brownian volatility shrinkage to a clearly labeled stylized multiplicative-shift special case; replaced Prop. 2(ii)'s unconditional calendar-overlap step with the exact positive exit-overlap covariance condition and fixed level-vs-slope/sign claims; updated Lemma C.1(a)–(b) (single-first-call $U_{v,0}\to 0$ and $2\tau^\ast > s_v^{\min}$ coercivity), restricted A2's $\underline{d} > 0$ bound to $n_v \ge 2$ vintages, removed the erroneous Merton (1980) precision lower bound, and tightened §2.1–2.2 adapted measurability, $(a,I_a)$ conditioning in (X), selected-pool (P), and T7 sign-aligned proportionality + Hessian spanning.
- **P1–P2 attribution, framing & hierarchy:** Fixed Appendix A Table A.1 FW units and qualified DLP (2012) dollar-sum ($S_v^2$) vs. inverse-FoF-size normalization plus $U_{v,H}$ pair weights $\omega_i(1-\omega_i)$; cited Tausch & Pietz (2024, eq. 4) in third person where eq. (8) defines horizon averaging; credited DLP (2012, fn. 8 & p. 522) in the main text; narrowed Prop. 2(iv) to *uninstrumented* cross-vintage pairs and added matched-assumptions comparator discussion citing Domínguez & Lobato (2004) and Muandet, Jitkrittum & Kübler (2020); qualified T8 per-vintage centering at $H=12.5$~yr, $S/V \to 0$ fixed-$b$ finite-sample approximation, and HI specification diagnostic; replaced VC "band" language with side-by-side sensitivity estimates, labeled all single-factor results "$\alpha$ fixed at zero", highlighted matched 23-vintage VC L2–$U$ contrasts ($+0.435$/$+0.551$), separated AR sets from Wald tests, qualified Fix-$A$ counterfactuals and held-out L2 CV; added the 2-fund lead-in and new Table 1 (`tab:criteria_comparison`, moving `tab:regime_summary` to Appendix C), added question-based MC headings, and tightened the abstract to 184 words.
- **P3 exhibits, PDFs & archive:** Updated auto-generated notes in `mc_paper.R`, `main_exhibits_u.R`, `horizon_invariance_tables.R`, `table_published_estimates.R`, and `step5_diagnostics.R` (zero numerical changes); recompiled `jbes_main_blind.pdf` (**40 pages**, 3 figures, 5 tables, 0 warnings/undefs) and `jbes_supplement_blind.pdf` (**67 pages**, 0 warnings/undefs); synced `jbes_cover_letter.tex/.md/.pdf` and rebuilt `jbes_replication.zip` (clean leak scan).

### 8.11 Second-round focused referee revision pass (2026-11-27, after `a4d444e`)

Full item-by-item log appended in `2025 sdf paper/jbes/jbes_referee_revision_plan.md` §5, and follow/pushback decision record archived alongside both referee reports in `2025 sdf paper/jbes/archive/jbes_referee_response_decisions.md` (addressing `2025 sdf paper/jbes/archive/jbes_second_round_referee_report.pdf`, *"A focused revision before submission"* — zero re-estimation required).

- **Table 1, Intro, Conclusion, Lemma C.2, Prop. 2(iv), Appendix B:** Tightened all three population targets in Table 1 (L2 interior displacement via nonzero noise-contribution gradient; $U_{v,H}$ shared-shock displacement via nonzero vintage-average gradient of $q_v A_{v,H}^2 c_v$ + T7 global separation; inception mean roots vs. nominal pointwise AR coverage); revised Intro/Lemma C.2/Table B.1 to require diffuse weights ($h_v = \sum_i \omega_{v,i}^2 \to 0$, sufficient route to L2 consistency absent shared shocks); updated Appendix B to cite Theorem 1(T7); qualified Prop. 2(iv) to "can have multiple roots" since $A_{v,H}(\beta)$ at $H>0$ can shift or remove $\beta^\ast$.
- **Worked example & T3 proof:** Stated both per-period $R^\dagger = R_0(1+\sigma_N^2/n)^{1/T}$ and cumulative holding-period conventions in §2.2; conditioned on $\mathcal{F}^{\mathrm{pub}}$ ($\mathcal{F}_{\max(T,K-1)}^{\mathrm{pub}}$); aligned Supplement §C.4 T3 proof with eq. (9)'s $n=2$ coefficient $\sigma_N^2$.
- **Horizon 12.5 yr, Figure 1B, Prop. 2(ii)–(iii):** Relabeled Figure 1B vertical axis (`figure_optimal_H_theory.tex` L126) and Supplement Tables C.2–C.3 to *"Zero-score horizon"*; distinguished the stylized $12.65$-yr duration calculation (motivation) from the VC shared-shock MC $-0.008$ bias (calibration evidence); refined Prop. 2(ii) (independent $e_u$ / equal outlays vs. $\Gamma_e$, nonnegative common-entry covariance, linear-component contribution to $\Lambda_H^{\mathrm{cross}}$); stated eq. (17) as a **local** MSE optimum near $H_0^\ast$ with uniform differentiability in Prop. 2(iii); restored $q_v$ in $2q_v A_{v,H}(\nabla_\theta A_{v,H})U_{v,0}$ (Prop. 2(iii) & eq. (C.3)); fixed random-weight expectation $\mathbb{E}[\widehat{U}_H\mid\mathcal{W}]$ in Appendix C.11(ii)(a).
- **Attribution & remaining alignments:** Added $(1-h_v)^{-1}$ reweighting remark after eq. (11) and phrased $+0.435$/$+0.551$ in §§5.2–5.3 as the *"matched 23-vintage difference between pooled L2 and the normalized pair estimator $U_{v,0}$"*; defined pool mean $\mu_v(\theta) := \mathbb{E}[PV_v(\theta)\mid\mathcal{G}_v]$ in T4/§C.5, marked $\mathbb{E}[N\mid\mathcal{G}_v,a,I_a]=1$ as sufficient for (X), flagged future-reporting eligibility as an additional selected-cohort restriction in (P); qualified §5.3/Conclusion coercivity via Lemma C.1(b) and distinguished fixed eq. (18) instrument $z_{v,H}$ from parameter-dependent $A_{v,H}(\alpha,\beta)$ in Remark C.1(iii); added finite-dependence/orthogonality condition after Lemma C.3; restricted DLP (2012) inverse-FoF-size equivalence strictly to FW in §2.2, §3.4, Table A.1, and `main_exhibits_u.R`.
- **QA & archive:** Recompiled Figure 1, re-ran `main_exhibits_u.R` (zero diff in `exhibits_data.csv`), rebuilt `jbes_main_blind.pdf` (**40 pages**, 3 figures, 5 tables, 0 warnings, 0 >10pt overfull hboxes, 184-word abstract, clean blind scan) and `jbes_supplement_blind.pdf` (**67 pages**, 0 >10pt overfull hboxes, clean blind scan), and refreshed `jbes_replication.zip`.

### 8.12 Round-3 referee revision pass (2026-11-28, after `971680e`)

- **Source:** Third-round AI referee report (`2025 sdf paper/jbes/archive/jbes_third_round_referee_report.pdf`, GPT 6.2 Astra; closes the main objection to $U_{v,H}$, retains the VC application as written, and accepts the unified T1–T8 structure) plus author editorial TODO on fixed-$b$ positioning.
- **Verification & full record:** All three substantive mathematical findings (with explicit counterexamples) and both bookkeeping catches were independently verified from first principles and implemented cleanly with zero re-estimation / zero R changes (full algebra & rationale in `2025 sdf paper/jbes/archive/jbes_referee_response_decisions.md` §5; completed checklist in `2025 sdf paper/jbes/jbes_referee_revision_plan.md` §6).
- **Implemented edits (4 workstreams):**
  1. **R3.1 (Fixed-$b$ positioning):** Replaced `Fixed-$b$ inference` with `Horizon selection` in Keywords; revised Abstract to *"and establish valid small-sample (fixed-$b$) and identification-robust inference"*; aligned Cover Letter editorial-track phrasing; kept fixed-$b$ intact in §3.3, §4, §5, and Appendix D.
  2. **R3.2 (Table 1 L2 row + Abstract/Conclusion wording):** Separated T3's $\mathcal{F}_{t_v^0}$-conditional variance under (P) from T4–T6's $\mathcal{G}_v$-conditional $A_{v,H}^2\sigma_{\mathrm{eff},v}^2$ under (X), both vintage-averaged (with matching note in Supplement §C.4 and T3); changed "tilts/shifts" $\to$ "can shift" in Abstract, Conclusion, and Cover Letter.
  3. **R3.3 (Prop. 2(ii)(b)):** Restored *"with nonnegative calendar autocovariances"* ($\gamma_e \ge 0$, matching Supplement §C.11(ii)(b) L433) into the common-entry parenthesis of Main Prop. 2(ii)(b) L408.
  4. **R3.4 (Finite-$V$ bias centering + two bookkeeping fixes):** Added Bao & Ullah (2007, *JoE* 140(2), 650–669, `bao2007timeseries`); centered §C.11(ii)(c)'s stochastic $O(1/V)$ expansion on $\theta_{V,H}^\dagger$ plus $(\theta_{V,H}^\dagger - \theta_H^\dagger)$ with explicit $V^{-1}\sum\mathbb{E}[s_{v,H}(\theta_H^\dagger)] = \gamma_1(H)/V + o(V^{-1})$ (or $o(V^{-1})$) rate condition + higher-order regularity across Prop. 2(ii)(c), Prop. 2(iii), and §C.11; noted competing-root ties in $\mathcal{V}$ unresolved at $O(V^{-1})$; harmonized $\widehat{U}_H \xrightarrow{p} Q_H/\mu_q$ (using $\mu_q$ to avoid clashing with A2's $\bar{q}$) and moved eligibility indicator $E_v = \mathbf{1}_{\{n_v\ge 2\}}$ inside the expectation in §C.11(ii)(a); used *"idiosyncratic diagonal contribution"* ($h_v\sigma_{\nu,v}^2$, $1/n_v$ under EW) at Main L109/L319/L407 & Table C.1 and changed "fails" $\to$ "can fail" at L484 (softening Remark D.1 L745).
- **QA & archive:** Rebuilt `jbes_main_blind.pdf` (**40 pages**, 32 lines on p. 40, 3 figures, 5 tables, **186-word abstract**, 0 warnings, 0 >10pt overfull hboxes, 0 `??`, clean blind scan), `jbes_supplement_blind.pdf` (**68 pages**, 0 >10pt overfull hboxes, 0 `??`, clean blind scan), and `jbes_cover_letter.pdf` (**2 pages**); archived `jbes_third_round_referee_report.pdf` in `2025 sdf paper/jbes/archive/`.
