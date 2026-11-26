# JBES Fit Analysis: Benchmarking Against 10 Recent JBES Papers (2025–2026)

*Date:* November 26, 2026  
*Target Journal:* *Journal of Business & Economic Statistics* (JBES)  
*Manuscript Audited:* `2025 sdf paper/tex/sdf_paper.tex` (`sdf_paper.pdf`), `sdf_paper_online_appendices.tex`, `sdf_paper_appendix_sdf.tex` (`sdf_paper_appendix.pdf`)  
*Benchmark Corpus:* 10 recent JBES papers (2025–2026 / accepted working versions) stored in this folder (`2025 sdf paper/jbes/`).

---

## Bottom-Line Verdict: **Yes — Structurally ~75–80% There (Much Closer Than RAPS/RFS)**

Benchmarking `sdf_paper.tex` and `sdf_paper_online_appendices.tex` against all 10 recent JBES papers confirms that **JBES is a dramatically better structural and cultural fit than RFS or RAPS.**

Unlike RAPS/RFS—where recent papers keep main texts almost equation-free and judge primarily on novel economic findings and fresh data—all 10 JBES benchmarks share our paper's exact DNA: formal identification/asymptotic theory in the main text, proofs in an Online Supplement, extensive Monte Carlo validation, and a focused empirical *illustration*. Moreover, **Dacheng Xiu** (Chicago Booth) is currently one of JBES's three Joint Editors (alongside Yingying Fan and Michal Kolesár), with **Andrew Patton** on the AE board, and Massacci (2025) and Liu (2025) prove JBES actively publishes empirical asset-pricing and fund-performance econometrics citing the exact JFEC literature we cite (e.g., Burnside 2016).

That said, **we are not submission-ready for JBES as-is.** Five concrete gaps—two substantive econometric ones and three structural/length ones—stand between the current draft and a realistic shot at an R&R.

---

## 1. Head-to-Head Scorecard: Our Manuscript vs. the 10 JBES Benchmarks

| Dimension | 10 Uploaded JBES Papers (2025–2026) | Our Current Draft (`sdf_paper.tex`) | Closeness |
| :--- | :--- | :--- | :--- |
| **Main-text length** | 11–16 printed 2-col pp. (working papers: **34–36 pp.**, Massacci longest at ~40–44 double-spaced pp.); JBES rule: **15–30 pp., ideally ≤40 double-spaced** | **37 pp. text + 4 pp. refs single-spaced** (~16.3k tokens ≈ **58–62 double-spaced pp.**) | **Too long (~1.5× ceiling)** |
| **Online Supplement** | Standard & unlimited (up to **85 pp.** in Mehrabani & Parsaeian); holds **all proofs** (9/10 papers) + extra MC/empirics | Already split cleanly in `sdf_paper_appendix.pdf` (**39 pp.**); only 2 short proofs remain inline | **Very close** |
| **Main-text floats** | **4–6 total** (2–4 tables + 1–3 figures; even Massacci's 16-pp. paper has 4 tables + 2 figures) | **14 total** (7 tables + 7 figures) | **Too float-heavy** |
| **Section order** | **8/8 empirical papers:** `Intro → Theory → Monte Carlo → Empirical App → Conclusion` | `Intro → Theory → Data/Imputation/NAVs + MC → Empirical App → Conclusion` | **Easy reorder** |
| **Math & theorem density** | 1–6 Theorems / Propositions, 4–8 numbered Assumptions, 1–2 Lemmas, Notation paragraph (6/10 papers) | 1 Theorem (T1–T8), 3 Propositions, 1 Lemma (main) + 3 Lemmas (App), 2 Assumptions | **Spot-on density**, but Ass. 2 is too high-level |
| **Role of empirical app** | **Illustrative** (2.5–6 pp., 1–2 tables/figs; 2/10 papers have **zero** real data) | 7.5 single-spaced pp. (+ 4 pp. of Preqin cleaning in Sec 3.1–3.3) | **Strong (vintage cutoff is NOT fatal here)** |
| **Framing** | General econometric method/test/ID framework first; dataset is an application | Currently reads like an audit/repair of one PE estimator (`DLP12`) | **Needs reframing** |
| **Title & Abstract** | Strictly descriptive methodological titles; ≤200-word plain abstract (no `\citet`, minimal math) | Title already close to JBES style (refine & drop "SDF"); Abstract is **303 words** with `\citet{DLP12}` & heavy math | **Easy polish** |

---

## 2. Where We Are Already Very Close

1. **Right architecture & math appetite:** Every single one of the 10 papers follows the exact `Motivation/Toy Model → Identification & Asymptotics → Monte Carlo → Illustrative Empirics → Online Proof Supplement` blueprint we already built. Even our 1-deal analytical example (`sdf_paper.tex`, lines 300–324) mirrors Wei (2025, Sec. 2 Toy Model) and Harding, Lamarche, & Muris (2025, Fig. 1 Motivating Example).
2. **Proposition 2 is genuinely JBES-grade theory:** Reading our proof of Proposition 2 in `sdf_paper_online_appendices.tex` (lines 279–472), this is our true theoretical centerpiece for JBES: the exact multi-deal duration-matching equation, exponential-family monotonicity and global bounds, FKG/Harris inequality proof of lognormal volatility shrinkage, Chebyshev active-life overlap shift for calendar shocks, second-order $O(1/V)$ bias expansion, finite-$V$ MSE expansion, and the matrix-spectrum proof of rank-1 collapse and finite-$V$ indefiniteness for cross-vintage leave-band-out pair products.
3. **Our Preqin sample cutoff (vintages ≤2010, April 2022 export) is NOT a fatal flaw at JBES:** Unlike RFS/RAPS/JFQA—where referees routinely reject over data vintage—JBES treats the empirical section as a proof-of-concept illustration (Hafner et al. 2026 re-examine Kilian's 2009 oil dataset; Liu 2025 uses $N=177$ CRSP mutual funds; Wei 2025 and Wang et al. 2022 use no real data at all).
4. **Companion software / reproducibility ethos:** Virolainen (2025) highlights his CRAN package `gmvarkit` right in the main text; ASA/JBES strongly rewards clean, reproducible R code and simulation suites—something we already have built in `r_project/`.

---

## 3. The 5 Real Gaps Between Our Current Draft and a JBES R&R

Ranked from **highest referee-rejection risk** to **easiest mechanical fix**:

### Gap 1 (Highest Risk): Asymptotic Proof Rigor & the Unresolved $V = 400$ Studentization Puzzle
A tough JBES econometrician (especially from Xiu's or Kolesár's pool) will go straight to **Section 2.4** (`sdf_paper.tex`, lines 585–617) and **Appendix D** (`sdf_paper_online_appendices.tex`, lines 532–570) and spot three specific vulnerabilities immediately:
- **High-level "black-box" proof:** Right now, `Assumption 2` and the entire proof of Proposition 3 occupy **one 10-line paragraph** (`sdf_paper_online_appendices.tex`, lines 536–566), assuming *directly* that the nonlinear criterion, score, and Hessian satisfy Jenish–Prucha (2009, 2012) NED/mixing and uniform $(2+\delta)$-moment bounds. Even JBES papers that cite existing CLTs (like Harding et al. 2025 citing Hansen–Lee 2019, or Massacci 2025 building on Bailey et al. 2021 and Qu–Perron 2007) state **primitive conditions on the underlying factor and cash-flow processes** (bounded fund lifetimes $L_{\max}$, mixing/moment bounds on $(F_u, \nu_i, \eta_v)$, gross-return floor $\underline{G}$) and *verify* that the compounded nonlinear products $A_{v,H}^2 PV_i PV_j$ and their derivatives inherit bounded moments, stochastic equicontinuity, and $m$-dependence / NED.
- **1D vintage calendar vs. 2D spatial random fields + fixed-$b$ FCLT mismatch:** Our vintages $y_v \in \mathbb{Z}$ live on a **1D calendar-time axis** (with a few skipped years: 1989 in BO, 1983–84 in VC), yet we invoke 2D spatial random-field theorems (Jenish–Prucha 2009/2012; Kim–Sun 2011) alongside 1D regular-lattice Kiefer–Vogelsang (2005) fixed-$b$ critical values—**without a formal statement or proof showing the fixed-$b$ limit distribution holds for our studentized $t$-statistic** under irregularly spaced/heteroscedastic vintages and estimated plug-in scores.
- **The smoking gun at `sdf_paper.tex` line 803:** Our own MC notes that **even at $V = 400$**, equal-vintage $U_{v,H}$ still rejects in **6.8%–10.8%** of draws with $N_{\mathrm{eff}} = 13\text{--}27$, attributing it tentatively to "heavy-tailed compounded weights (Martin 2012) that we have not tested." A JBES referee will ask: *If $N_{\mathrm{eff}} \le 27$ at $V=400$, does your compounded 4th-order product $A_{v,H}^2 PV_i PV_j$ actually satisfy your finite $(2+\delta)$ score-moment condition (which requires finite $(8+4\delta)$ moments of compounded returns), or is your CLT/fixed-$b$ approximation breaking down because of heavy tails?* Addressing that cleanly (either via explicit tail/lifetime conditions or diagnosing whether it is a finite-horizon kurtosis effect that vanishes as $V$ grows further) is essential for JBES.

### Gap 2 (High Risk): Framing & Connecting to the Econometrics Literature (Not Just DLP12)
- **Too narrow as currently framed:** Right now the Abstract and Section 1 read like *"we diagnose and fix Driessen–Lin–Phalippou (2012) for private equity."* All 10 JBES benchmarks pitch a **general econometric problem** first (here: *identification, finite-pool bias, horizon selection, and weak-ID inference in nonlinear cash-flow discounting / M-estimation for irregular, pooled non-traded cash-flow panels*—covering PE, VC, private credit, real estate, infrastructure, loans/CLOs, project finance) and treat the specific dataset (Preqin Buyout/VC) and baseline estimator (DLP12) as the flagship application.
- **Missing econometrics link on $U_{v,H}$:** A JBES econometrician will immediately recognize dropping own-unit diagonal squares $\sum_{i \neq j} \bar\epsilon_i \bar\epsilon_j$ as related to **jackknife/leave-one-out quadratic-form debiasing** (e.g., Angrist–Imbens–Krueger 1999 JIVE; Kline, Saggio, & Sølvsten 2020; Mikusheva & Sun 2022; U-statistic M-estimation). Right now we don't cite that literature at all. Citing it—and sharply distinguishing *what is genuinely new here* (nonlinear parameter-dependent compounding $A_{v,H}(\theta)^2$, exact duration-matching $H_0^\ast$ under shared/calendar-time shocks in Prop. 2, non-coercive exploding-$\alpha$ boundary geometry, and why cross-vintage LOO collapses to rank 1 and becomes indefinite)—preempts the *"isn't $U_{v,H}$ elementary?"* referee objection.

### Gap 3 (Moderate–High Risk): Manuscript Length & Float Overload
- At **~58–62 double-spaced pages** and **14 main-text floats**, the manuscript substantially exceeds JBES's *"normally 15–30 pp., ideally ≤40 pp. double-spaced"* target (where longer manuscripts face an explicitly higher novelty bar).
- **How to reach ~34–38 double-spaced pages cleanly without losing content:**
  1. Reorder sections to the standard JBES template: `1. Intro (+ Notation block) → 2. Model & Finite-Pool Identification → 3. Optimal Horizon Theory & Increasing-Domain Inference → 4. Monte Carlo Study → 5. Empirical Application: Private Equity Cash Flows → 6. Conclusion`.
  2. Move Preqin-specific data engineering (`Sec 3.2` fund-size imputation details, `Sec 3.3` terminal-NAV Figure 4 & $hv$-block CV details) to the Online Supplement, retaining one tight summary paragraph in Section 5.
  3. Move the two inline proofs (Lemma 1 and Prop. 1) to the Online Supplement (matching 9/10 JBES papers) and trim main-text floats from **14 down to 6–7**.
  4. De-duplicate prose across Intro, Sec 3.4, Sec 4, and Conclusion (right now the same Buyout/VC numbers and MC figures are stated 3–4 times).

### Gap 4 (Easy, Mandatory): JBES Formatting, Abstract, Double-Blind & Replication Rules
- **Abstract:** Trim from 303 words to **≤200 words**, remove `\citet{DLP12}` and dense numerical/LaTeX clutter, and lead with the general econometric contribution.
- **Notation paragraph:** Add a crisp 4–6 line `Notation.` block at the end of Section 1 (present in 6/10 JBES papers).
- **Double-blind anonymization:** Remove/anonymize the identifying GitHub URL `https://github.com/quant-unit/Fundwise_SDF/tree/master/r_project` at `sdf_paper.tex` line 812.
- **ASA / JBES Data & Code Policy:** JBES requires a **Supplementary Materials** section at the end of the paper plus a public replication archive (and, because Preqin is proprietary, an explicit **exemption request to the Editor** paired with public R simulation/synthetic-panel replication scripts—ideally organized cleanly as a lightweight R package or self-contained bundle).

### Gap 5 (Easy, High-Value): Title Refinement (Avoiding "SDF")
All 10 JBES papers use strictly descriptive methodological titles (zero RFS/RAPS question titles and no narrow finance acronyms). Avoiding "SDF" in the title also eliminates the KN16 vs. DLP12 semantic trap upfront and broadens appeal across JBES's general econometrics readership (see `jbes_revision_plan.md` for final non-SDF title candidates).

---

## 4. Strategic Takeaway

- **Compared to RFS / RAPS:** JBES is a **much higher-probability, lower-friction fit**—we do *not* need new Preqin data or an economic-story rewrite; our existing math, Proposition 2, and 72-design Monte Carlo are assets here rather than liabilities.
- **Honest effort estimate to make JBES genuinely competitive:** Roughly **1–2 weeks of focused econometric tightening** (upgrading Appendix D's asymptotic proof from high-level to primitive conditions + clarifying the 1D fixed-$b$ / $V=400$ tail behavior; reframing Abstract/Intro around general cash-flow panel econometrics + LOO/U-statistic literature; reordering and trimming main text to ≤38 double-spaced pages with 6–7 floats; anonymizing & packaging replication code).
- **Fallback synergy:** Every single one of those upgrades also strengthens the paper for **JFEC**, **Journal of Empirical Finance**, or **Econometric Reviews / Econometrics Journal** if JBES passes—nothing is wasted.

---

## Appendix: The 10 JBES Benchmark Papers Audited in This Folder

1. `202516.pdf` — Mehrabani & Parsaeian (2025 WP, JBES): *Shrinkage Estimation and Identification of Latent Group Structures in Panel Data with Interactive Fixed Effects* (36 pp. main + 85 pp. Online Supplement; 6 Theorems/Corollaries, all proofs online; 2.5-pp. PWT empirical app).
2. `A Statistically Identified Structural Vector Autoregression with Endogenously Switching Volatility Regime.pdf` — Virolainen (2025, JBES 43(1), 44–54): 11 printed pp.; all proofs + MC online; highlights CRAN R package `gmvarkit`.
3. `Combining Instrumental Variable Estimators for a Panel Data Model with Factors.pdf` — Harding, Lamarche, & Muris (2025, JBES 43(3), 684–695): 12 printed pp.; Motivating Example + Fig. 1 in Intro; builds on Hansen & Lee (2019) while deriving new joint distribution & combination theorems.
4. `ConsiderationSetsProject_july_upload.pdf` — Aguiar & Kashaev (2023 WP, JBES): *Identification and Estimation of Discrete Choice Models with Unobserved Choice Sets* (35 pp. main; proofs + entire MC online).
5. `Estimating Latent-Variable Panel Data Models Using Parameter-Expanded SEM Methods.pdf` — Wei (2025, JBES 43(2), 324–337): 14 printed pp.; Toy Model in Sec. 2; zero real-data empirical application.
6. `Estimation-and-Inference-for-a-Semiparametric-Time-Varying-Panel-Data-Model.pdf` — Fei Liu (2025, JBES 43(4), 956–967): 12 printed pp.; U.S. CRSP mutual fund empirical application ($N=177, T=143$, 3 printed pp.).
7. `IdentificationofLatentSubgroupsforTime-VaryingPanelDataModels1.pdf` — He, Luo, Liu, Mao, & Zhou (2026, JBES 44(3), 991–1003): 13 printed pp.; Notation block; all proofs online.
8. `Instability of Factor Strength in Asset Returns.pdf` — Massacci (2025, JBES 43(4), 910–925): **Closest asset-pricing benchmark** (16 printed pp., 4 tables + 2 figs; Banking & Finance Dept, King's College London; cites Burnside 2016 JFEC; builds on Bailey et al. 2021 & Qu–Perron 2007).
9. `ISBA_RP_2026-05.pdf` — Hafner, Herwartz, & Wang (2026 reprint, JBES): *Statistical Identification of Independent Shocks with Kernel-Based Maximum Likelihood Estimation and an Application to the Global Crude Oil Market* (34 pp. main; re-examines Kilian 2009 oil data).
10. `ssrn_id4126360_code678999.pdf` — Wang, Weiergraeber, & Xiao (2022 WP, JBES): *Identification of Dynamic Discrete Choice Models with Hyperbolic Discounting Using a Terminating Action* (35 pp. main + 9 pp. App; zero empirical dataset).
