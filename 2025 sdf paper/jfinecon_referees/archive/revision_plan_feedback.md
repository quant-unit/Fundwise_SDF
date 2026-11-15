# Feedback on `revision_plan.tex` (draft of 15 Nov 2026)

**Sources checked:**
- `revision_plan.tex`
- both referee reports
- the working paper `tex/sdf_paper.tex` and `tex/sdf_paper_online_appendices.tex`
- the cited R code lines
- the cited papers: DLP12, KN16, Harris–Jenkinson–Kaplan 2014, Kaplan–Schoar 2005, Robinson–Sensoy 2016

Line numbers below refer to those files. Appendix numbering follows the separately compiled online appendix (JFEC build), which is the version the referees read.

**What is right:** the diagnosis and the order of work (code repairs B1–B9 before any rerun). The exact rescaling identity $\bar\epsilon_{v,H}=A_{v,H}PV_v$ is the right centerpiece.

**Plan claims verified against the paper and code:**
- B3 and B4 bugs: `estim_model_optimized.R:511`, `:600`; `prep_preqin.R:229`.
- Lemma 1(3) uses $\Psi_{t^{\mathrm{Future}}}/\Psi_\tau$, which discounts rather than compounds (`sdf_paper.tex:255`).
- Sample counts for vintages ≤ 2010 are 473 Buyout funds over 25 vintages and 528 VC funds over 27 vintages, versus the "894/931" stated at `sdf_paper.tex:122`.
- "Mathematically immune" (`:888`), the NAV co-movement sentence (`:997`), the monthly-frequency claim (`:903`), and the "heuristic" passage and equicorrelation remark in Appendix B (appendix `:526`, `:558`).

Items are ordered by importance. Items 1–9 are correctness issues; items 10–14 are judgment calls, labeled as such.

## 1. Positioning of $U_{v,H}$ under shared vintage shocks

**Plan says:** the proposed abstract (Section 8.2) claims $U_{v,H}$ "recovers the structural SDF parameters with as few as two funds per vintage."

**Problem:** T7 gives consistency for $\theta_0$ only when $c_v\equiv 0$. By T3 and T6(b), with $c_v>0$ the limit criterion $\E[A_{v,H}^2(\mu_v^2+c_v)]$ generally has a nonzero gradient at $\theta_0$. So both infinite-pool least squares and $U_{v,H}$ are inconsistent for $\theta_0$ even when (P) and (X) hold. DGP (b) is built to satisfy (P) and (X) with $c_v>0$, so the planned simulation will show $U_{v,H}$ biased, which contradicts the abstract. Section 1 of the plan already says shared vintage errors are the binding empirical limit.

**Change:**
- Describe $U_{v,H}$ as removing the $O(1/n_v)$ fund-level noise and reaching the infinite-pool least-squares target with $n_v\ge 2$. Section 5's "Why this framing" item 1 already uses this accurate wording.
- Name the estimand $\theta_H^\dagger$ explicitly.
- State that under $c_v>0$, only zero-mean moment or IV criteria stay consistent for $\theta_0$ as $V\to\infty$, subject to root selection.

## 2. The β-bias direction is contradicted by the paper's own Appendix D; cite DLP12 p. 522

**Plan says:** the $H=0$ noise tilt biases "$\hat\alpha_0$ and $\hat\beta_0$ upward," and this is "the exact mathematical root of the 'exploding alpha' problem in DLP12" (Section 1, item 2.2; T4; T5).

**Evidence:**
- The paper's Online Appendix D, `tab:driessen_comparison` (appendix `:1293–1326`), simulates $H=0$ with α free on vintage-year portfolios. In all ten scenarios S1–S10, mean $\hat\beta$ lies between 0.62 and 0.93 (true value 1) and mean $\hat\alpha$ between +0.13 and +0.71 %/month (true value 0).
- DLP12 report the same direction (JFQA, p. 522): with very high idiosyncratic volatility, "an α that is too high on average and a β that is too low." They attribute it partly to the fact that "a high α helps to decrease the absolute size of the NPVs." The bias shrinks as the number of projects per fund rises from 15 to 50 to 100.
- The main-text single-fund simulation (α + MKT, `sdf_paper.tex:1024`) instead shows $\hat\beta\approx 1.5$ at short horizons.
- Caveat: both of the paper's simulations were run with the simulator that B1 says is defective, so they are suggestive, not conclusive. DLP12's result is unaffected by B1.
- The robust point is analytic: a negative β-score at $\theta_0$ does not determine the sign of the β bias when α is estimated jointly (see change (b)). Rerun both simulations after B1 to settle the sign.

**Changes:**
- (a) Claim an upward β tilt only for β estimated alone (α fixed).
- (b) For the joint case, use the Taylor expansion the paper already has in Appendix B (`eq:theta_shift`, appendix `:400–407`): $\theta^\dagger-\theta_0\approx-[\nabla^2Q]^{-1}\nabla Q(\theta_0)$. The sign of the β component depends on the α–β off-diagonal of the Hessian. Compute it from the data; the gradient correlation is already planned in Section 10, item 4.
- (c) The paper already credits DLP12 for the degeneracy itself (`sdf_paper.tex:424`, `:914`; appendix `:1283` cites their footnote 8). What it does not cite is the p. 522 passage, which states the noise-tilt mechanism qualitatively. Cite it.
- (d) Frame the contribution as the exact decomposition (T3/T4) and the $H>0$ coercivity result, not as the discovery of the mechanism. Present the $H=0$ limit in T5 as a lemma.

## 3. Rerun Appendix D (the DLP replication) after B1 and link the two

Appendix D states that the paper "was unable to reproduce the favorable small-sample properties reported in DLP12": $\hat\beta$ is 0.62–0.93 against DLP's 0.98 (appendix `:1289–1290`). The plan's B1 identifies two simulator defects that bias exactly this kind of exercise:
- deal payoffs compound over `start:end`, which creates a pricing error of $G_a(\theta_0)-1$ even with zero noise;
- Gaussian shocks followed by default truncation raise the conditional mean.

The plan never links B1 to Appendix D.

**Change:** rerun `tab:driessen_comparison` after B1. If the gap to DLP closes, rewrite the "unable to reproduce" paragraph; it is currently an unforced negative claim about a published paper. If the gap remains, B1 is ruled out as the explanation and the paragraph can stay.

## 4. $U_{v,0}$ keeps the exploding-alpha degeneracy at $H=0$

**The argument:**
- Under the predetermined origin, every vintage has $CF_{v,t_v^0}=0$ (plan, evidence row "Inception convention").
- As α→∞, every fund's $PV_{v,i}\to 0$, so $U_{v,0}\to 0$.
- Since $A_{v,0}=1$, $\E[U_{v,0}(\theta_0)\mid\C_v]=c_v(\theta_0)\ge 0$.
- So the α→∞ limit beats $\theta_0$ whenever $c_v>0$. When $c_v=0$ (so the population criterion is $\E[\mu_v(\theta)^2]$), $\theta_0$ stays the unique minimizer on a compact $\Theta$ under the usual identification condition ($\mu_v(\theta)\neq 0$ for $\theta\neq\theta_0$). But the criterion tends to the same value 0 toward large α, so α is weakly identified near the bound.
- This is the same degeneracy as squared loss: only $H>0$ provides coercivity.

**Changes:**
- Say this explicitly.
- T5 covers $Q_H$ only. Add a coercivity result for the sample $U_{v,H}$ with $H>0$. Sketch: if every fund's first nonzero flow is a capital call, the dominant pair products are positive and $U_{v,H}\to+\infty$. This needs a proof, because $U_{v,H}$ can be negative.
- Appendix §3 of the plan (VC hits the +2%/mo bound under both criteria at $H=0$) is consistent with this, but does not state which origin convention it used; state it.
- Under the sample-origin convention, the α→∞ limit of $U_{v,0}$ is the average off-diagonal product of origin-date flows. It is zero unless at least two funds have a flow on the origin date.

## 5. The argument for rejecting the mean moment

**Plan says:** reject the mean moment as the main estimator because with α=0 it gives $\hat\beta=1.86$ in Buyout and a second root at 3.00 (triage row R2-3.1; T8).

**Problem:** suppose the model were correct and the only difference between the criteria were squared loss's upward variance tilt (upward for β with α fixed at 0, by the plan's own T4). Then the mean-moment root should lie *below* the least-squares estimate on the same sample (1.21), not above it. A root of 1.86 points to misspecification with α=0 imposed and/or sampling noise with $V=25$, not to a defect of mean moments. Pooled Buyout NPVs are positive for 1996–2005 at α=0 (plan, Appendix §2). Multiple roots are a valid concern, but they can be handled by root selection.

**Change:** present 1.86 as a diagnostic of misspecification, not as evidence that the estimator fails.

## 6. What in the paper is actually false at $H=0$; Corollary 1 is conditional, not false

**Plan says:** Corollary 1 "falsely claimed $H=0$ consistency for squared loss under bounded $n_v$" (Assumption A4(iii)); W2 says "delete bounded-$n_v$ Corollary 1."

**Problem:** `cor:consistency_structural` (appendix `:993–1001`) is stated *conditional on* $Q_0$ being uniquely maximized at $\theta_0$. As stated it is true. It is vacuous here because that condition generically fails with bounded $n_v$, which is R2's point. The statements that are false as written are the unconditional ones in the main text:
- `sdf_paper.tex:435`: at $H=0$ the estimator "targets the structural parameter $\theta_0$";
- `:454`: $\theta_H^\dagger$ "coincides with the structural parameter $\theta_0$ only when $H=0$";
- `:1126`: $H=0$ "yields unbiased estimates under the exact-pricing restriction."

**Change:** fix those three sentences. Either drop the corollary, or keep it and state that its identification condition fails generically under bounded $n_v$ (T3/T4 show why). Don't describe it as a false theorem.

## 7. Proof and labeling issues the plan misses

- **Consistency proof (appendix `:974–985`).** The step "there exists a further subsequence … such that $\hat\theta_H\to\bar\theta$ almost surely" is not valid for a random sequence as written. Replace it with the standard identifiable-uniqueness argument: continuity, compactness and a unique maximizer give $\sup_{\|\theta-\theta^\dagger\|\ge\varepsilon}Q_H(\theta)<Q_H(\theta^\dagger)$. Then $P(\|\hat\theta-\theta^\dagger\|\ge\varepsilon)\le P(Q_H(\theta^\dagger)-Q_H(\hat\theta)\ge\eta)\to 0$. Alternatively, cite NM94 Theorem 2.1 directly.
- **"Theorem 2 in Appendix C" does not exist** (plan, evidence row 3 and W2). The result missing the score-centering condition and the integrated Hessian is the CLT, `prop:clt_H`, stated in the main text and proved in the appendix.
  - The proof uses a single intermediate point $\tilde\theta_H$ for a vector-valued score (appendix `:1031–1040`). That is exactly where the integrated Hessian $\int_0^1\nabla^2Q_{n,H}(\theta^\dagger+s(\hat\theta-\theta^\dagger))\,ds$ is needed.
  - Score centering is needed because `as:population_limit_H` gives only $o(1)$ convergence of the mean criterion, not $o(n^{-1/2})$.
- **"Remark 2" is ambiguous.** In the main text, Remark 2 is `rem:cv_diagnostic`. The equicorrelation remark is Remark 2 only in the separately compiled appendix (`rem:equicorr`). In the working-mode build, where the appendix is `\input` into the main file, it is Remark 4. Refer to all results by label so the drafting step deletes the right remark.

## 8. Table 1 of the paper contains false entries; the KN16 entry is correct

- `sdf_paper.tex:814–822`: the "Our approach" column says "Both (tested empirically)" for the SDF specification (linear and exponentially affine) and for the cross-sectional unit (single fund and vintage-year portfolio). The paper estimates only the linear SDF (Eq. `eq:linear_sdf`) and, empirically, only vintage-year portfolios; single funds appear only in the simulation. Correct both cells.
- Keep "exponentially affine" for KN16. KN16 Eqs. (1) and (8) use $M_{t+1}=\exp(a-b\,r_{m,t+1})$ and call it exponentially affine. Referee 1's description "$M=a+bR$" is inaccurate, so the Table 1 revision should not adopt it.
- On "without strong distributional assumptions" (`sdf_paper.tex:112`): KN16 use almost exactly this phrase about their own method in their introduction. R2's objection applies to KN25, which adds i.i.d. joint lognormality. Attribute the phrase to KN16 and state that KN25 adds structure; don't simply delete it.

## 9. The plan's empirical numbers come from a different sample than the paper's

The paper's sample is North America, USD (`sdf_paper.tex:847`). The plan's probes (Section 1 bullets, the evidence rows tagged [D], Appendix §§1–3) use the stored **global** fund-level file. Examples:
- the "1.23 → 0.62 versus 1.21 → 0.96" Buyout horizon profiles;
- the $U_{v,H}$ estimates;
- the vintage-residual runs;
- the 1.86 mean-moment root.

The paper's own North America tables (`tab:empirical_bo_mkt`, `tab:empirical_vc_mkt`, still under the B4 pooling bug) show Buyout EW 1.24 → 0.68 and VC EW 1.88 → 1.05.

**Change:** label every probe number as global-sample. Don't carry any of them into the revised text, abstract or literature table until Gates G1/G4 have rerun them on the North America constant-composition pools. The qualitative conclusions (pooling fix flattens the Buyout profile; shared vintage residuals dominate) are hypotheses for the North America sample until then.

## 10. Gate G3 (judgment)

**Plan says:** include the IV-GMM estimator only if its β RMSE is no worse than $U_{v,H}$ **and** its bias is materially smaller.

**Problem:** in DGP (b) the IV estimator is the consistent one, and with $V=25$ it will probably have higher variance. The RMSE condition is then likely to exclude the only consistent estimator by design.

**Change:** drop the "RMSE no worse" requirement, or report the IV estimator in any case as the consistent comparator, judged on bias and coverage.

## 11. SDF framing (judgment, strong recommendation)

**Plan says:** keep "SDF" in the title and text. It says the referees "objected only because the submitted text blurred three distinct objects," and that the distinction "answers both referees completely" (Section 1).

**Problems:**
- (a) Being an SDF "on the PE cash-flow space" holds by construction for any strictly positive deflator that satisfies the PE pricing moment the estimator imposes. For a restricted payoff space, the set of such deflators is large (Hansen–Richard). The label adds no restriction beyond the moment condition being estimated.
- (b) The proposed abstract calls α "abnormal performance." If $\Psi(\theta_0)$ prices PE cash flows, PE has zero abnormal performance relative to $\Psi$ by construction. α is abnormal performance relative to the benchmark return $1+r_f+\beta^\top F$, which is DLP's return-model α. The plan's own Result 3.2 shows that an economy-wide SDF in this family forces α=0.
- (c) Both referees made the point explicitly. R1: "The DLP-2012 approach does not use the SDF approach." R2 quoted KN25: "$1/R^b_h$ is not an SDF, as it depends on an asset-specific parameter (β)." The paper's Table 1 also lists DLP12 under "SDF specification: linear factor model," but DLP12 never use the terms "SDF," "stochastic discount factor," "discount factor" or "pricing kernel." They present a return model for abnormal performance and risk exposure. The referee pool at the next journal will overlap.

**Change:** pick one of two options.
- (i) Call the estimated object a benchmark deflator, use SDF language for motivation and the literature, and consider removing "SDF" from the title.
- (ii) If "SDF" stays, add a precise scope statement in Section 2 of the revised paper: the deflator prices PE cash flows only; it does not price the T-bill or the market; and α is read relative to the benchmark return. Also delete "answers both referees completely."

"Stability" can stay in the title either way: T5's coercivity result is a stability result.

## 12. Wording and overclaims

- **"Convexity":** the proposed abstract says horizon averaging "restores global convexity in α." T5 proves coercivity ($Q_H\to+\infty$ as α→∞), not convexity. The paper uses the same loose wording ("not globally convex," `sdf_paper.tex:116`, `:425`). The actual problem is a degenerate infimum toward large α, not a lack of convexity. Use "coercive in α" or "rules out the degenerate α→∞ solution" in both places.
- **"Stable betas":** the proposed abstract claims "stable single-factor market betas." VC is not stable on any evidence available now:
  - paper, North America, pre-fix: EW 1.05–2.00 and FW 0.82–1.68 across horizons;
  - plan probes, global: 0.83–2.32 under squared loss and 0.95–2.10 under $U_{v,H}$, falling to 0.01 and −0.50 under the predetermined origin at $H=30$.
  - Restrict the claim to Buyout, or write the abstract after Gate G1.
- **Tone:** drop "bulletproof" (Section 5), "100% mathematically correct" (Step 2) and "answers both referees completely" (Section 1). None of these should carry into the paper or a cover letter.

## 13. Interpretation of the vintage residual cycle (judgment)

**Plan says:** the residual runs (global sample: Buyout + for 1996–2005 and − for 2006–2010; VC + for 1992–1998 and − for 1999–2010) are read as "$c_v>0$ and/or $\mu_v^2>0$."

**Add the literature:**
- Harris, Jenkinson and Kaplan (2014) find that both absolute and PME performance are negatively related to aggregate capital commitments, for Buyout and VC.
- Kaplan and Schoar (2005): funds raised in boom times appear to perform poorly.
- Robinson and Sensoy (2016) find that hot-market underperformance is sharply attenuated relative to public equities. HJK also note that Robinson and Sensoy find no such relation for buyout PMEs.
- So the evidence is mixed after risk adjustment.

**Implication:** aggregate fundraising at the vintage date is in $\F_{t_v^0}$. If it predicts risk-adjusted NPVs, then (P) fails for a constant-α model. That is misspecification, not a shared shock that (P) permits.

**Suggested diagnostic:** if vintage-level aggregate commitments are available, regress $PV_v(\hat\theta)$ on them.

## 14. Inference with $V\approx 25$ (judgment)

The paper's own SHAC theory requires $b_V\to\infty$ and $b_V/V\to 0$ (appendix `eq:bandwidth_conditions`), yet it uses $b_V=12$ years with about 25 vintages (`sdf_paper.tex:1144`), so $b_V/V\approx 0.5$. The plan itself says there are fewer than four effective independent macro cycles.
- Keep the CLT repairs from item 7: a theorem that is stated cannot have a broken proof. If the theorem is not used, dropping it is the alternative.
- Make simulated coverage at $V=25$ the stated basis for inference. It is already listed among the reported metrics.

## 15. Minor

- **DGP (c):** now that the wedge is only a diagnostic, consider moving it to the appendix or cutting it, to limit scope.
- **Public code:** the paper links a public code repository (`sdf_paper.tex:1220`). The B1–B9 repairs should be pushed there once verified, and any posted working-paper version updated, because its empirical horizon profile was partly driven by the B4 pooling bug.
