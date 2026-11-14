**Independent audit of the review of the two referee reports**

14 November 2026. Main focus: Appendix B of the review, including its implications for Appendix B of the paper.

**Assessment.** The review identifies a real and important problem: a zero expected pricing error does not generally make the average squared pricing error stationary at the structural parameter. Its rescaling identity and its new gradient calculation are correct under appropriate assumptions. However, I would not use Appendix B as the revision blueprint in its present form. It overstates what those calculations establish, applies a proposed correction beyond the assumptions of its proof, and incorrectly certifies both the simulation implementation and parts of the paper's existing formal analysis. Appendix A also describes an SDF estimand that the review's own Appendix B says the actual estimator need not recover.

The most consequential new finding is a timing mismatch in the checked-in simulation code. This is distinct from the review's already identified portfolio-aggregation issue, which I reproduced.

I left the original review, paper, and estimation code unchanged. The accompanying [independent checks](</Users/christausch/github/Fundwise_SDF/2025 sdf paper/jfinecon_referees/independent_audit_checks.py>) enumerate finite-support examples and counterexamples. They do not reproduce the paper's empirical estimates.

**Reading guide and scope.** References to “review lines” refer to [the review source](</Users/christausch/github/Fundwise_SDF/2025 sdf paper/jfinecon_referees/review_of_referee_reports.tex>). Review PDF pages 39–48 contain Appendix B. I read the full review and both reports, checked the relevant paper text and proofs against the supplied April 26 PDF and the current LaTeX sources, inspected the relevant R implementation, independently checked the central formulas, reran selected replica calculations, and checked the principal disputed literature passages. This is not a complete replication of the paper or a certification of its SHAC implementation. Section and proposition numbering differ between the combined April PDF and the referee-facing sources; the source links below avoid that ambiguity.

| Priority | Finding | Consequence |
|---|---|---|
| Essential | The simulation includes the investment month's return in the exit payoff; the estimator's investment-to-exit ratio excludes it. | The claim that the actual simulation satisfies the exact zero-wedge model is unsupported and, for this code, false at the stated parameters. |
| Essential | The actual simulation estimator extends the observation grid beyond the last cash flow; the replica does not. | The replica cannot establish the mechanism or optimal horizon for the published long-horizon results. |
| Essential | The pair-product correction needs equal conditional means, not just independent fund shocks. | Different fund schedules can make the proposed population objective negative away from the truth, even with perfect pricing and no noise. |
| Essential | The “all formal results remain correct” verdict overlooks conditioning, compounding-ratio, and CLT-centering problems. | Some formal statements or assumptions need repair, not just surrounding prose. |
| Essential | Appendix A equates the estimated object with a member of the PE-pricing set, despite Appendix B's target shift. | Keeping “PE-specific SDF” does not resolve the estimand problem for the current estimator. |
| Important | “No horizon is consistent,” monotonic beta decline, and an order-\(1/n_v\) parameter shift are stated too broadly. | Replace these with conditional and local results; separate score, parameter bias, and finite-sample MSE. |
| Important | A flat product-criterion horizon profile is proposed as a specification test without sampling theory. | At most it is a population implication and a possible diagnostic pending valid joint inference. |
| Confirmed | Changing composition in empirical vintage portfolios, quarterly cash-flow dates, and horizon saturation. | Keep these findings, but do not promise the empirical conclusions remain unchanged before rerunning the corrected construction. |

**1. What in Appendix B is correct and worth preserving**

The identity is exact. For a fixed cash-flow series and an admissible nonzero discount-factor process,

\[
\epsilon_{\tau,i}(\theta)=PV_i(\theta)\frac{\Psi_{t_i^0}(\theta)}{\Psi_\tau(\theta)},
\qquad
\bar\epsilon_{i,H}(\theta)=PV_i(\theta)A_{i,H}(\theta).
\]

It applies to all revaluation dates because the whole cash-flow stream is used at each date. The implementation in [getNPVs.R](</Users/christausch/github/Fundwise_SDF/r_project/helper/getNPVs.R:12>) has exactly this structure. The strict positivity of \(A\), however, needs positive benchmark gross returns; the statement at review line 1193 that this “needs no assumption” should distinguish the algebra from positivity and integrability.

The conditional-mean argument is also correct. If \(A_{i,H}(\theta_0)\) is measurable with respect to an environment \(\mathcal G\), and

\[
E[PV_i(\theta_0)\mid\mathcal G]=0,
\]

then \(E[A_{i,H}PV_i]=0\). Independence across deals is unnecessary for this particular conclusion. It becomes relevant for pooling rates and the simplified variance formula.

Result 18, the vector gradient formula at review lines 1210–1228, is correct for the stated unit-investment, single-exit, conditionally independent deal model, provided differentiation and expectation may be interchanged and the horizon window is fixed with respect to the trial parameter. In particular,

\[
\nabla Q_H(\theta_0)=\frac{2}{n_v^2}
E\left[A_{v,H}^2\sum_j c_j^2\nu_j
\left(\bar S_{v,H}-S_{t_v^0,b_j}\right)\right].
\]

I independently enumerated all 16 outcomes of four independent two-point deal shocks, with staggered investments and exits and both alpha and beta parameters. Numerical derivatives agreed with the formula to less than \(2\times10^{-10}\) at three horizons. This is stronger evidence for the algebra than approximate agreement in a noisy Monte Carlo.

At \(H=0\), with a free alpha, strictly positive gross returns, positive holding durations, and nonzero idiosyncratic variance, the alpha component of this gradient is strictly negative. That is an especially clean result to retain. It establishes a local incentive to increase alpha at the structural parameter. It does not, by itself, establish a global exploding-alpha optimum.

The pair-product identity is exact:

\[
\frac{1}{n(n-1)}\sum_{i\ne k}e_ie_k
=\bar e^2-\frac{s_e^2}{n}.
\]

Under conditional independence and a common conditional mean, its expectation is the square of that mean. I verified both assertions by exact enumeration. The problem is the scope of the application, not the identity.

**2. The checked-in simulation does not implement the review's exact DGP**

Review lines 178, 1162, 1179, 1392 and 1406 repeatedly assert that the paper's simulation satisfies Tier 2 and has no wedge. Two implementation differences need to be resolved first.

**2a. One extra return enters the payoff.** Both [simulation_improved.R](</Users/christausch/github/Fundwise_SDF/r_project/simulation/simulation_improved.R:113>) and the generator used by the current simulation runner, [simulation_multifactor.R](</Users/christausch/github/Fundwise_SDF/r_project/simulation/simulation_multifactor.R:225>), compound `df$deal[start:end]`. They then place the investment at `start` and the distribution at `end`.

The estimator builds cumulative benchmark returns and values the two cash flows using their ratio; see [estim_model_optimized.R](</Users/christausch/github/Fundwise_SDF/r_project/estim_model_optimized.R:326>). The ratio from investment date \(a\) to exit date \(b\) contains returns \(a+1,\ldots,b\). Ignoring default truncation, let \(G_u=1+\alpha_0+r_{f,u}+\beta_0 F_u\). The generated payoff has conditional mean

\[
E[\mathrm{Div}\mid\mathcal G]=\prod_{u=a}^{b}G_u,
\]

whereas the paper's pricing convention requires \(\prod_{u=a+1}^{b}G_u\). Therefore, for a unit investment,

\[
E[\delta(\theta_0)\mid\mathcal G]
=-1+\frac{\prod_{u=a}^{b}G_u}{\prod_{u=a+1}^{b}G_u}
=G_a-1,
\]

which is generally nonzero. This persists as the number of funds grows.

A noise-free example makes the mismatch visible: monthly gross returns 1.1, 1.2, 1.3 generate a payoff of 1.716 under the inclusive convention. Its inception NPV under the estimator is +0.1. Using the subsequent two returns instead gives 1.56 and NPV zero. The independent script checks this exactly.

The review's Python replica uses `(k > start) & (k <= end)` and therefore implements the subsequent-period convention. It silently repairs the mismatch rather than reproducing it. This means its successful zero-wedge checks cannot validate the R simulation.

There is also a distinct inception convention: `max.month == 0` returns the raw discounted sum, while `max.month == 1` multiplies it by the first cumulative return. These are different, parameter-dependent normalizations. The configured simulation grid includes both 0 and 1. Before comparing tables labeled “0 years,” establish which code branch each table uses. For larger values, `getNPVs` averages the first `max.month` observations, while the review's displayed \(0,\ldots,12H\) window has \(12H+1\) observations when uncapped.

**Action:** settle the cash-flow timing convention; align generation and valuation; run a noise-free validation that every deal and pooled vintage has zero inception NPV at the true parameters; then reproduce the affected simulation results. I verified the current source, not the provenance of every cached table. Thus this is a demonstrated problem in the checked-in implementation and a required provenance check for the published numbers, not a claim that I have already measured their numerical error.

**2b. The horizon extends beyond the last cash flow in the actual estimator.** [The monthly-grid constructor](</Users/christausch/github/Fundwise_SDF/r_project/estim_model_optimized.R:186>) sets `max.date <- as.Date("2020-01-01")` when `use.simulation` is true and fills missing monthly flows with zeros. The grid thus normally extends through December 2019. `getNPVs` caps the horizon at the length of this padded series, not the last nonzero cash flow.

By contrast, `pieces` and `pieces_ab` in [check_appendix_horizon.py](</Users/christausch/github/Fundwise_SDF/2025 sdf paper/jfinecon_referees/check_appendix_horizon.py:41>) cap it at the last nonzero pooled cash flow, at most 15 years in that replica. Consequently, the replica's 15- and 20-year cases coincide, while early vintages in the actual simulation can compound errors well beyond their last exit. This directly affects the sign and size of \(\bar S-S_b\).

The review does disclose that it had not checked this endpoint (line 1460). That caveat should be resolved before interpreting the paper's long-horizon decline using the replica. In particular, “largely a tilt signature” and “part ... is a finite-sample effect and not the tilt” are not established for the paper's actual design.

**2c. Default truncation is another departure from exact mean-one noise.** For a Gaussian shock, \(E[(G+e)_+]>G\). Setting a deal's payoff to zero whenever a monthly gross return becomes negative therefore changes its conditional mean. The review acknowledges an approximation at the baseline noise level, but later statements treat it as an exact theorem about the simulation. Separate the idealized no-truncation result from the implemented process, especially in the higher-volatility scenarios.

**3. The cross-fund product theorem is conditional on more than Tier 2**

Result 19 explicitly adds exchangeability and independence given \(\mathcal G\); its proof is valid under those added conditions. But the short answer, consistency table, and empirical recommendation repeatedly reduce the requirement to Tier 2 or independent fund shocks.

For conditionally independent but heterogeneous fund errors, write

\[
\mu_i(\theta)=E[PV_i(\theta)\mid\mathcal G].
\]

Then the actual conditional expectation is

\[
E[U_{v,H}(\theta)\mid\mathcal G]
=\frac{A_{v,H}^2}{n_v(n_v-1)}\sum_{i\ne k}\mu_i(\theta)\mu_k(\theta)
=A_{v,H}^2\left[\bar\mu(\theta)^2-\frac{s_\mu^2(\theta)}{n_v}\right].
\]

This need not be nonnegative. All \(\mu_i(\theta_0)=0\) only establishes a zero value at the truth; it does not establish minimization there.

Here is a counterexample within the benchmark framework. Take two funds, no noise, zero risk-free return, and two factor realizations +0.2 and -0.2. The true beta is one. Fund 1 invests 1 at date 0 and receives 1.2 at date 1. Fund 2 invests 1 at date 1 and receives 0.8 at date 2. Both are priced exactly by the true benchmark. Measured at the common date 0,

\[
PV_1(\beta)=\frac{0.2(1-\beta)}{1+0.2\beta},\qquad
PV_2(\beta)=\frac{0.2(\beta-1)}{(1+0.2\beta)(1-0.2\beta)}.
\]

The cross-fund product is zero at the truth and strictly negative away from it on the positive-return domain. At beta 0.5 it is -0.009183; at beta 1.5 it is -0.008453. Minimization favors a wrong beta. Multiplying by a common positive \(A_H^2\) does not repair this.

This does **not** refute Result 19 with its explicit exchangeability assumption. It refutes treating that assumption as automatic when using heterogeneous schedules or claiming that Tier 2 alone gives the product estimator consistency.

The conditioning field in Assumption A includes the realized deal timings. Once those are conditioned on, otherwise identically generated funds usually have different conditional mean functions away from the truth. Dropping timing from the conditioning field can restore conditional iid structure if timing and shocks are independently generated across funds. But then a horizon window based on the sample's earliest investment and latest exit is not simply a fixed function of that smaller field. The proof has to account for the dependence introduced by that window.

A workable route is to use a predetermined common vintage date and horizon window and explicitly assume that the full fund cash-flow vectors, including their schedules, are conditionally iid given the common factor environment. Then the product proof applies. Alternatively, derive and verify a more general identification condition for heterogeneous schedules. Either route is additional modeling work.

Common vintage shocks also survive the correction. Without conditional independence,

\[
E[U\mid\mathcal G]
=\frac{A^2}{n(n-1)}\sum_{i\ne k}
\left\{\mu_i\mu_k+\operatorname{Cov}(PV_i,PV_k\mid\mathcal G)\right\}.
\]

An alternative criterion should therefore be presented as a promising research extension or a carefully qualified robustness estimator. It is not yet an off-the-shelf fix for Preqin. Fund-size weighting, outcome-dependent reporting, and fund eligibility must enter the assumptions as well. Subtracting the ordinary within-vintage variance is not the weighted correction.

**4. Appendix B's interpretations are stronger than its gradient result**

**4a. “No horizon is consistent” is false as a universal statement.** Review lines 1162, 1272–1273 and the proposed replacement paragraph at 1442 should say that squared loss is **not generally** consistent under bounded pooling. A nonzero noise variance is not itself a proof that every horizon has a nonzero population score. The appendix itself explains how terms can cancel.

For an explicit example, take a one-period payoff \(N\) with mean one and positive variance \(v\), true gross benchmark return one, and average over dates 0, 1 and 2. With candidate gross return \(r>0\),

\[
Q(r)=\frac{(r+1+r^{-1})^2}{9}\left[(r-1)^2+v\right].
\]

Both factors are minimized at \(r=1\), so the structural parameter is the unique minimizer despite positive noise and a single fund per vintage. This example uses a window extending past exit, as the actual simulation implementation permits. Even if one restricts windows to end at exit, special beta-only models can have a stationary variance at the truth. The correct general conclusion is lack of a consistency guarantee, not universal inconsistency.

**4b. The beta sign argument in Section 2 is wrong.** Review line 235 infers the sign of \(E[S\operatorname{Var}(N\mid\mathcal G)]\) from the sign of \(E[S]\). The conditional variance is not constant in the additive-noise model.

For one period, take \(F=0.2\) or \(-0.1\), each with probability one half, \(r_f=0\), \(\beta_0=1\), and additive shock variance 0.04. Then

\[
E[S]=0.027778>0,\quad
E[S\nu]=-0.00042867<0,\quad
Q'_0(\beta_0)=+0.00085734.
\]

The growth-optimal beta is 2.5, so the true beta is below it, yet the score has the opposite sign to the review's claim. Use bounded two-point idiosyncratic shocks of ±0.2 and every gross payoff is positive; the example does not rely on Gaussian tails.

Similarly, \(\partial_\beta A_H\) need not be positive, and \(S_{t^0,t}\) need not increase with time in the beta direction. A simple negative-factor example in the check script makes the fitted beta **rise**, from 0.91 to approximately 0.95287, when the averaging window extends from date 0 to dates 0 and 1.

**4c. A score of order \(1/n_v\) does not automatically imply a parameter bias of that order.** Result 18 is a statement about \(\nabla Q_H(\theta_0)\). To obtain

\[
\theta_H^\dagger-\theta_0=O(1/n_v),
\]

one needs a locally identified limiting signal criterion, a suitably bounded inverse Hessian, uniform moment bounds, and a minimizer staying in that neighborhood. Weak curvature can turn a small score into a large displacement. This qualification is central in a paper whose empirical finding is weak identification of alpha and additional factors.

The “sum has order \(n_v\) terms” argument also presumes a bounded deal count per fund and bounded relevant moments. For unequal portfolio weights, effective diversification depends on concentration, for example \(\sum_i p_i^2\), not solely the head count.

**4d. Do not infer empirical attribution from a falling profile alone.** Review line 1295 and revision items 11–13 suggest that the empirical slope is largely the tilt and that a remaining slope under the product criterion is a specification test. Neither is established. Sampling variation, heterogeneous schedules, changing portfolio composition, NAV errors, parameter constraints, and identification strength can all affect the profile. Comparing estimates from different horizons requires their joint sampling distribution. A non-flat finite-sample profile is not, by itself, rejection of Tier 2; a flat one does not establish it either.

**4e. “Nothing to average out” needs a conditional qualifier.** The dates contain no independently observed residual conditional on the factor path. Nevertheless, unconditional errors \(A_\tau PV\) need not be perfectly correlated, and averaging their random weights can change unconditional variance. Thus the identity invalidates an explanation based on independent date-specific noise, but does not prove that averaging can never reduce variance. The earlier distinction in review line 156 between “no new residual” and “no new moment” is the right one. Carry it through to the abstract draft at line 1447.

Likewise, changing parameter-independent weights across heterogeneous units can change estimator efficiency. The statement at review line 525 that the effect runs through the parameter-dependence of \(A\) alone is too strong. A common constant rescaling cancels; unit-specific weights generally do not.

**5. The review should not certify all of the paper's formal results unchanged**

There are at least three substantive qualifications, in addition to routine regularity conditions.

**5a. The paper's conditional wedge lemma needs an extra assumption.** The paper defines \(b_{H,i}\) conditional on fund structure: deal count, investment dates and the horizon set. Its proof then uses the original investment-date pricing restriction to discard the conditional mean of the deal error. See [the definition](</Users/christausch/github/Fundwise_SDF/2025 sdf paper/tex/sdf_paper_online_appendices.tex:175>) and [the lemma](</Users/christausch/github/Fundwise_SDF/2025 sdf paper/tex/sdf_paper_online_appendices.tex:298>).

From \(E[\delta\mid\mathcal F_a]=0\), one cannot generally infer \(E[\delta\mid\mathcal F_a,\mathcal D]=0\) when \(\mathcal D\) includes future endogenous timing or selection information. Conditional covariance decomposition then contains the extra term

\[
E[X\delta\mid\mathcal D]
=\operatorname{Cov}(X,\delta\mid\mathcal D)
+E[X\mid\mathcal D]E[\delta\mid\mathcal D].
\]

For example, let \(\Psi\equiv1\), invest one, and pay 1.5 at date 1 in the high state or 0.5 at date 2 in the low state, with equal probabilities. Inception pricing holds. For a window capped at the exit, conditioning on the horizon set reveals the state. The conditional mean error is ±0.5, while every covariance with \(\Psi_a/\Psi_\tau=1\) is zero. The displayed conditional lemma fails in this setting.

It is valid for a fixed exogenous design with the pricing restriction maintained conditional on that design. State that explicitly, or keep the unconditional formulation and handle random windows inside expectations. Result 17 in the review uses an **unconditional** \(b\), which is a useful reformulation but not literally the same object as the paper's conditional \(b_{H,i}\).

The review's two-tier table has a related gap: \(\mathcal G\), generated by factors and timing, need not contain all investment-date information. Mean zero conditional on \(\mathcal G\) is not automatically stronger than mean zero conditional on \(\mathcal F_a\). Appendix A notices this in Result 15; Appendix B should preserve the additional conditioning or independence condition. Random investment amounts also need to be included in the relevant environment or covered by a conditional moment assumption.

Nor does the inverse-benchmark functional form in Eq. (19) itself imply Tier 2. That is a restriction on the joint distribution of payoffs, timing and factors, not a consequence of calling the discount factor a benchmark. A zero wedge is a consequence of Tier 2, but is not equivalent to it: covariance can vanish through cancellation even when conditional pricing fails.

**5b. Lemma 1(3) is an identity for a differently weighted variable, not the stated compounding operation.** In the supplied paper, PDF page 6, and [the current source](</Users/christausch/github/Fundwise_SDF/2025 sdf paper/tex/sdf_paper.tex:258>), the NFV-type object is defined as

\[
m^{\mathrm{NFV}}(\tau,f)=\frac{\Psi_f}{\Psi_\tau}\delta.
\]

But \(\Psi\) is a discount process, and \(\delta\) is measured at the investment date \(a\). Revaluing that error forward to \(f\) uses

\[
\mathrm{NFV}_f=\frac{\Psi_a}{\Psi_f}\delta.
\]

The covariance proof is algebraically correct for the variable the paper defines. Its description as compounding the investment-date error, and the claimed direct connection to the estimator, use the opposite ratio. For the actual forward revaluation, the covariance at inception is \(\operatorname{Cov}(\Psi_a/\Psi_f,\delta\mid\mathcal F_a)\). The review's own identity makes this distinction transparent. Rewrite this part or explicitly present it as a different weighted moment; do not preserve the current interpretation unchanged.

**5c. Population stabilization alone does not justify the stated centered CLT.** The paper assumes uniform convergence of the mean criterion to a limiting \(Q_H\), but supplies no rate. Its CLT proof then asserts a zero-mean score limit at the limiting maximizer. See [population stabilization](</Users/christausch/github/Fundwise_SDF/2025 sdf paper/tex/sdf_paper_online_appendices.tex:763>) and [the score-CLT step](</Users/christausch/github/Fundwise_SDF/2025 sdf paper/tex/sdf_paper_online_appendices.tex:1052>).

Consider independent bounded errors \(\varepsilon_i=\pm1\), deterministic \(a_i=i^{-1/4}\), and criterion contributions

\[
Z_i(\theta)=(\theta-a_i-\varepsilon_i)^2
\]

on a compact parameter interval containing zero. The limiting mean criterion is \(\theta^2+1\); it has a unique minimizer, bounded moments and nonsingular Hessian. The variables are independent, and squared loss is Lipschitz on their bounded relevant support. Yet

\[
\widehat\theta=\bar a_n+\bar\varepsilon_n,\qquad
\sqrt n\,\bar a_n\sim \tfrac43 n^{1/4}\to\infty.
\]

Consistency for zero holds; a centered \(\sqrt n\) normal limit at zero does not. The check script shows the deterministic centering term increasing from 4.15 at \(n=100\) to 42.16 at \(n=10^6\).

Add a score-centering rate such as \(n^{-1/2}\sum_iE[s_i(\theta_H^\dagger)]\to0\), or center at a finite-array population target and control its distance from the limit. Also state the score differentiability, score moments, mixing/NED rates, interiority and Hessian conditions actually needed. The current prose invokes several of these without listing them as explicit assumptions. This does not invalidate standard extremum theory; it invalidates an unqualified certification of the current statement.

There are two smaller formal repairs. The vector Taylor expansions in the paper's target-shift and CLT proofs should use an integrated Hessian rather than assert one common intermediate point for a vector-valued gradient. Also, once expectations are taken over random fund structures, the residual term in the conditional-mean gradient decomposition is an expected conditional covariance, not necessarily the unconditional covariance printed underneath it.

**6. Appendix A does not settle the SDF objection for the actual estimator**

The Hilbert-space representation in Result 14 is essentially correct for unrestricted square-integrable linear pricing representations, assuming the price functional is continuous. But:

* The law of one price is not itself continuity on an infinite-dimensional payoff space.
* The unrestricted Riesz representation does not establish a strictly positive discount factor. Appendix A switches from an unrestricted \(L^2\) definition to positive adapted processes in the multiperiod discussion.
* Public-asset pricing is conditional in the multiperiod definition. The unconditional, one-factor, constant-risk-free calculation in Result 6 does not establish Result 16(c)'s general “if and only if” without corresponding assumptions on conditional factor distributions and an admissible interior growth-optimal portfolio.
* Non-uniqueness relative to a restricted payoff set does not imply that a particular finite-dimensional family contains a pricing member, nor that the estimator recovers it.
* The statement that a common discount factor for buyout and venture cannot be in the one-factor family is too strong. In the review's own zero-alpha setting, an admissible public growth-optimal member can price both; its beta need not equal either segment's exposure.

The final point is decisive here. Appendix B says the current squared-loss estimator generally targets \(\theta_H^\dagger\ne\theta_0\), even under the favorable DGP. Appendix A nevertheless says that the estimated discount factor prices PE at zero error and that squared loss selects \(\theta_0\) (review lines 1029, 1094–1096 and 1143). Those assertions cannot both be maintained.

A one-period example is enough: let the payoff be \(R_0N\), with \(E[N]=1\) and variance \(v>0\), and estimate a positive gross discount return \(R\) by minimizing \(E[(R_0N/R-1)^2]\). The minimizer is \(R=(1+v)R_0\). Its expected pricing error is \(-v/(1+v)\), not zero. Thus the estimated population limit need not be in the PE-pricing set whose definition is used to defend the name.

The safe distinction is between (i) a structural benchmark parameter that prices cash flows under stated assumptions and (ii) the criterion-specific parameter estimated from noisy pooled cash flows. “PE-specific pricing kernel” can describe the former if carefully defined. It does not automatically describe the latter. My preference is to lead with estimation of benchmark exposures and horizon-dependent criterion targets. If “SDF” is retained, the title and abstract still need this distinction, not only three terminology concessions.

For citation accuracy, [Hansen and Jagannathan's original paper](https://larspeterhansen.org/wp-content/uploads/2016/11/ImplicationsofSecurityMarketData.pdf) treats restrictions from asset prices and explicitly exploits positivity. It is not a citation that, by itself, establishes the review's particular dynamic PE-only interpretation. The Riesz proof can stand on its own with its exact assumptions; do not claim that the terminology question has been mathematically settled.

**7. Literature and referee-response judgments**

Several of the review's defenses are justified; others need narrower wording.

| Comment | Audit judgment and suggested response |
|---|---|
| R1-1 / R2-1a: SDF versus benchmark | Accept the substantive distinction. Do not rely on Appendix A's naming argument to dispose of it; separate the structural pricing object from the actual criterion limit. |
| R1-2: conditions for consistency | The review's all-horizon infill argument is useful, with a conditional LLN, identification, positivity, and uniform convergence. Call its assumptions sufficient DLP-type conditions, not necessary conditions. |
| R1-3: compare methods | Different estimands rule out a naive beta-bias comparison, but do not rule out a common performance or valuation comparison under an explicitly specified DGP. The proposed KN16-style pricing-error readout is constructive. Do not dismiss the request as mathematically impossible. |
| R1-4: empirical comparison | Keep the recommendation for a compact literature comparison, verifying samples, gross/net returns, and definitions. The review has not yet supplied those verified numbers. |
| R1-5: exposition | Valid. Consolidate the repeated results and remove contradictory recommendations. |
| R2-1b: KN16 | The review is right about the disputed sentence: KN16 p. 1447 discusses the need for many non-overlapping windows when estimating GPME. It is not a contradiction of public-side kernel calibration. Still repair the manuscript's ambiguous use of “identification.” |
| R2-1c: distributional assumptions | The review is right that the referee's KN25 example does not directly refute a statement about DLP12/KN16. Still distinguish distributional assumptions from conditional mean, independence, timing, and homogeneity assumptions. |
| R2-1d: literature engagement | Accept. KN25 is relevant, but its constrained dispersion method, benchmark construction and SDF-based performance level make it more than the same criterion with one constraint removed. Describe the relationship as a useful structural analogy. |
| R2-2a / 3.1 / 3.4 | The square-versus-mean distinction and rescaling identity are strong answers. Neither establishes a general finite-sample efficiency improvement. That empirical claim still needs matched simulations. |
| R2-2b: horizon choice | Treat horizon as a sensitivity parameter until a selection rule is justified. If cross-validation is used, score competing horizons with a common loss and an appropriate dependence-aware design. |
| R2-2c: Merton | The general difficulty of estimating means is relevant. The review's chosen 15–25% residual volatility and 30 effective years are illustrative, not an empirically established lower bound. Fold dispersion is not automatically an estimator standard error. |
| R2-3.2: GSW | The two mechanisms are not identical. Under conditional period-by-period pricing, the compounded discounted payoff has the appropriate martingale property, whereas revaluation by an inverse kernel can have a nonzero covariance. State martingale-difference or conditional restrictions; pairwise uncorrelated errors alone are insufficient to eliminate all higher-order products. |
| R2-3.3: Martin/GvN | The benchmark can remove the systematic mismatch term under correct specification. It does not escape multiplicative-tail problems in general; see the clarification below. The referee's “64 months” is indeed a unit error: GvN uses quarters. |
| R2-4.1: NAVs | The review correctly distinguishes smoothing of return series from the level of a terminal NAV. Its single-root comparative static is conditional, not a proof of monotonicity for the pooled least-squares estimator. A common 0–100% haircut does not span all heterogeneous, state-dependent or negative appraisal biases. |
| R2-4.2: frequency/weights | The review's quarter-end finding is confirmed. Describe the actual construction; it does not verify the referee's separate claim about universal Preqin mid-quarter imputation. |
| R2-4.3: simulation realism | The simulated estimator observes aggregated fund flows, so the referee overstates deal observability. Common shocks, cash-flow timing, and terminal valuation nevertheless deserve a robustness design. Resolve the newly found timing mismatch first. |
| R2-4.4: follow-on investments | Reallocating exits to separate investment tranches is an accounting possibility, but arbitrary allocation does not preserve a conditional pricing restriction for every tranche. State the condition or work with an aggregate dynamic pricing equation. |
| Minor comments | Keep the useful corrections on ratio definitions, saturation, labels, and references, but check numbering against the exact submission version. The supplied combined PDF already numbers some propositions continuously where the referee-facing appendices differ. |

**The Martin qualification matters.** Review lines 737–749 describe a benchmark at the true beta as Martin's “non-generic case.” Martin defines that case using the full one-period risk-adjusted return approaching one. Canceling the systematic component leaves the idiosyncratic multiplier \(N\); if it remains nondegenerate, the full process still has the rare-event problem. The review recognizes this in its next bullet, but the conclusion that the pathology “does not reach” the estimator is too strong. In addition, its squared-loss limit need not equal the true beta, so exact systematic cancellation at the true beta does not establish cancellation at the estimated limit. The defensible point is that diversifiable idiosyncratic noise and a common realized-kernel shock are different statistical problems. See the definition and discussion in the local Martin paper, PDF page 6 / printed p. 350, and [Martin's original working paper](https://www.nber.org/system/files/working_papers/w16219/w16219.pdf).

The two-root result also needs existence and domain qualifications. The algebra

\[
E\left[\frac{1+r_f+\beta_{PE}F}{1+r_f+\beta F}\right]-1
=(\beta_{PE}-\beta)g(\beta)
\]

is correct where the expectations exist. The roots are \(\{\beta_{PE}\}\cup\{\beta:g(\beta)=0\}\), not necessarily exactly two distinct admissible roots. They may coincide or the second may lie outside the positive-return parameter domain. For Gaussian simple factor returns, a reciprocal affine return has a pole and the required population inverse moments generally do not exist for nonzero beta. The replica's finite grid and ordinary sample paths do not cure this population issue. Use a bounded-support factor design and a uniform positive-return margin for formal simulation claims, or clearly separate finite-sample numerical evidence from population assertions.

Finally, the review calls its Assumption A “DLP's own assumption.” DLP's published Assumption 1 specifies independence from market returns at all leads and lags and serial independence of project shocks; it explicitly distinguishes cross-project correlation assumptions from the assumptions used for inference (printed pp. 519–520). Conditioning additionally on the entire realized timing schedule, plus conditional independence across all deals, is a stronger sufficient formulation. It should be labeled as such. Endogenous stopping can preserve an unconditional stopped-martingale pricing equation while destroying conditional mean one after conditioning on the realized stopping time.

**8. Data and numerical evidence: confirmed findings and limits**

I reproduced the empirical aggregation issue using the stored global fund-level and vintage-level EW files. The stored vintage cash flow equals the date-specific mean over reporting funds, to floating-point accuracy. A fixed portfolio instead requires constant fund weights through time, with missing data distinguished from genuine zero flows.

| Segment/vintage | Funds | Reported-flow mean: positive/negative flow ratio | Constant-composition pool ratio |
|---|---:|---:|---:|
| Buyout 2000 | 42 | 2.760795 | 2.553055 |
| Buyout 2005 | 69 | 2.442632 | 1.975459 |
| Venture 2000 | 79 | 1.231860 | 0.886674 |
| Venture 2005 | 43 | 5.576438 | 2.181242 |

For Buyout 2000 the number of reporting funds is 8 at the first date, peaks at 41, and is 2 at the last date, matching the review. The offending operation is the date-specific `mean` in [prep_preqin.R](</Users/christausch/github/Fundwise_SDF/r_project/empirical/prep_preqin.R:229>). These numerical examples use the global file; the North America construction uses the same operation. A precise effect on the paper's selected North America estimates requires a reconstruction and rerun with the exact eligibility rules.

The direction of the estimated-beta change is not established merely by the cash-flow ratios. Both early contributions and late distributions are reweighted, and the estimator is nonlinear. Keep “could materially change the results”; remove any unconditional prediction of its direction.

All 49,626 nonzero BO/VC fund-level cash-flow observations I checked occur at quarter ends. The selected North America vintage portfolios have maximum date spans of 285 months for buyout and 354 months for venture, under both EW and FW. This confirms the saturation explanation.

The review's declaration that “all simulation and empirical numbers, and the headline empirical finding” stay unchanged should be withdrawn. A prose correction leaves the existing printed numbers as they are, but a construction or timing correction may change their validity and their values.

For the replica:

* I reran D1 and reproduced its reported gradients, including the seed-21 \(H=15\) numerical derivative of about +0.935 versus the conditional-expectation calculation of -1.710. The independent exact-enumeration check confirms the analytic formula despite this Monte Carlo noise.
* I also reran all seven D2/D3 configurations and both D5 configurations. They reproduce the displayed table entries at the reported precision. This verifies the replica tables, not their transfer to the R pipeline. I did not rerun D4's 500-block table or all of the earlier C-series checks.
* The D1 standard errors are computed as `sd/sqrt(V)`, even though adjacent vintages share overlapping market paths. For formal uncertainty statements, use independent replications or a dependence-aware estimate. A paired difference between the numerical and conditional calculations is also more informative than comparing two separate marginal standard errors.
* A 0.01 beta grid and a 0.1 percentage-point monthly alpha grid cannot establish exact unbiasedness or a perfectly flat profile. In D5, the positive-noise product estimate is actually beta 1.05 at \(H=15\), not exactly one. Say “close to the truth at this grid resolution.”
* Sample criterion values can be negative under the product correction; unrestricted optimization can exploit negative regions. Compact admissible parameter sets, identification, and a uniform LLN are material conditions, not optional implementation details.

**9. A replacement argument that is safe to build on**

The following is a compact, sufficient framework. It deliberately separates conditional pricing, pooling, identification and asymptotic regularity.

Let \(\mathcal C_v\) contain the common factor environment and a predetermined vintage valuation window. Suppose that, conditional on \(\mathcal C_v\), the normalized fund NPV functions \(X_{vi}(\theta)\) are independent with a common mean \(\mu_v(\theta)\) and variance \(\sigma_v^2(\theta)\). Assume \(\mu_v(\theta_0)=0\), a positive common weight \(A_{v,H}(\theta)\), finite moments, and an admissible compact parameter set. Then

\[
E[(A_{v,H}\bar X_v)^2\mid\mathcal C_v]
=A_{v,H}^2\left[\mu_v(\theta)^2+\frac{\sigma_v^2(\theta)}{n_v}\right].
\]

This gives the exact source of a possible squared-loss tilt: at the truth, the signal component has zero gradient while the noise component generally does not. The horizon changes both terms through \(A_{v,H}\). There is no general sign or monotonicity theorem for beta without further restrictions.

For the pair-product criterion,

\[
E[U_{v,H}(\theta)\mid\mathcal C_v]=A_{v,H}^2\mu_v(\theta)^2.
\]

If the expected signal criterion is uniquely minimized at \(\theta_0\), and the sample criterion converges uniformly to it as the vintage domain grows, its minimizer is consistent for every fixed admissible horizon. Independence and equal conditional means are sufficient conditions; if they are relaxed, the additional covariance and heterogeneity terms must be analyzed.

Under an infill asymptotic sequence, with suitable conditional LLNs and stable weights, the squared-loss noise term disappears. Its limiting zero set is then shared across positive horizon weights. This is the useful content of R1's conjecture. With bounded pooling, equality of the squared-loss target and the structural parameter is not guaranteed, including at inception.

In this formulation, the product criterion removes a population noise term, not finite-sample estimation error, weak identification, misspecification, or common shocks. At the population truth its horizon profile is flat under the maintained assumptions. Turning that implication into a formal test requires joint inference.

**Suggested replacement for the review's short answer:**

> Horizon averaging is an exact, parameter-dependent reweighting of each inception NPV. Under conditional mean-zero pricing given the relevant common environment, it does not introduce a mean pricing error. Nevertheless, with a bounded number of funds per portfolio, squared loss generally has a noise-dependent population target, including at the inception horizon. The target becomes structural under suitable infill and identification conditions. A cross-fund product criterion can remove this noise contribution under additional conditional independence and equal-mean assumptions; its applicability to heterogeneous fund schedules and the empirical construction requires separate verification. The direction and finite-sample value of a horizon change are model-dependent.

**10. Recommended order of work**

1. Resolve the simulation timing convention, the 0-versus-1-month labeling, and the padded horizon endpoint. Check zero-noise pricing before interpreting any bias plot.
2. Correct or justify the empirical portfolio construction, distinguishing reporting gaps from zero flows, then rerun the selected sample.
3. Retain the rescaling identity and the exact tilt formula; rewrite the consistency map using “not generally consistent” and explicit identification/regularity conditions.
4. Repair the conditional wedge lemma, the NFV interpretation, and the CLT-centering assumptions. Do not describe this as only a prose revision.
5. Reconcile Appendix A's pricing-set language with Appendix B's actual estimand. Lead the response to R1/R2 with the mathematical distinction rather than a terminology defense.
6. Treat the product criterion as a separate, testable extension. Start with predetermined common windows and conditional iid funds, then stress heterogeneous timing, common shocks, weights, and terminal NAVs before promising an empirical robustness column.
7. Reproduce the core comparison in the matched paper pipeline and report bias, dispersion, and a common loss separately. Only then decide how much of the stability claim survives.

The review's core insight is worth pursuing. The strongest revised paper would show precisely when reweighting preserves a valid pricing moment, why squared loss can nevertheless target the wrong parameter, and which additional assumptions make a correction work. The current review reaches those ideas, but its strongest conclusions need the qualifications and implementation checks above.
