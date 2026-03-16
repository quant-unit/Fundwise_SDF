# Conservative implementation plan for applying `refinement_plan.md` to `sdf_paper.tex`

## Guiding principle

Implement the reframing as a **light-touch editorial pass**, not as a substantive rewrite of the paper. The current draft already contains much of the desired identification/stability framing, especially in the abstract, introduction, methodology, empirical interpretation, and appendix. The safest plan is therefore to preserve the existing structure, equations, propositions, tables, and figures, and to focus on tightening section framing and removing the remaining overstatements.

Any new sentence should be traceable to one of the following already in the paper:
- the abstract and introduction on horizon choice, target shift, and stability ([sdf_paper.tex](/Users/christausch/Library/Mobile Documents/com~apple~CloudDocs/Projects/Fundwise_SDF/2025 sdf paper/tex/sdf_paper.tex#L91), [sdf_paper.tex](/Users/christausch/Library/Mobile Documents/com~apple~CloudDocs/Projects/Fundwise_SDF/2025 sdf paper/tex/sdf_paper.tex#L114));
- the methodology on structural vs. horizon-specific targets ([sdf_paper.tex](/Users/christausch/Library/Mobile Documents/com~apple~CloudDocs/Projects/Fundwise_SDF/2025 sdf paper/tex/sdf_paper.tex#L546), [sdf_paper.tex](/Users/christausch/Library/Mobile Documents/com~apple~CloudDocs/Projects/Fundwise_SDF/2025 sdf paper/tex/sdf_paper.tex#L655));
- the simulation results on finite-sample behavior and horizon sensitivity ([sdf_paper.tex](/Users/christausch/Library/Mobile Documents/com~apple~CloudDocs/Projects/Fundwise_SDF/2025 sdf paper/tex/sdf_paper.tex#L1457), [sdf_paper.tex](/Users/christausch/Library/Mobile Documents/com~apple~CloudDocs/Projects/Fundwise_SDF/2025 sdf paper/tex/sdf_paper.tex#L1608));
- the empirical interpretation on weak identification of richer models and stronger stability of parsimonious market specifications ([sdf_paper.tex](/Users/christausch/Library/Mobile Documents/com~apple~CloudDocs/Projects/Fundwise_SDF/2025 sdf paper/tex/sdf_paper.tex#L1995));
- the appendix on target shift and finite-sample stabilization under horizon averaging ([sdf_paper.tex](/Users/christausch/Library/Mobile Documents/com~apple~CloudDocs/Projects/Fundwise_SDF/2025 sdf paper/tex/sdf_paper.tex#L2088)).

## What should not change

To keep the revision robust and conservative, do **not** change:
- equations, propositions, corollaries, assumptions, proofs, or appendix derivations;
- simulation design, empirical sample construction, reported estimates, or figure/table content;
- the core section structure, unless a later pass shows a strong need for local paragraph reordering.

This means the revision should be mostly sentence-level and paragraph-level, with no changes to the paper's underlying results.

## Recommended implementation order

### Pass 1: keep the existing framing that already works

Preserve the current identification/stability language in:
- the abstract ([sdf_paper.tex](/Users/christausch/Library/Mobile Documents/com~apple~CloudDocs/Projects/Fundwise_SDF/2025 sdf paper/tex/sdf_paper.tex#L91));
- the introduction, especially the discussion of horizon choice, target shift, stability, and weak identification ([sdf_paper.tex](/Users/christausch/Library/Mobile Documents/com~apple~CloudDocs/Projects/Fundwise_SDF/2025 sdf paper/tex/sdf_paper.tex#L114), [sdf_paper.tex](/Users/christausch/Library/Mobile Documents/com~apple~CloudDocs/Projects/Fundwise_SDF/2025 sdf paper/tex/sdf_paper.tex#L150));
- the methodology and asymptotic framing around \(H=0\) versus \(H>0\) ([sdf_paper.tex](/Users/christausch/Library/Mobile Documents/com~apple~CloudDocs/Projects/Fundwise_SDF/2025 sdf paper/tex/sdf_paper.tex#L546), [sdf_paper.tex](/Users/christausch/Library/Mobile Documents/com~apple~CloudDocs/Projects/Fundwise_SDF/2025 sdf paper/tex/sdf_paper.tex#L682));
- the empirical interpretation subsection ([sdf_paper.tex](/Users/christausch/Library/Mobile Documents/com~apple~CloudDocs/Projects/Fundwise_SDF/2025 sdf paper/tex/sdf_paper.tex#L1992)).

These parts already align closely with the reframing memo and should serve as anchor text for later local edits elsewhere.

### Pass 2: tighten the first three pages without changing the paper's substance

Goal: make the introduction read more cleanly as an identification/stability paper while preserving the current content.

Safe edits:
- move the identification question slightly earlier within the introduction, but keep the existing private-equity motivation paragraph;
- tighten the contribution wording so that contribution 1 is about how horizon choice changes the estimand, contribution 2 is about stability and dependence-robust inference, and contribution 3 is about what the empirical panel supports stably;
- keep the current caveats that the paper does not show alpha is zero and does not claim a definitive population factor structure.

Because the introduction already says most of this, this pass should mainly compress and sharpen, not add new ideas.

### Pass 3: reframe the simulation section as finite-sample identification analysis

Goal: keep all simulation evidence, but make the interpretation more consistently about recoverability and weak identification.

Priority locations:
- Horizon discussion and simulation conclusion ([sdf_paper.tex](/Users/christausch/Library/Mobile Documents/com~apple~CloudDocs/Projects/Fundwise_SDF/2025 sdf paper/tex/sdf_paper.tex#L1608), [sdf_paper.tex](/Users/christausch/Library/Mobile Documents/com~apple~CloudDocs/Projects/Fundwise_SDF/2025 sdf paper/tex/sdf_paper.tex#L1624)).

Edits to make:
- replace stronger language such as "most critical", "drastically improves", and "all other measures combined" with sample- and design-specific language;
- describe horizon averaging as a regularization or stabilization device that improves finite-sample recoverability while shifting the target, consistent with the methodology and appendix;
- make sure statements about the 12.5--15 year range remain explicitly tied to the reported simulation designs rather than presented as universal.

This pass should not alter any numerical results or figure references.

### Pass 4: soften the strongest empirical overclaims and align the whole empirical section with the existing interpretation subsection

Goal: bring the tone of the empirical results section into line with the already well-calibrated interpretation subsection.

Priority locations:
- NAV sensitivity language ([sdf_paper.tex](/Users/christausch/Library/Mobile Documents/com~apple~CloudDocs/Projects/Fundwise_SDF/2025 sdf paper/tex/sdf_paper.tex#L1383), [sdf_paper.tex](/Users/christausch/Library/Mobile Documents/com~apple~CloudDocs/Projects/Fundwise_SDF/2025 sdf paper/tex/sdf_paper.tex#L1417));
- baseline single-factor discussion ([sdf_paper.tex](/Users/christausch/Library/Mobile Documents/com~apple~CloudDocs/Projects/Fundwise_SDF/2025 sdf paper/tex/sdf_paper.tex#L1694));
- two-factor \(q^5\) summary ([sdf_paper.tex](/Users/christausch/Library/Mobile Documents/com~apple~CloudDocs/Projects/Fundwise_SDF/2025 sdf paper/tex/sdf_paper.tex#L1970));
- comparison-to-\cite{DLP12} empirical subsection ([sdf_paper.tex](/Users/christausch/Library/Mobile Documents/com~apple~CloudDocs/Projects/Fundwise_SDF/2025 sdf paper/tex/sdf_paper.tex#L2689)).

Edits to make:
- replace phrases such as "strong empirical backing", "immediately collapses", "driven by finite-sample noise rather than genuine systematic risk exposures", and "severely under-identified" with narrower sample-specific formulations;
- reuse the paper's own safer language: "comparatively stable", "weakly identified", "sensitive to horizon choice", "diagnostic evidence", "does not imply absence in population";
- where possible, shift conclusions from "best model" language to "most stably estimated specification in the current sample/design" language.

The empirical interpretation subsection already provides good target wording and should be used as the benchmark for this pass.

### Pass 5: rewrite the conclusion around what the data can and cannot identify

Goal: make the conclusion consistent with the introduction and empirical interpretation, while reducing statements that sound universal or too strong.

Priority location:
- the opening and middle of the conclusion ([sdf_paper.tex](/Users/christausch/Library/Mobile Documents/com~apple~CloudDocs/Projects/Fundwise_SDF/2025 sdf paper/tex/sdf_paper.tex#L2038)).

Edits to make:
- foreground the central tradeoff already proved in the paper: inception-only pricing preserves the structural restriction, while horizon averaging can improve finite-sample stability but shifts the target;
- keep the sample-specific empirical takeaway that current panels are more informative about broad market exposure than about alpha or richer decompositions;
- preserve the existing caution that the results do not imply alpha is zero and do not prove that private equity is literally one-factor in population.

The current conclusion already contains much of this material, so this should be a tightening pass, not a full replacement.

## Optional items that should be deferred unless needed

To keep the revision conservative, postpone the following unless the core wording pass leaves a clear gap:
- changing the title;
- adding new tables or figures, including an "identification taxonomy" table;
- restructuring sections or moving large blocks of text;
- expanding criticism of prior papers beyond what the current draft already documents.

These may be useful later, but they are not required for the current reframing and they increase the risk of overshooting the paper's evidence.

## Editing rules to enforce during implementation

- Prefer replacing or trimming existing sentences over adding new paragraphs.
- If a sentence sounds broader than the evidence, narrow it with one of these anchors already used in the draft: "in our sample", "under our semiparametric design", "in the available cash-flow panels", or "in the reported simulations".
- When discussing \(H>0\), always keep the target-shift caveat nearby.
- When discussing single-factor results, avoid wording that implies the population factor structure has been discovered.
- When discussing alpha or additional factors, describe instability as weak identification in the current data rather than proof of economic irrelevance.

## Validation checklist after the edit pass

After each revision block:
- compile the LaTeX file and confirm there are no broken references or formatting regressions;
- inspect the diff to verify that changes are rhetorical rather than technical;
- confirm that no new claim exceeds what is explicitly supported by theorems, simulations, tables, or figures;
- confirm that the paper still reads consistently from abstract through conclusion.

## Practical recommendation

Execute the revision in this order:
1. abstract + introduction;
2. simulation framing;
3. empirical results tone;
4. conclusion;
5. final cleanup pass for repeated phrases and consistency.

That sequence is safest because the abstract and introduction establish the framing standard, and the later sections can then be aligned to that standard with minimal disruption.
