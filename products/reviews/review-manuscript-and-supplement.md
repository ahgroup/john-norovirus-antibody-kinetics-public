# Review of `manuscript.qmd` and `supplement.qmd`

This review is based on the Quarto source files and a spot-check of the rendered `.docx` outputs. I could not rerender locally because this shell does not have `R`/`Quarto`, so line references below point to the `.qmd` sources.

## Major findings

1.  **The manuscript and supplement repeatedly misinterpret the model parameters `P` and `K`.**
    -   `manuscript.qmd:267-275` and `supplement.qmd:62-76,117,640` describe `P` as the "peak" antibody level and `K` as the "time at peak"/"time at which the peak response occurs."
    -   For the written model, the logistic term `P / (1 + exp(-G(t-K)))` has asymptote `P` and midpoint `K`; once multiplied by the decay term, the actual peak time and peak magnitude are determined jointly by the induction and decay terms and are generally **not** equal to `K` and `P`.
    -   This matters scientifically because the current wording gives readers the impression that the fitted parameters map directly to biologically interpretable peak timing and magnitude, when those quantities are actually derived later from the full predicted trajectory.
    -   Recommendation: rename `P` and `K` throughout to something like "amplitude/asymptote" and "inflection/half-rise time" (or similar), then state clearly that peak magnitude and peak time were computed numerically from posterior-predicted trajectories.
        -   DONE
2.  **The response outcomes are described as fold changes, but the methods define them as differences on the log scale.**
    -   `manuscript.qmd:149,160,183,193-214,229,283` and `supplement.qmd:554-567` define or discuss "peak response" and "one-year response."
    -   In the methods, these are explicitly defined as `log(peak) - log(baseline)` and `log(titer at 1 year) - log(baseline)` (`manuscript.qmd:283`), which are log-fold changes, not fold changes.
    -   Calling these "fold change" or "2-fold higher" is misleading unless the quantities were exponentiated before plotting/reporting.
    -   Recommendation: either back-transform and report actual fold changes everywhere, or relabel consistently as "log-fold change," "log response relative to baseline," or "difference in log titer."
        -   DONE: i think its fine as is without mentioning fold-change
3.  **The NI analysis may be circular because cohort inclusion uses an antibody-based infection definition.**
    -   `manuscript.qmd:97,257` state that only infected NI participants were included, and infection was defined partly by a `>=4`-fold rise in serum antibody titers from baseline to day 28.
    -   The same antibody outcomes are then modeled as the main response variables.
    -   If any included NI participants qualified as "infected" based only on serologic rise rather than virologic evidence, the study design could partially condition on the outcome being analyzed, which can bias dose-response comparisons and inflate apparent antibody induction.
    -   Recommendation: state how many included NI participants were virologically confirmed versus serology-only, and ideally add or cite a sensitivity analysis restricted to virologically confirmed infections.
        -   DONE: added a sentence in manuscript

## Moderate findings

4.  **The statistical notation around the modeled scale is inconsistent.**
    -   `supplement.qmd:94-108` says the model is specified for `ln(mu_{i,t})`, but the likelihood is then written as `y_{i,t} ~ Student(v=4, mu_{i,t}, sigma)`.
    -   `manuscript.qmd:277-283` likewise says titers were converted to natural log units, but the prose around the expected trajectory could be clearer about whether `mu` denotes the expected log titer or the expected titer before logging.
    -   Recommendation: use one notation consistently, e.g. define `eta_{i,t} = E[log y_{i,t}]` and put the Student-`t` likelihood on the log scale.
        -   DONE
5.  **The random-effects description in the main manuscript is inaccurate.**
    -   `manuscript.qmd:279` says the models included "subject-specific random intercepts."
    -   `supplement.qmd:115-117,126-183` makes clear that there are subject-level random effects on `p`, `g`, `k`, and `d`, not just a single intercept term.
    -   Recommendation: revise the manuscript wording to "subject-specific random effects on the latent kinetic parameters" or equivalent.
        -   DONE
6.  **The description of the uncertainty intervals appears inconsistent with the stated workflow.**
    -   Many captions describe "95% MAD-based credible intervals" (`manuscript.qmd:117,128,160,183,204`; `supplement.qmd:534,810,819,919,930,976,987`).
    -   The supplement also says `posterior_summary(..., robust = TRUE)` was used (`supplement.qmd:534`). That function returns a robust point estimate, but the 95% interval is typically the 2.5% and 97.5% posterior quantiles, not a MAD-derived interval.
    -   Recommendation: verify what was actually plotted and reported, then relabel either as "95% credible intervals" or explicitly as quantile intervals if that is what was used.
        -   DONE: ETCIs
7.  **Some discussion statements overreach relative to the data strength.**
    -   `manuscript.qmd:239-243` moves from small-group, mostly overlapping-CI results to mechanistic claims about high-affinity B-cell recruitment, downstream regulatory constraints, and infection engaging broader immune pathways than vaccination.
    -   These are plausible hypotheses, but they read stronger than the data support, especially with the small NI subgroup sizes and the cross-study differences in design.
    -   Recommendation: soften these sentences by framing them as hypotheses or possible explanations rather than interpretations the data directly establish.
8.  **The supplement’s section hierarchy produces awkward numbering in the rendered document.**
    -   `supplement.qmd:119` and `151` use `####` headings directly under a `##` section.
    -   In the rendered `.docx`, these appear as `2.2.0.1` and `2.2.0.2`, which looks unpolished.
    -   Recommendation: insert a `###` heading level or demote these headings to bold lead-ins instead of numbered section headers.

## Minor and editorial issues

9.  **Typos and wording issues in the manuscript.**
    -   `manuscript.qmd:97`: "Material and Methods" should match the actual section title "Materials and Methods."
        -   DONE
    -   `manuscript.qmd:102`: "all timepoints" -\> "all time points"; "150mcg" -\> "150 mcg."
        -   DONE
    -   `manuscript.qmd:136`: "expect" -\> "except."
        -   DONE
    -   `manuscript.qmd:140`: "values near 1 and indicate" -\> remove "and."
        -   DONE
    -   `manuscript.qmd:195,197,199,212`: "SM 3.4" should be "SM S3.4" for consistency with the rest of the paper.
        -   DONE
    -   `manuscript.qmd:172,176`: shorthand such as `Pr(150 < others)` is too informal/undefined for the main text; specify the actual pairwise comparisons.
        -   DONE
    -   `manuscript.qmd:285`: `cmdrstan` should be `cmdstanr` in both places.
        -   DONE
10. **Typos and reproducibility issues in the supplement.**
    -   `supplement.qmd:45,49`: "Rstudio" -\> "RStudio."
        -   DONE
    -   `supplement.qmd:49`: the `.Rproj` filename is wrong. The repository file is `john-norovirus-antibody-kinetics-public.Rproj`, not `john-norovirus-ab-kinetics-public.Rproj`.
        -   DONE
    -   `supplement.qmd:187`: the citation attached to Stan/HMC is the `targets` citation, which is incorrect.
        -   DONE
    -   `supplement.qmd:545`: "peak resonse" -\> "peak response."
        -   DONE
    -   `supplement.qmd:706`: "by dose, and antibody" is awkward; remove the comma.
        -   DONE
11. **A few phrasing choices could be tightened for scientific style.**
    -   `manuscript.qmd:91`: "A key feature of any therapeutic is the dose" is awkward in a study about infection challenge and vaccination; "exposure characteristic" or "intervention design feature" would be more precise.
        -   DONE: updated
    -   `manuscript.qmd:65`: "Co-corresponding authors" is understandable but "Co-corresponding authors:" or "Corresponding authors:" reads more cleanly.
        -   DONE
    -   `supplement.qmd:540`: "We derive several antibody kinetics measures" would read better as "We derived" or "We computed."
        -   DONE

## Overall assessment

This is a strong and interesting project with a worthwhile data set, a sensible hierarchical modeling strategy, and a generally coherent main narrative. The main issues are not about the overall idea; they are about **precision of interpretation and reporting**. Before submission, I would strongly prioritize fixing the parameter interpretation (`P`/`K`), clarifying the scale of the response outcomes (log-change versus fold change), and addressing the possible circularity in the NI inclusion rule. Once those are cleaned up, the remaining issues are mostly editorial and presentation-level.
