# Codebase Review

This review is based on static inspection of the repository and cross-checking against `products/manuscript.qmd` and `products/supplement.qmd`. I could not execute the R pipeline in this shell because `R`/`Quarto` are not installed here, so I focused on reproducibility, code correctness, and consistency with the written methods/results.

## Findings

1.  **The current pipeline is not reproducible as checked in because several production functions still contain active `browser()` calls.**
    -   `R/05-tables.R:7,36,97`
    -   `R/04-measures-and-draws.R:405`
    -   `R/06-figures.R:205,815,1253`
    -   These functions are invoked by the main manuscript and supplement targets, including `tab_model_metrics` and `plot_measures_vs_1year_draws` (`_targets.R:399-405,519-533`), so even `targets::tar_make("reproduce_manuscript")` would stop in an interactive debugger rather than run end to end.
    -   This is the highest-priority issue because it blocks reproduction outright.
        -   DONE
2.  **Core-count handling is inconsistent and likely oversubscribes the machine, contradicting the reproducibility text.**
    -   `_targets.R:25-45` launches a local `crew` controller with `availableCores(omit = 2)`, which implies multiple targets may run concurrently.
    -   `R/02-model-fits.R:15-17,123-129` then gives each `brms::brm()` call `cores = detectCores() - 1`.
    -   In practice, that creates nested parallelism: multiple worker processes can each try to run several Stan chains in parallel. The supplement and README say the pipeline “leaves two cores free” (`products/supplement.qmd:45`, `README.md:19`), but the checked-in code does not guarantee that.
    -   I would either serialize model-fitting targets or set `cores = 1` inside `brm()` when using `crew`.
        -   DONE
3.  **The measure-processing code is out of sync with the measure-generation code and appears brittle or nonfunctional.**
    -   `R/03-model-predictions.R:11-13` still predicts out to day `1826` (5 years).
        -   DONE
    -   `R/04-measures-and-draws.R:237-239,277-279` expects `response_5year` and `time_response`.
        -   DONE
    -   But `R/04-measures-and-draws.R:166-185` only stores `peak`, `peak_time`, `growth_rate`, `peak_response`, `baseline`, `decay_rate`, `titer_1year`, and `response_1year`.
    -   That leaves the downstream processing helper referencing objects that are never created anywhere in the repo. At minimum this is dead/stale code; depending on how `tibble()` handles those `NULL`s at runtime, it may also break `process_*_group_measures` targets (`_targets.R:315-329` and the analogous supplement targets).
4.  **The repository does not contain the code needed to reconstruct key preprocessing decisions described in the paper.**
    -   `R/01-read-data.R:5-16` just loads preprocessed `.rds` files.
    -   The manuscript and supplement describe nontrivial preprocessing decisions: infection-based cohort restriction (`products/manuscript.qmd:257-261`), censoring rules (`products/manuscript.qmd:279-283`, `products/supplement.qmd:189-191`), assay-specific bounds, and timepoint harmonization.
    -   None of that provenance is reconstructable from source data in the current repo. Users can reproduce the modeling stage from processed inputs, but not the upstream data assembly that underpins the scientific claims.
        -   DONE: not relevant
5.  **The absolute-measures supplement figures do not match the checked-in plotting code or their own captions.**
    -   The captions say the plots include “Smaller open points” for individual-level posterior estimates and use a “pseudo-log scale” (`products/supplement.qmd:1047-1058`).
    -   But the pipeline never computes individual-level measure targets for these figures; it only feeds `process_*_group_measures` into `func_plot_absolute_measures()` (`_targets.R:690-711`).
    -   `R/06-figures.R:1134-1228` plots only one `geom_pointrange()` layer from the group-level data and uses ordinary `scale_y_log10()` / `scale_y_continuous()`, not a pseudo-log transformation.
    -   So the checked-in code is inconsistent with the documented figure content.
        -   DONE
6.  **Several manuscript/supplement descriptions of implemented measures do not match the code.**
    -   `R/04-measures-and-draws.R:109-114,159-164` defines `peak_response` and `response_1year` as differences on the log scale.
    -   The manuscript repeatedly calls those quantities “fold change” or “fold response” (`products/manuscript.qmd:149,160,183,283`), which is not what the code computes unless those values are later exponentiated.
        -   DONE: got rid of fold-change
    -   `R/02-model-fits.R:31,48` labels `k` as “half-peak time” in code comments, while the manuscript/supplement describe `K` as the time of peak (`products/manuscript.qmd:275`, `products/supplement.qmd:72,117,640`).
        -   DONE
    -   `R/04-measures-and-draws.R:56-60` computes growth based on reaching half of peak titer on the original scale, while the supplement text describes the natural log-transformed titer reaching 50% of its peak value (`products/supplement.qmd:558`), which is a different quantity.
        -   DONE
7.  **The figure code applies additional filtering that is not fully disclosed in the captions.**
    -   `R/04-measures-and-draws.R:380-395` already trims each measure to its middle 90% and samples 500 draws.
    -   `R/06-figures.R:293-298` then applies another 95% truncation for growth-rate panels.
    -   `R/06-figures.R:300-311` additionally removes the bottom 10% of the `response_1year`/`titer_1year` axis variable for baseline-relative plots.
    -   The main-text caption for Figure 6 says the smaller points are “sampled after removing the upper and lower 5% of draws for each measure” (`products/manuscript.qmd:204`), which understates what the plotting code is actually doing.
        -   DONE
8.  **The main trajectory plotting code hardcodes an infection-specific x-axis label even for vaccination figures.**
    -   `R/06-figures.R:74-75` labels the x-axis “Days post-inoculation.”
    -   The same function is used for both NI and NV trajectories (`_targets.R:467-486`).
    -   For the vaccination study this should be “Days post-vaccination” or a neutral label such as “Days since exposure/intervention.”
        -   DONE
9.  **The supplement decay-day combination plot appears to discard the legends that the caption says should be present.**
    -   `R/06-figures.R:1236-1240` removes the legend from each sub-plot before combining them.
    -   But the caption for `fig-peak-decay-days-1year` says color represents the predicted one-year response “using antibody-specific color scales shown in panel order” (`products/supplement.qmd:1024`).
    -   Unless those scales are added elsewhere after the current code path, the checked-in helper does not preserve them.
        -   DONE: scales are fine
10. **The schematic figure is not faithful to the implemented measure definitions.**
    -   `R/06-figures.R:733-741` uses `half_peak_day <- peak_day * 0.5` to locate the “growth rate” segment, which is not the same as the actual code definition in `R/04-measures-and-draws.R:54-90`.
    -   `R/06-figures.R:738-793` labels decay using days 85 and 95, whereas the measure code defines decay over one-day intervals at 60, 90, 180, and 365 days (`R/04-measures-and-draws.R:116-147`).
    -   The inline comment at `R/06-figures.R:735` also says growth rate is “technically” using 70% peak titer, which conflicts with the implemented 50% threshold.
    -   Since the supplement presents this as a schematic of the derived measures (`products/supplement.qmd:543-545`), it should reflect the actual definitions more closely.
        -   DONE
11. **One plotting helper has a latent bug because it hardcodes `decay_rate_D90` despite accepting a variable `decay_day`.**
    -   In `R/06-figures.R:210-212`, the function constructs `decay_name` dynamically.
    -   But `R/06-figures.R:243` still pivots `decay_rate_D90` explicitly instead of using `decay_name`.
    -   Current `_targets.R` calls this helper only with `decay_day = 90`, so the bug is latent right now, but the function signature suggests it is intended to be reusable and it is not actually generic.
        -   DONE
12. **The checked-in documentation has multiple reproducibility mistakes that should be fixed alongside the code.**
    -   `README.md:14` says the pipeline file is `targets.R`, but the actual file is `_targets.R`.
        -   DONE
    -   `README.md:23` and `products/supplement.qmd:49` use the wrong `.Rproj` filename.
        -   DONE
    -   `README.md:19,23` and `products/supplement.qmd:45,49` use `Rstudio` instead of `RStudio`.
        -   DONE
    -   `products/supplement.qmd:187` cites Stan/HMC with the `targets` citation rather than a Stan/cmdstanr citation.
        -   DONE

## Open questions / assumptions

-   I could not inspect the contents of `data/dat_ni.rds` and `data/dat_nv.rds` directly in this shell, so I could not verify whether the NI infection subset is based purely on virologic confirmation or partly on the antibody-based criterion discussed in the manuscript.
-   Some generated outputs already exist under `results/`, so a subset of these issues may have been introduced after those artifacts were built. The review above is about the **current checked-in source code**, not necessarily the exact source state that produced the existing figures/tables.

## Overall assessment

The repository has a sensible high-level structure and the analysis logic is readable, but the current checked-in state is not submission-grade from a reproducibility standpoint. The immediate priorities are:

1.  Remove all `browser()` calls.
2.  Resolve the stale `response_5year` / `time_response` scaffolding.
3.  Fix the nested parallelism/core-allocation strategy.
4.  Bring the manuscript/supplement text and figure captions into alignment with the implemented measures and plotting code.
5.  Decide whether the repo is intended to reproduce the analysis from processed inputs only or from raw study data, and document that boundary explicitly.
