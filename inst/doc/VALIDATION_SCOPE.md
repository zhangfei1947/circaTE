# Validation scope of the first circaTE package

This document separates software correctness from scientific qualification. The package makes the selected, existing research implementations reusable; it does not establish a new model or a universal false-discovery guarantee. See `MODELING_HISTORY.md` for the three development rounds and their source locations.

## Intended analysis

The first version targets matched RNA-seq and Ribo-seq counts for the same genes and biological samples, with known sampling times and culture/experimental batch labels. Separate assay library-size offsets are used. The primary model is a fixed-period first harmonic (H1): a joint test of the TE cosine and sine coefficients, with batch and linear time nuisance terms when requested and identifiable.

The development design consisted of 20 matched sample pairs, two continuing culture batches, and ten times over DD12–48. The primary period was 24 h; 22 h was a separate sensitivity analysis. Supporting a user-supplied period argument is a software capability, not validation of arbitrary periods or permission to select each gene's minimum p across periods without additional multiplicity control.

Estimated amplitude is the log2 half peak-to-trough amplitude of the fitted periodic component. Phase is its maximum modulo the specified period relative to time zero; it is not the maximum of a trend-plus-harmonic trajectory. Flat signals have no meaningful phase. Counts and offsets define relative TE, not absolute protein synthesis or global translational shifts.

## Status of methods

| Method family | What is provided | Qualification |
|---|---|---|
| Paired edgeR QL H1 | Joint count model with matched-sample effects and Ribo-specific harmonic contrast | Favored count-model development route; the raw test did not control error adequately in every autocorrelation/low-trough stress test |
| Precision-weighted AR GLS | Separate assay voom precisions combined for paired log-TE; culture-block AR1 plus nugget; ordinary or moderated F | Improves power relative to matched unweighted GLS in v3, but estimated precision/correlation uncertainty and low-count log bias remain |
| Matched unweighted AR GLS | Same H1 covariance-estimation design with all precision weights set to one | Controls an implementation comparison; it is not a generally calibrated fallback |
| Weighted independent GLS | Precision weighting with independent working errors | Diagnostic comparison; substantially inflated all-null event rate under the v3 mixture |
| Paired DESeq2 LRT | Raw count-model likelihood-ratio test | Raw extreme-tail inflation observed in all three rounds; reported for transparency, not as a universally reliable discovery engine |
| DESeq2 finite-F reference | Heavy-tailed heuristic reference for the LRT statistic with model degrees of freedom | Often conservative in the tested simulations; not an exact small-sample DESeq2 F test |

The older v2 unweighted AR code estimated covariance under an H2 mean design. The first-version matched H1 comparator should not be mislabeled as numerically identical to that historical implementation. Preserve method names, covariance-estimation design, and relevant diagnostics when comparing outputs.

H2, adaptive H1/H2 testing, and finite-grid covariance envelopes are historical experiments, not default first-version methods. H2 improved double-peaked signal detection in v2 but showed excessive all-null discoveries in its high-power form. It tests first and second harmonics together, so a pure 12 h component can yield a significant test at nominal period 24 h. Non-sinusoidal TE was not crossed with every deep-RNA-trough setting in v3. The unregularized single-gene joint NB-plus-AR prototype was unstable and excluded before the v3 holdout freeze; it is not a supported primary model.

## Empirical calibration is research evidence, not a real-data certification

The v2/v3 empirical maps were learned from independent simulated null data, frozen, and evaluated on independently seeded holdout simulations. They estimate tail probabilities for specified mixtures of expression, dispersion, timing, phase, and covariance. Pooled maps represent the training mixture; maximum-across-scenario maps are more conservative across the enumerated scenarios, not arbitrary future data-generating processes.

Do not apply bundled or historical simulation maps to real data and label the resulting q-values confirmed or calibrated at 5% FDR. The package's ordinary fitting functions must not silently perform that transfer. A reproducible analysis may report the raw working-model p-values and explicit sensitivity diagnostics, but must preserve the qualification that their tail behavior depends on the data and design.

Calibration requires its own analysis identity: period, model and version, design, normalization, filtering/testing universe, simulation mechanism, calibration seeds and strata, held-out validation, and a documented applicability assessment. A map trained at 24 h does not automatically apply to 22 h. A map trained on a 720-gene panel is not automatically justified for a full real transcriptome, particularly when weights and dispersions are estimated across genes. Mixing filtered and unfiltered universes also changes both the target population and BH behavior.

Empirical tail probabilities must not be zero or extrapolated beyond observed calibration support without separately validated tail modeling. In v3, the pooled minimum p was 1/64801; the cross-scenario maximum map's minimum was 1/7201. At 720 tests, the latter cannot produce a single isolated BH discovery at 0.05. Finite resolution is part of the power comparison, not a numerical defect to bypass with zero p-values.

Zero all-null discovery events in 270 simulated experiments does not prove a zero error rate; its binomial reference 95% upper bound is about 1.36%. Likewise, an overall mean FDP near 5% does not establish conditional validity in every low-expression subgroup. V3's calibrated weighted AR retained elevated descriptive FDP in some low-depth/deep-trough subgroups.

## Software checks and scientific checks have different purposes

Package tests should verify input validation, correct sample/gene alignment, identifiable designs, finite offsets, preservation of paired deletions, reproducible coefficients and contrasts, covariance algebra, invariance to jointly permuting samples, valid p-value ranges, declared failure handling, BH family definitions, and propagation of convergence/boundary flags. Reference comparisons against frozen research outputs can detect implementation drift. These checks do not establish FDR control for new biological data.

Invalid inputs or unidentified designs must produce actionable errors rather than silently switching models. All-zero/noninformative genes, failed fits, missing uncertainty estimates, and boundary covariance estimates require explicit flags and a documented output policy; they must not be silently dropped to improve the displayed results. Where a failed fit is conservatively assigned p=1, the failure status remains visible. A successful optimizer return does not mean the covariance parameters are well identified.

First-version package build/check and unit-test results should be recorded separately by the package release workflow. This document does not assert that a particular release passed merely because the research implementations did. The research provenance includes successful frozen-output audits, a conservatively handled DESeq2 nonconvergence in v3, and GLS optimizer warnings that were independently rechecked without changing the frozen outputs.

## What an applied result should retain

An interpretable result should preserve the tested gene universe and multiplicity family; assay feature definitions; sample matching, batch and time metadata; normalization and offsets; period and nuisance terms; method/version; p and q alongside effects; convergence and covariance-boundary flags; and low-count summaries for both assays. Report whether evidence survives the predetermined period and QC sensitivity analyses. Sensitivity outputs are not an uncorrected search for the smallest p-value.

Low-count observations should not be automatically removed on the basis of the observed TE trough. Such deletion selects phases and changes the estimand and precision. For NCU02302, the weighted amplitude persisted after downweighting low Ribo observations, whereas its q-value was sensitive to paired deletion, pseudocount, and period. This is evidence about influence and uncertainty, not confirmation or proof of artifact.

Filtering may be appropriate for a declared analysis domain, but conditional power among retained genes must be distinguished from recovery of all original true positives. A model's inability to detect a gene in a low-information trough does not show absence of TE rhythm. Arbitrary interpolation, pseudocount selection by significance, winsorization, or post hoc removal of trough phases was not validated as a remedy.

## Untested or incompletely tested extensions

The three rounds do not establish performance for arbitrary irregular sampling, missingness mechanisms, single-culture inference, many unequal batches, nonstationary or strong negative residual dependence, general cross-gene dependence, free-period search, arbitrary non-sinusoidal TE combined with deep RNA troughs, or robust estimation of uncertainty on phase/amplitude under those settings. Mathematical support for a design matrix or covariance calculation is not evidence that each such setting has acceptable operating characteristics.

Both biological studies contributed to development. New independent random seeds are independent simulated realizations, not external biological validation. Only 36 h of observation and two culture batches limit circadian interpretation. Without spike-ins, global shifts cannot be recovered; Ribo occupancy is not protein production rate. Cross-study comparisons, platform-confounded genotype differences, and eL31_KO without a same-study WT are outside the package's validated causal scope.

Further model optimization must use new development and validation partitions. The existing holdouts may be inspected for diagnosis, but once used to choose a correction they cannot serve again as untouched confirmation of that correction.
