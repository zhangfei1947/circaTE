# circaTE

**Fixed-period translational-efficiency modeling for paired RNA-seq and Ribo-seq.**

circaTE is an R package that tests harmonic changes in relative translational efficiency (TE), accounts for biological sample pairing, and provides diagnostics for low-count troughs and temporal dependence. It includes a command-line interface for workstation and cluster use.

Version **0.1.0** is a research development prerelease. The software is licensed under the **MIT License** and maintained by **Phil Zhang** (<philzhang@tamu.edu>). Model p-values remain conditional and exploratory: software correctness does not establish universal false discovery rate (FDR) control or prove endogenous circadian rhythmicity.

## Installation

R >= 4.4 is required. Install the dependencies first:

```r
install.packages(c("BiocManager", "statmod", "jsonlite", "digest"))
BiocManager::install(c("edgeR", "limma"))

# Optional, for the DESeq2 comparison methods:
BiocManager::install("DESeq2")
```

Download `circaTE_0.1.0.tar.gz` from the [v0.1.0 release](https://github.com/zhangfei1947/circaTE/releases/tag/v0.1.0), then run:

```sh
R CMD INSTALL circaTE_0.1.0.tar.gz
```

Alternatively, install from the public source repository:

```sh
git clone https://github.com/zhangfei1947/circaTE.git
R CMD INSTALL circaTE
```

The package has been validated with R 4.4.2, edgeR 4.4.2, limma 3.62.2, DESeq2 1.46.0, and statmod 1.5.0. Other versions may change numerical results; outputs record the versions used. circaTE is not currently distributed through CRAN or Bioconductor.

## Quick start

```r
library(circaTE)

data <- example_te_data()
fit <- fit_te(data, method = "edger", period = 24)
fit
head(fit$results)

qc <- te_qc(data)
head(qc$genes)
plot_te_gene(fit, "gene001")
write_te_results(fit, "example_results")
```

The example contains 120 synthetic genes and demonstrates the API. Its discovery count is not a benchmark of power or FDR. The output directory must be new or empty.

## Input requirements

Each analysis represents one biological group. Provide:

- Two raw, nonnegative integer count matrices, `rna` and `ribo`, with identical genes in rows and paired biological sample IDs in columns, in the same order.
- A metadata data frame containing `sample_id`, numeric `time_hours`, and `batch`, in count-column order. A batch identifies a culture followed over time; duplicate batch/time combinations are rejected.
- Consistent, biologically appropriate feature definitions across assays, such as RNA CDS counts and Ribo CDS/P-site counts.

If RNA and Ribo libraries have different identifiers, explicitly map them to biological sample IDs using a reliable pairing manifest. circaTE does not infer pairing from column positions, silently intersect genes, impute observations, or remove low-count observations.

```r
data <- te_data(rna, ribo, metadata)

# If normalization has already been calculated:
data <- te_data(
  rna, ribo, metadata,
  offsets = list(
    rna = log_rna_effective_library_sizes,
    ribo = log_ribo_effective_library_sizes
  )
)

# TSV count files: gene_id first, then one column per biological sample.
data <- read_te_data("RNA.tsv", "Ribo.tsv", "samples.tsv")
```

By default, TMM normalization is estimated separately for RNA and Ribo. There is no automatic abundance filter, and all attempted genes remain in the BH family. Genes with zero total counts in either assay receive `p = 1`, missing effect estimates, and an uninformative status.

If count matrices have already been filtered, supply offsets calculated from an appropriate full-gene normalization universe. Do not re-estimate full-library normalization from a small candidate list. Offsets are **natural logarithms of effective library sizes**, not raw library sizes.

## Available methods

| `method` | Purpose and qualification |
|---|---|
| `edger` | Default paired joint edgeR v4 quasi-likelihood model. Shared RNA variation is absorbed by sample-pair effects; serial TE dependence is not explicitly modeled. |
| `weighted_ar` | Observation-precision weighted, culture-block AR1 plus nugget GLS with variance moderation. Low-count log bias and estimated-covariance uncertainty remain. |
| `unweighted_ar` | Matched H1 covariance-estimation comparator with unit observation weights. It is not the historical v2 H2 covariance-estimation implementation. |
| `weighted_independent` | Precision-weighted, moderated GLS with independent working errors. |
| `weighted_ar_ordinary` | Weighted AR GLS without variance moderation. |
| `deseq2` | Paired DESeq2 LRT, retained as a research comparator because stress tests identified excessive extreme-tail false positives. |
| `deseq2_finite` | A conservative, heuristic finite-degrees-of-freedom F reference for the LRT; not an exact small-sample test. |

The default includes TE batch effects and a linear trend. Specify the period in advance. Development used 24 hours as the primary period and 22 hours as a separate sensitivity analysis; the package does not select the smallest p-value across periods or methods.

H2/adaptive models are not part of the v0.1.0 analysis API. Their improved detection of some non-sinusoidal signals came with unresolved false-positive behavior.

### Design requirements

AR models require at least two independent culture batches. Weighted models additionally need residual replication for the categorical-time plus batch voom precision design. Completely staggered sampling times across cultures may leave no residual degrees of freedom for estimating precision.

The package checks the H1 design actually used for fitting, including identifiable six-hour sampling schedules. A relative singular-value threshold of `1e-8` rejects numerically aliased designs, such as attempting a full 24-hour harmonic from 12-hourly observations. Negative AR correlation is supported only on a regular time lattice; irregular timing restricts the AR parameter to nonnegative values. Mathematical support for a design does not imply that its statistical operating characteristics have been validated.

## Results and diagnostics

`fit$results` includes `gene_id`, `p`, `q`, `cos`, `sin`, `amplitude_log2`, `phase_hours`, `informative`, `fit_failed`, and `status`.

- `amplitude_log2` is the harmonic half peak-to-trough range of log2 TE.
- `phase_hours` is the harmonic peak modulo the specified period, relative to time zero, excluding the linear trend.
- `q` uses BH correction over all attempted genes in one fit. There is no automatic combined correction across biological groups, methods, or periods.
- Failed gene tests remain in the family at `p = 1`. Covariance convergence and boundary diagnostics are retained; the software does not silently switch models after a failure.

`te_qc()` reports assay-specific minimum counts, numbers of low-count observations, and the longest within-culture run of low counts. A gap larger than that culture's minimum sampling interval breaks a run. This is an operational count of observations, not an estimate of trough duration; interpret irregular sampling using the actual times.

`plot_te_gene()` shows descriptive RNA, Ribo, and log2 TE trajectories. Its `pseudocount` argument changes the displayed ratio only; it does not refit the GLS model, which uses a fixed pseudocount of 0.5.

`write_te_results()` exports results and QC TSV files, the complete fit RDS, a JSON manifest containing configuration, versions and diagnostics, and session information. It refuses to overwrite a nonempty directory.

## Command-line and cluster use

Locate the installed script:

```sh
Rscript -e 'cat(system.file("scripts", "circaTE.R", package="circaTE"))'
```

Use the printed path in place of `/path/to/circaTE.R`:

```sh
Rscript /path/to/circaTE.R \
  --rna RNA.tsv --ribo Ribo.tsv --metadata samples.tsv \
  --method edger --period 24 --out results_24h
```

Run `Rscript /path/to/circaTE.R --help` for available options. The CLI treats inputs as observed data and provides no switch to bypass the simulation-only calibration restriction. On an HPC system, run full analyses through the site's scheduler on compute nodes.

## Experimental simulation calibration

`build_te_calibrator()` constructs empirical null maps from independent, entirely null simulation fits. `calibrate_te()` accepts simulation fits only and checks matching method, period, trend, timing design, gene count, normalization, software versions, and simulation population context. Dataset IDs and content fingerprints prevent accidental reuse of training observations as evaluation data.

```r
# Supply independently generated, entirely null simulation fits.
calibrator <- build_te_calibrator(
  null_fits, scenario = scenario_labels, all_null = TRUE
)
calibrated <- calibrate_te(heldout_fit, calibrator, mode = "max")
```

These checks cannot prove that user-supplied null labels or independence claims are true. Correct simulation design and independent held-out evaluation remain necessary.

Historical v2/v3 calibration maps are not bundled or automatically transferred to observed data. Empirical calibration uses inclusive ties and an add-one correction, with no zero p-values or tail extrapolation. Finite p-value resolution can limit BH discoveries. Acceptable average performance in a simulation mixture does not guarantee conditional FDR control within low-expression subgroups or in new biological data.

## Evidence and validation

The [modeling history](inst/doc/MODELING_HISTORY.md) documents three development rounds and changes in the evidence; that detailed research summary is currently in Chinese. The [validation scope](inst/doc/VALIDATION_SCOPE.md) and [software validation summary](VALIDATION.md) are in English.

The simulation rounds used different difficulty, filtering domains, and signal shapes. Their power estimates must not be compared directly as measures of improvement. In the v3 stress distribution, scenario-max calibrated edgeR provided a useful power/error tradeoff, but this is not a universal real-data calibration guarantee. NCU02302 remains an exploratory candidate sensitive to analysis choices.

Validation includes input alignment, zero-information handling, design checks, historical numerical regression, sample-order invariance, irregular-time handling, calibration leakage prevention, and inclusive empirical probabilities. Two checkout-specific regression tests explicitly skip when the non-distributed historical research files are absent; core synthetic tests run independently.

The repository contains software, tests, and research summaries. It does not contain the original count matrices, private machine configuration, or large local validation artifacts, and runtime analysis does not depend on the historical project directories.

## License and support

circaTE is distributed under the [MIT License](LICENSE.md), copyright (c) 2026 Phil Zhang. Third-party dependencies retain their own licenses.

Report reproducible software problems through [GitHub Issues](https://github.com/zhangfei1947/circaTE/issues). Include the method, input design, error message, and package versions, while excluding sensitive data. Maintainer: Phil Zhang, <philzhang@tamu.edu>.
