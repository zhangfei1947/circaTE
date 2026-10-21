# Validation of circaTE 0.1.0

Validated locally with R 4.4.2, edgeR 4.4.2, limma 3.62.2, DESeq2 1.46.0 and statmod 1.5.0.

- 189 local test assertions passed with no failures, warnings or skips in the development workspace.
- `R CMD check --no-manual --no-build-vignettes` passed with no errors, warnings or notes. The PDF manual was not built.
- Installed-package edgeR and weighted AR analyses across five real groups (10 fits) exactly reproduced the frozen reference p-values and amplitudes.
- Installed CLI success and invalid-option paths passed.
- Renaming the package to circaTE was followed by a fresh build, installation, unit tests, package check and CLI run.

The research data and historical model directories are deliberately not distributed. Two checkout-specific regression tests skip when those historical files are absent. Synthetic input, fitting, calibration and design checks are self-contained. Tests reproduce software behavior; they do not establish universal error control in new biological data.

See [modeling history](inst/doc/MODELING_HISTORY.md) and [scientific scope](inst/doc/VALIDATION_SCOPE.md) for empirical results, uncertainty and limits. Historical source references in those documents identify the development evidence and are not runtime dependencies.
