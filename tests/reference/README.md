# R oracle fixtures

`rOracle.json` is generated from the R packages gofCopula 0.4-3, copula
1.1-7, and VineCopula 2.6.1 (all available from CRAN; this repository does
not contain the R sources). The fixture freezes copula CDF/PDF values,
Rosenblatt transforms, raw test statistics, parameter estimates, and frozen
bootstrap chains. The MATLAB tests never require R at runtime.

To regenerate it, install the R packages and run the generator from the
repository root:

```bash
Rscript -e 'install.packages(c("gofCopula","jsonlite"), repos="https://cloud.r-project.org")'
Rscript tests/reference/generateROracle.R
```

To keep your default R library untouched, point `R_LIBS` at any scratch
directory (for example `R_LIBS="$(mktemp -d)"`) for both commands.

Regeneration must keep existing top-level JSON keys byte-identical; only
append new keys. The oracle tests' tolerances are calibrated to the package
versions above, and the generator records the versions it used in the
fixture's `metadata`.
