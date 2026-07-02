# Numerical conventions and precision

All observations, parameters, statistics, and probabilities are calculated as IEEE 754 binary64 (`double`). Inputs are rejected when NaN, infinity, invalid dimensions, or family parameter domains would make a statistic undefined.

Copula formulas use stable elementary functions such as `log1p` and `expm1` where cancellation is material, and the extreme-value Pickands functions with their first two derivatives are evaluated analytically (log-domain where powers such as `w^-theta` would overflow). Elliptical calculations use factorizations and linear solves rather than explicit matrix inverses. Correlation matrices must be symmetric positive definite. Probability inputs are kept away from exact zero and one when an inverse distribution or log-density requires an open unit interval.

## NumericMode

One switch selects between two complete, internally consistent regimes. `"corrected"` is the **default**; `"rCompatible"` reproduces gofCopula 0.4-3 (R) decision-for-decision for cross-language comparison.

| Behavior | `rCompatible` (R 0.4-3 parity) | `corrected` (default) |
|---|---|---|
| Bootstrap p-value | `count/M` (can be exactly 0) | `(count+1)/(M+1)` |
| Replicate margins | statistics computed on the raw copula draws, never re-ranked (as in R) | every replicate passes through the same margins pipeline as the observed data (pseudo-observations for `"ranks"`, inverse-then-refit for parametric margins) |
| Replicate refitting | parameters are always refitted, even with `ParamEst=false` (R hardcodes this) | refit only what the user asked to estimate |
| Student df in statistics | df ceiled to an integer for `gofCvM`/`gofKS` and capped at 60 for `gofPIOSTn` (R's computational restrictions) | fractional df used as fitted |
| White statistic | exact reproduction of `VineCopula::BiCopGofTest(method="white", B=0)` | generalized information-matrix sandwich with the nuisance-parameter influence adjustment |
| Gaussian PIOS parameterization | information matrices over all d(d-1)/2 pairwise correlations; Tn refits by Pearson correlation of the normal scores (as in R) | the fitted dispersion's own parameterization |

The selected mode is stored in every `gofcopula.GofResult`.

## Calibration (why the default changed)

The R package computes the observed statistic on pseudo-observations but the bootstrap replicate statistics on raw copula draws that are never re-ranked (`internal_bootstrap.R` in the original R package). Raw draws carry marginal sampling noise that ranked data do not, so replicate statistics are stochastically much larger than the observed one. Consequence, measured with this implementation: under a correctly specified null with `Margins="ranks"`, `rCompatible` p-values pile up near one (median > 0.9) and power against strongly misspecified alternatives collapses. The Genest & Remillard (2008) parametric bootstrap -- and `copula::gofPB` -- transforms every simulated sample exactly like the data. `corrected` mode implements that procedure; its null p-values are uniform and its power is restored (see `tests/calibrationTest.m`, which enforces both properties, and pins the rCompatible behavior so it cannot drift silently).

Use `rCompatible` only to reproduce results computed with the R package.

## White's test and VineCopula

For non-Student families, `rCompatible` reproduces VineCopula's published computation exactly, including its unconventional Hessian-like term `mean(-c'/c^2 + c''/c)` and uncentered variance. For the Student copula, VineCopula's source contains a per-observation influence correction that is multiplied by a gradient vector initialized to zero and never assigned -- dead code -- so the effective computation is the uncorrected covariance of the vech discrepancies, which is what this implementation computes (validated against VineCopula 2.6.1 output at n=8 and n=100). `corrected` mode provides the full influence adjustment for every family.

## Reproducibility

Set `Seed` to an integer for reproducible bootstrap sampling. Each replicate uses a deterministic derived stream (including statistic-internal sampling such as the kernel test's model draws), so serial and parallel execution produce identical results. A replicate whose estimation fails is redrawn deterministically from the next substream, as in R. Parallel Computing Toolbox is optional.

Sample paths for a given seed are stable within a source release but may change between versions (the 0.5.0 samplers replaced bisection inversion with exact closed-form and Marshall-Olkin methods). Reproducing an R seed's sample path is impossible by design; cross-language agreement is validated at the level of statistics and distributions instead.

Bootstrap p-values have Monte Carlo uncertainty even when deterministic. For an estimated p-value `p` based on `M` replicates, a useful standard-error estimate is `sqrt(p*(1-p)/M)`. Increase `M` until this uncertainty is adequate for the intended inference.

## Reference validation

The frozen fixtures in `tests/reference/rOracle.json` were generated from gofCopula 0.4-3, copula 1.1-7, and VineCopula 2.6.1. They cover every copula CDF/PDF, supported Rosenblatt transforms, parameter estimation (bivariate and trivariate), raw test statistics including the Kendall processes, White (two sample sizes), kernel integration, PIOS, rotation conventions, parametric margin fits, and two bootstrap chains with frozen samples (raw and rank-margins pipelines). Analytic derivations (extreme-value Pickands derivatives, generator-derivative representations, closed-form conditional inverses) are verified by WolframScript scripts committed under `tools/wolfram/` with their outputs, and cross-checked against finite differences in the test suite. Neither R nor Wolfram Engine is required at runtime.
