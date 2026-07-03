# gofCopula for MATLAB — User Guide

A beginner-to-intermediate manual for goodness-of-fit (GoF) testing of copulas.
It assumes no prior copula experience. If you just want to copy-paste something
that works, read [Quick start](#2-quick-start) and [Choosing a test](#7-choosing-a-test-and-family).

**Contents**

1. [What this package does](#1-what-this-package-does)
2. [Quick start](#2-quick-start)
3. [Core ideas in 5 minutes](#3-core-ideas-in-5-minutes)
4. [Three ways to run a test](#4-three-ways-to-run-a-test)
5. [Reading the result](#5-reading-the-result)
6. [Supported families and tests](#6-supported-families-and-tests)
7. [Choosing a test and family](#7-choosing-a-test-and-family)
8. [Full parameter reference](#8-full-parameter-reference)
9. [Parameter tuning guide](#9-parameter-tuning-guide)
10. [The Power-Exponential copula](#10-the-power-exponential-copula)
11. [Time series / autocorrelated data](#11-time-series--autocorrelated-data)
12. [Margins](#12-margins)
13. [Troubleshooting](#13-troubleshooting)
14. [Recipes](#14-recipes)
15. [Low-level building blocks](#15-low-level-building-blocks)

---

## 1. What this package does

A **copula** describes the *dependence* between variables, separately from each
variable's own distribution. This package answers one question:

> *Is the dependence in my data consistent with a particular copula family
> (Gaussian, t, Clayton, …)?*

It returns a **p-value** for that hypothesis, computed by a parametric bootstrap.
Small p-value ⇒ the family is a poor fit. Typical uses: choosing a dependence
model for financial returns, risk aggregation, or validating an assumed copula.

**Requirements:** MATLAB R2025b, Statistics and Machine Learning Toolbox,
Optimization Toolbox. Parallel Computing Toolbox is optional (faster bootstrap).

**Setup** — add the `src` folder to your path (there is nothing to install):

```matlab
repositoryRoot = "/path/to/gofCopula_MATLAB";
addpath(fullfile(repositoryRoot,"src"));
```

Every public function lives in the `gofcopula.` namespace, e.g.
`gofcopula.gofCvM(...)`.

---

## 2. Quick start

```matlab
% Load a bundled dataset: 100 days of DAX/SMI log returns (rows = days).
s = load(fullfile(repositoryRoot,"data","IndexReturns2D.mat"));
x = s.IndexReturns2D;                       % 100-by-2 numeric matrix

% Test whether the dependence is Gaussian. M = number of bootstrap replicates.
result = gofcopula.gofCvM("normal", x, M=1000, Seed=42);
disp(result)

% The p-value:
p = result.Tests.PValue
```

`x` must be a numeric matrix: **rows are observations, columns are variables**
(≥ 2 columns). `Seed` makes the run reproducible. If `p ≥ 0.05` you have *no
evidence against* the Gaussian copula; if `p < 0.05` you reject it.

> Start with `M=1000`. Smaller values (e.g. `M=19`) are only for quick trials —
> their p-values are too coarse for real conclusions (see [§9](#9-parameter-tuning-guide)).

---

## 3. Core ideas in 5 minutes

- **Pseudo-observations.** The copula ignores each variable's marginal scale, so
  the data is first rank-transformed column by column to `rank/(n+1) ∈ (0,1)`.
  This is done automatically (`Margins="ranks"`).
- **Fit.** The chosen family's parameters are estimated from the data (the
  correlation matrix, plus a second parameter such as the t degrees of freedom
  or the Power-Exponential shape β).
- **Statistic.** A number measuring the gap between the *empirical* copula of the
  data and the *fitted* family copula. Different tests use different gaps.
- **Parametric bootstrap.** Because the parameters were fitted to the data, the
  observed statistic is systematically small. So the reference distribution is
  built by simulation: draw `M` fresh samples from the fitted copula, **re-fit**
  each one, and recompute the statistic. The p-value is the fraction of
  simulated statistics at least as extreme as the observed one.
- **Decision.** Reject the family when `p < α` (commonly `α = 0.05`). The test
  can only *falsify*: `p ≥ α` means "not contradicted", never "proven correct".

---

## 4. Three ways to run a test

**(a) One test, one family** — the simplest entry point:

```matlab
result = gofcopula.gofCvM("t", x, M=1000, Seed=42);      % Cramer–von Mises, Student-t
```

Every test has a wrapper: `gofCvM`, `gofKS`, `gofKendallCvM`, `gofKendallKS`,
`gofRosenblattSnB/SnC/Gamma/Chisq`, `gofArchmSnB/SnC/Gamma/Chisq`, `gofKernel`,
`gofWhite`, `gofPIOSTn`, `gofPIOSRn`. All accept the same options ([§8](#8-full-parameter-reference)).

**(b) Many families and/or many tests at once** — `gofcopula.gof`:

```matlab
result = gofcopula.gof(x, ...
    Copulas=["normal","t","clayton","gumbel"], ...
    Tests=["gofCvM","gofKendallCvM"], ...
    M=1000, Seed=42);
disp(result)                    % one GofResult per family, incl. combined rows
```

Omit `Copulas` to test every family that supports your dimension; omit `Tests`
to run every compatible test. Unsupported combinations are skipped with a
warning; a failing test degrades to a `NaN` row instead of aborting the batch.
`M` may be a vector (one bootstrap count per test).

**(c) From an explicit model** — when you want to fix or bound parameters,
choose a dispersion structure, or apply a rotation, build a `CopulaModel`:

```matlab
model  = gofcopula.CopulaModel("t", Theta=0.5, DegreesOfFreedom=6, ...
             EstimateTheta=true, EstimateDegreesOfFreedom=true, ...
             Dispersion="unstructured");
result = gofcopula.gofco(model, x, Tests=["gofPIOSTn","gofKendallCvM"], M=1000);
```

---

## 5. Reading the result

Every call returns a `gofcopula.GofResult`. The useful fields:

| Field | Meaning |
|---|---|
| `.Tests` | a table with columns `Name`, `PValue`, `Statistic` (one row per test) |
| `.Theta` | the fitted copula parameter(s) — a correlation or correlation vector |
| `.DegreesOfFreedom` | the fitted second parameter (t df, or PE shape β) |
| `.Copula`, `.Method`, `.NumericMode` | bookkeeping |

```matlab
result.Tests                     % table of p-values and statistics
result.Tests.PValue(1)           % the first test's p-value
result.Theta                     % fitted correlation
result.DegreesOfFreedom          % fitted df / beta
plot(result)                     % scatter of p-values vs the 5% line
```

**Interpreting `PValue`:** reject the family if it is below your level `α`.
A `NaN` p-value means the test could not be computed (see the warning text).

---

## 6. Supported families and tests

Fourteen families: `normal` (Gaussian), `t`, `clayton`, `gumbel`, `frank`,
`joe`, `amh`, `galambos`, `huslerreiss`, `tawn`, `tev`, `fgm`, `plackett`, and
`powerexp` (Power-Exponential — see [§10](#10-the-power-exponential-copula)).
(`"gaussian"` is accepted as an alias for `"normal"`.)

Which test works for which family and dimension depends on all three. **Query it
before a big run** instead of guessing:

```matlab
gofcopula.CopulaTestTable()                 % full family x test table
gofcopula.gofTest4Copula("powerexp", 3)     % tests available for PE in 3-D
gofcopula.gofCopula4Test("gofWhite")        % families supported by the White test
```

See [CapabilityMatrix.md](CapabilityMatrix.md) for the printed table. Rough summary:

- **Elliptical** (`normal`, `t`, `powerexp`): all tests, most in any dimension ≥ 2.
- **Archimedean** (`clayton`, `gumbel`, `frank`, `joe`, `amh`): empirical, Kendall,
  Rosenblatt, PIOS, plus the Archimedean-transform tests.
- **Extreme-value / one-parameter** (`galambos`, `huslerreiss`, `tawn`, `tev`,
  `fgm`, `plackett`): mostly bivariate only.

---

## 7. Choosing a test and family

**If you don't know which test to use**, these are good defaults:

| Situation | Recommended test | Why |
|---|---|---|
| General 2-D check | `gofCvM` | classic, compares the whole copula |
| Dimension `d > 2` | `gofPIOSTn` or `gofKendallCvM` | keep power as `d` grows |
| You want a fast, CDF-free test | `gofKendallCvM` | uses only a sampler |
| Density-based / higher power | `gofPIOSTn` | leave-block-out likelihood |

**What each test family does:**

- **`gofCvM`, `gofKS`** — distance between the empirical copula and the fitted
  copula's CDF (Cramér–von Mises = integrated squared, Kolmogorov–Smirnov = max).
- **`gofKendallCvM`, `gofKendallKS`** — same idea on the *Kendall function*
  `K(t)=P(C(U) ≤ t)`; needs only a sampler, so it is robust and dimension-friendly.
- **`gofRosenblatt*`** — apply the Rosenblatt transform (which turns the copula
  into independent uniforms under the null), then test that the output is
  uniform/independent (`SnB`, `SnC`) or map it to a Gamma/Chi-square score.
- **`gofArchm*`** — the Archimedean analogue of the Rosenblatt tests (Archimedean
  families only).
- **`gofKernel`** — compares kernel density estimates of the data and a simulated
  model sample (bivariate only).
- **`gofWhite`** — White's information-matrix equality test (bivariate).
- **`gofPIOSTn`, `gofPIOSRn`** — "pseudo in-and-out-of-sample": `Tn` uses
  leave-block-out refitting, `Rn` an information-matrix ratio. Density-based and
  strong in higher dimensions.

**Comparing families:** run `gofcopula.gof(x, Tests="gofCvM", M=1000)` across all
families and pick the one that is *not* rejected (largest p-value). See
[Recipe 14.2](#14-recipes).

---

## 8. Full parameter reference

Every `gof*` wrapper, `gof`, `gofco`, and `runTest` share these name–value
options (defaults in parentheses).

| Option | Default | Meaning |
|---|---|---|
| `Param` | `0.5` | starting value of the copula parameter θ (a correlation for elliptical; the family parameter otherwise). Also the *fixed* value when `ParamEst=false`. |
| `ParamEst` | `true` | estimate θ from the data. Set `false` to test a fully specified copula. |
| `DF` | `4` | the family's **second** parameter: Student-t degrees of freedom, or Power-Exponential shape **β**. Ignored by families without one. |
| `DFEst` | `true` | estimate `DF` (t df or PE β). Set `false` to fix it at `DF`. |
| `Margins` | `"ranks"` | how columns become uniforms. `"ranks"` (recommended), `"none"` (data already in `[0,1]`), or a parametric fit (`"norm"`, `"t"`, `"gamma"`, …). See [§12](#12-margins). |
| `Flip` | `0` | copula rotation in degrees: `0/90/180/270` (bivariate only). |
| `M` | `1000` | number of bootstrap replicates (⇒ p-value resolution). In `gof`, may be a vector (one per test). |
| `MJ` | `100` | Monte-Carlo sample size used inside the kernel test. |
| `Dispersion` | `"exchangeable"` | elliptical correlation structure: `"exchangeable"` (one shared ρ) or `"unstructured"` (a full matrix). Matters only for `normal`/`t`/`powerexp` with `d > 2`. |
| `BlockSize` | `1` | leave-block-out block length for `gofPIOSTn` (must divide `n`). |
| `KernelScale` | `0.5` | bandwidth of the `gofKernel` density estimate. |
| `IntegrationNodes` | `12` | Gauss–Legendre nodes for the kernel-test integral. |
| `Lower`, `Upper` | `[]` | bounds on θ during estimation. |
| `Seed` | `[]` | RNG seed. **Set it for reproducible p-values.** |
| `Processes` | `1` | parallel workers for the bootstrap (needs Parallel Computing Toolbox). |
| `NumericMode` | `"corrected"` | `"corrected"` = calibrated p-values (recommended). `"rCompatible"` = reproduce R gofCopula 0.4-3. |

`gof` adds `Copulas`, `Tests`, `Priority` (`"copula"`/`"tests"`), and
`CustomTests`. `runTest` adds `CustomTest`, `ModelSample`, and `BootstrapSamples`
(inject your own `n×d×M` replicate array — the hook for a block bootstrap, [§11](#11-time-series--autocorrelated-data)).

---

## 9. Parameter tuning guide

### 9.1 `M` — the single most important knob

The p-value is `(#{replicate ≥ observed} + 1) / (M + 1)`. Its Monte-Carlo
standard error is about `sqrt(p(1−p)/M)`:

| `M` | approx. SE of a p-value near 0.5 | use for |
|---|---|---|
| 19–99 | 0.05–0.11 | quick trials only |
| 1000 | ~0.016 | **default for reporting** |
| 10000 | ~0.005 | borderline cases near α |

Rule of thumb: pick `M` so the SE is small relative to how close `p` is to `α`.
If `p ≈ 0.06` and you must decide at `α = 0.05`, use `M = 10000`. Runtime scales
linearly with `M`; use `Processes` to parallelize.

### 9.2 `Processes` — parallel speed-up

Set `Processes` to your number of physical cores (needs Parallel Computing
Toolbox). Results are **identical** to serial for a given `Seed` (each replicate
uses its own seeded substream), so parallelism is free accuracy-wise:

```matlab
gofcopula.gofPIOSTn("normal", x, M=10000, Seed=42, Processes=8);
```

### 9.3 Which test — power vs cost

- Cheapest: `gofKendallCvM`, `gofKS`, `gofCvM` (sampler / CDF only).
- Most expensive: `gofPIOSTn` (refits the model for every leave-out block *and*
  every replicate) and, for `powerexp`, the Monte-Carlo-CDF tests. They are also
  among the most powerful. Budget accordingly, or use `Processes`.

### 9.4 `Dispersion` (elliptical, `d > 2`)

`"exchangeable"` estimates a single correlation shared by all pairs — fast and
stable, good when pairwise dependence is similar. `"unstructured"` estimates the
full correlation matrix — more flexible, needs more data. For `d = 2` they are
equivalent.

### 9.5 `BlockSize` (`gofPIOSTn`)

Controls the leave-block-out jackknife. `BlockSize` **must divide `n`**. Larger
blocks are faster (fewer refits) but coarser. `1` (leave-one-out) is the most
thorough; a divisor giving ~10–20 blocks is a reasonable speed/quality trade.

### 9.6 `Lower` / `Upper`

Constrain estimation when a family's parameter is near a boundary or a fit
wanders (e.g. keep a correlation in `[-0.9, 0.9]`). Also the fix for a
"boundary estimate" error ([§13](#13-troubleshooting)).

### 9.7 Kernel-test knobs

`KernelScale` (bandwidth) and `IntegrationNodes` affect only `gofKernel`.
Increase `IntegrationNodes` (e.g. to 20–32) for a smoother integral; adjust
`KernelScale` if the density comparison looks over/under-smoothed.

### 9.8 `NumericMode`

Keep `"corrected"` unless you are specifically cross-checking against R
gofCopula 0.4-3, in which case use `"rCompatible"` (see [Numerics.md](Numerics.md)).

---

## 10. The Power-Exponential copula

`powerexp` is an **elliptical** family that generalizes the Gaussian copula with
one extra shape parameter **β** (stored in the `DF` slot):

- **β = 1** is exactly the Gaussian copula;
- **β < 1** gives heavier joint tails, **β > 1** lighter (more uniform-like).

```matlab
% Estimate R and beta, then test — beta is estimated automatically (DFEst=true):
result = gofcopula.gofPIOSTn("powerexp", x, M=1000, Seed=42);
Rhat    = result.Theta                 % fitted correlation
betaHat = result.DegreesOfFreedom      % fitted shape beta

% Test a copula with a FIXED beta (e.g. a heavy-tailed beta = 0.6):
gofcopula.gofCvM("powerexp", x, DF=0.6, DFEst=false, M=1000, Seed=42);
```

Estimation follows a two-step recipe: the correlation `R` by inverting Kendall's
τ (`R̂ = sin(π τ̂ / 2)`, projected to the nearest correlation matrix if needed),
then β by 1-D maximum pseudo-likelihood.

**Supported tests:** all elliptical tests — `gofPIOSTn/Rn` (recommended),
`gofKendallCvM/KS`, `gofCvM/KS` (via a Monte-Carlo CDF), `gofRosenblatt*`,
`gofKernel`/`gofWhite` (bivariate). See [PowerExponential.md](PowerExponential.md).

### 10.1 Tuning the Power-Exponential numerics

The PE margins are built numerically (a PE margin is not itself PE). The defaults
are calibrated to be far more accurate than the bootstrap noise, so **you rarely
need to touch them**. β is estimated on `[0.2, 5]` at roughly constant cost
thanks to a sinh-spaced grid. Advanced knobs (used by calling the internal
`peMarginals`/`peCopulaCDF` directly):

| Knob | Default | When to raise it |
|---|---|---|
| `GridPoints` (marginal grid) | 8000 | β below ~0.2, or you want tighter tails |
| `QuadraturePoints` (radial) | 8000 | same |
| `MaxRadius` | auto from `n`, β | force a wider tail than the data quantiles need |
| `MonteCarloSize` (CDF) | 5×10⁵ | tighter `gofCvM`/`gofKS` CDF (error ~1/√size) |

All are linear in cost. Accuracy is validated to ~1e-7 (mass to 5–7 digits) and
to machine precision at β = 1.

---

## 11. Time series / autocorrelated data

**Important assumption:** the parametric bootstrap treats the rows of `x` as
**independent**. For serially dependent data (e.g. raw high-frequency returns,
EEG), the p-values can be miscalibrated. Remedies:

1. **Use `gofcopula.runTestSerial`** (built-in, recommended): it decimates the
   data to near-independence and then runs the ordinary bootstrap on the
   sub-sample — family-agnostic, one call. For the `normal` copula,
   `Method="phase"` instead calibrates against coherent phase-randomized
   surrogates and keeps every row (no power loss). See `docs/SerialDependence.md`.
2. **Pre-whiten the margins** (common for financial data): fit an ARMA–GARCH
   model per column and feed the standardized residuals to the test. The bundled
   `Banks`/`CryptoCurrencies` datasets are already volatility-adjusted this way.
3. **Block bootstrap by hand:** generate your own block-resampled replicate array
   of size `n×d×M` and pass it via `runTest(..., BootstrapSamples=yourArray)`.
   `runTestSerial` automates the correct version of this for you.

Quick self-check with the Econometrics Toolbox (per column): `lbqtest` for serial
correlation and `archtest` for volatility clustering. If both are insignificant,
the plain bootstrap is fine.

---

## 12. Margins

`Margins` controls how each column is mapped to `(0,1)` before the copula is
fitted:

- **`"ranks"`** (default) — nonparametric pseudo-observations `rank/(n+1)`. Use
  this unless you have a reason not to.
- **`"none"`** — the data is *already* uniform in `[0,1]` (e.g. probabilities);
  it is used as-is.
- **Parametric** — fit a named marginal and apply its CDF: `"norm"`, `"t"`,
  `"gamma"`, `"lnorm"`, `"weibull"`, `"exp"`, `"beta"`, `"cauchy"`, `"chisq"`,
  `"f"`. May be a single name (all columns) or a string array (one per column):

```matlab
gofcopula.gofCvM("normal", x, Margins=["ranks","norm"]); % col 1 nonparametric, col 2 Normal
```

---

## 13. Troubleshooting

| Message / symptom | Cause and fix |
|---|---|
| `... is not supported for <family> in dimension d` | that test/family/dimension combination is unavailable. Query `gofcopula.gofTest4Copula(family,d)` and pick a supported test. |
| `BlockSize must divide the number of observations` | choose a `BlockSize` that divides `n` (e.g. `n=100` ⇒ 1, 2, 4, 5, 10, 20, 25, 50). |
| `The estimated <family> parameter lies at its domain boundary` | the fit hit a boundary (unstable). Supply `Lower`/`Upper`, try another family, or check the data. |
| `Data with Margins='none' must lie in [0,1]` | you passed raw data with `Margins="none"`. Use `Margins="ranks"` instead. |
| `At least two data columns are required` | `x` must have ≥ 2 columns (rows = observations). |
| p-value is `NaN` in a `gof` batch | that single test failed (see the warning); the batch continued. |
| Test is very slow | reduce `M`, set `Processes`, prefer a cheaper test ([§9.3](#9-parameter-tuning-guide)), or increase `BlockSize` for PIOS. |

---

## 14. Recipes

### 14.1 A reportable single test

```matlab
s = load(fullfile(repositoryRoot,"data","IndexReturns2D.mat"));
r = gofcopula.gofCvM("t", s.IndexReturns2D, M=10000, Seed=1, Processes=8);
fprintf("t-copula: p = %.3f, fitted df = %.1f\n", r.Tests.PValue, r.DegreesOfFreedom);
```

### 14.2 Compare several families and pick the best

```matlab
x = s.IndexReturns2D;
r = gofcopula.gof(x, Copulas=["normal","t","clayton","gumbel","frank","powerexp"], ...
        Tests="gofCvM", M=1000, Seed=7);
for k = 1:numel(r)
    fprintf("%-10s p = %.3f\n", r(k).Copula, r(k).Tests.PValue(1));
end
% The family with the largest p-value is the least contradicted by the data.
```

### 14.3 Test Power-Exponential vs Gaussian

```matlab
pG  = gofcopula.gofPIOSTn("normal",   x, M=1000, Seed=3).Tests.PValue(1);
rPE = gofcopula.gofPIOSTn("powerexp", x, M=1000, Seed=3);
fprintf("Gaussian p = %.3f | PE p = %.3f, beta = %.3f\n", ...
    pG, rPE.Tests.PValue(1), rPE.DegreesOfFreedom);
% beta far from 1 with Gaussian rejected but PE not rejected = a non-Gaussian,
% still-elliptical dependence.
```

### 14.4 Higher dimension

```matlab
x3 = load(fullfile(repositoryRoot,"data","IndexReturns3D.mat")).IndexReturns3D;   % 200-by-3
gofcopula.gofKendallCvM("normal", x3, Dispersion="unstructured", M=1000, Seed=9);
```

---

## 15. Low-level building blocks

For simulation or custom pipelines you can call the primitives directly:

```matlab
U = gofcopula.copulaRandom("t", 500, 0.5, Dimension=2, DF=5);     % sample a copula
c = gofcopula.copulaPDF("t", U, 0.5, DF=5, Log=true);            % (log) density
C = gofcopula.copulaCDF("normal", U, 0.5);                        % CDF
Z = gofcopula.rosenblatt("normal", U, 0.5);                       % Rosenblatt transform
m = gofcopula.CopulaModel("clayton", Theta=2);                    % a reusable model
```

These use the same `Dispersion`, `DF`, and `Rotation` options as the tests.

---

**See also:** [GettingStarted.m](GettingStarted.m) · [CapabilityMatrix.md](CapabilityMatrix.md) ·
[PowerExponential.md](PowerExponential.md) · [Datasets.md](Datasets.md) ·
[MigrationGuide.md](MigrationGuide.md) (R option names) · [Numerics.md](Numerics.md).
