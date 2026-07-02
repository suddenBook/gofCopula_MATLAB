# R-to-MATLAB migration guide

The MATLAB API is namespaced to avoid collisions. Required data/model inputs are positional; optional inputs use MATLAB name-value syntax. Family and test values are normalized to strings.

## API mapping

| R gofCopula 0.4-3 | MATLAB |
|---|---|
| `gof(x, copula=..., tests=...)` | `gofcopula.gof(x,Copulas=...,Tests=...)` |
| `gofco(copulaobject, x, tests=...)` | `gofcopula.gofco(model,x,Tests=...)` |
| `gofCvM(copula, x, ...)` | `gofcopula.gofCvM(copula,x,...)` |
| `gofKS` | `gofcopula.gofKS` |
| `gofKendallCvM`, `gofKendallKS` | same names in `gofcopula` |
| `gofRosenblattSnB`, `SnC`, `Gamma`, `Chisq` | same names in `gofcopula` |
| `gofArchmSnB`, `SnC`, `Gamma`, `Chisq` | same names in `gofcopula` |
| `gofKernel`, `gofWhite` | same names in `gofcopula` |
| `gofPIOSTn`, `gofPIOSRn` | same names in `gofcopula` |
| `gofCustomTest` | `gofcopula.gofCustomTest` with a function handle |
| `gofCheckTime(copula, x, tests=...)` | `gofcopula.gofCheckTime(copula,x,Tests=...)` (two-point timing extrapolation; returns a `duration`) |
| `gofGetHybrid(result, p_values=, nsets=)` | `gofcopula.gofGetHybrid(result,PValues=,PValueNames=,NumberOfTests=)` |
| `gofOutputHybrid(result, tests=, nsets=)` | `gofcopula.gofOutputHybrid(result,Tests=<name substrings>,NumberOfTests=)` |
| `CopulaTestTable` | `gofcopula.CopulaTestTable` |
| `gofCopula4Test` | `gofcopula.gofCopula4Test` |
| `gofTest4Copula` | `gofcopula.gofTest4Copula` |
| `gofWhichCopula` (deprecated in R) | `gofcopula.gofCopula4Test` |
| `gofWhich` (deprecated in R) | `gofcopula.gofTest4Copula` |
| `print.gofCOP`, `plot.gofCOP` | `disp`/`plot` on `gofcopula.GofResult` |

## Option mapping

| R argument | MATLAB name-value |
|---|---|
| `param`, `param.est` | `Param`, `ParamEst` |
| `df`, `df.est` | `DF`, `DFEst` (also estimated for the t-EV copula) |
| `margins` | `Margins` (`"none"` replaces R's `margins=NULL`) |
| `flip` | `Flip` or `CopulaModel.Rotation` (recorded in `GofResult.Rotation`) |
| `M` (scalar or per-test vector in `gof`) | `M` (same) |
| `MJ` | `MJ` |
| `dispstr="ex"/"un"` | `Dispersion="exchangeable"/"unstructured"` |
| `m` | `BlockSize` |
| `delta.J` | `KernelScale` |
| `nodes.Integration` | `IntegrationNodes` |
| `lower`, `upper` | `Lower`, `Upper` (honored on every estimation path) |
| `seed.active` | `Seed` (scalar or M+1 vector) |
| `processes` | `Processes` |
| `priority` | `Priority` |
| `customTests` (function names) | `CustomTests` (cell of function handles) |
| `p_values`, `nsets` | `PValues`/`PValueNames`, `NumberOfTests` |

MATLAB adds `NumericMode`. The default is `"corrected"` (statistically calibrated bootstrap); `"rCompatible"` reproduces R 0.4-3 exactly, including its miscalibrated rank-margins bootstrap. See [Numerics](Numerics.md) for the complete semantics matrix — this is the single most important difference from both R and pre-0.5 versions of this implementation.

## Models and results

R copula objects are represented by `gofcopula.CopulaModel`. Its public properties record family, parameter vector, degrees of freedom, correlation structure, rotation, and whether parameters should be estimated.

Tests return `gofcopula.GofResult` objects. The `Tests` table has variables `Name`, `PValue`, and `Statistic`; fitted dependence, margin information, the applied rotation, and the numeric mode are stored in the other object properties. A multi-family `gof` call returns an object array, one element per fitted family. `gofGetHybrid`/`gofOutputHybrid` also return `GofResult` objects with metadata preserved, as their R counterparts return full gofCOP objects.

## Behavioral notes

- Rows are observations and columns are margins, as in R.
- Rank margins use `rank/(n+1)` with R's `ecdf` convention (ties take the upper rank).
- Rotations 90/180/270 degrees are bivariate only and follow R's `.rotateCopula` convention: 90 maps `(u1,u2)` to `(1-u2,u1)`, 270 to `(u2,1-u1)`.
- Replicate refitting is mode-dependent: `rCompatible` always refits (as R hardcodes); `corrected` refits only the quantities marked for estimation.
- A failed pseudo-likelihood fit falls back to exact inversion of Kendall's tau; a replicate whose estimation fails is redrawn deterministically. In `gof`, a failing test degrades to an NaN row (and poisons its hybrid combinations), as in R.
- Estimation guards mirror R: a negative Clayton estimate, a negative Frank estimate above dimension two, and estimates exactly at a family's domain boundary are errors.
- MATLAB random streams differ from R's RNG. A fixed MATLAB `Seed` is reproducible across repeated MATLAB calls and across serial/parallel modes (including statistic-internal sampling), but it does not reproduce an R seed's sample path, and sample paths may change between source releases.
- Custom statistics are function handles with contract `(u,copulaModel) -> finite real scalar`, where `u` contains the transformed pseudo-observations.
- Documented margin-parameterization deviations: R fits the `t` margin as `dt(x, df, ncp)` and `chisq` as `dchisq(x, df, ncp)`; MATLAB fits a location-scale t and a one-parameter chi-square. All other margin families match R's parameterization (Weibull reported as scale/shape by `fitdist`).
