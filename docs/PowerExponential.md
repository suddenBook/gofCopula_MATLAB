# Power-Exponential copula (`powerexp`)

The power-exponential (PE) copula is the elliptical copula obtained from the
multivariate power-exponential distribution of Gómez, Gómez-Villegas & Marín
(1998). It generalizes the Gaussian copula with a single **shape parameter β**:
β = 1 is exactly the Gaussian copula, β < 1 gives heavier tails, and β > 1 lighter
(more uniform-like) tails. Like the Gaussian copula it is radially symmetric and
tail-independent, but β lets it depart from Gaussian dependence while keeping the
same correlation matrix.

## Parameters

| Symbol | Meaning | Stored in |
|---|---|---|
| `R` (θ) | correlation matrix (unit diagonal) | `CopulaModel.Theta` |
| β | shape (β = 1 ⇒ Gaussian) | `CopulaModel.DegreesOfFreedom`, read via `CopulaModel.Beta` |

β reuses the model's existing second-scalar slot (`DegreesOfFreedom`), so it flows
through the same estimation, bootstrap and public-API plumbing as the Student-t
degrees of freedom. The read-only `Beta` property is a clearer alias.

```matlab
sample = load("data/IndexReturns2D.mat");
% PIOS test (the recommended density-based test for PE):
r = gofcopula.gofPIOSTn("powerexp", sample.IndexReturns2D, M=1000, Seed=42);
% direct primitives:
U = gofcopula.copulaRandom("powerexp", 500, 0.5, Dimension=2, DF=0.7);  % DF carries beta
c = gofcopula.copulaPDF("powerexp", U, 0.5, DF=0.7, Log=true);
```

## Not marginally closed — the key caveat

A univariate margin of a d-dimensional PE distribution is **not** itself a
univariate PE (unlike the Gaussian and Student-t, whose margins stay in-family).
Consequently:

- The copula density cannot reuse the Student-t template (`tinv`); it needs the
  **true** elliptical margin, built numerically from the joint generator.
- The **d-dimensional** PE copula is not built from lower-dimensional PE copulas:
  its bivariate margins are elliptical but not bivariate PE copulas. **β is
  therefore dimension-specific and should not be compared across dimensions.**
- A goodness-of-fit test at a fixed dimension d is still perfectly well posed
  (fit `PE_d`, bootstrap from `PE_d`, refit `PE_d`); only cross-dimensional
  reuse is disallowed.

## How it is computed

The marginal generator, quantile `Q1` and CDF `F1` are built in
`+internal/+elliptical/peMarginals.m` by evaluating the generator
`exp(−s^β/2)` analytically inside a fine radial quadrature (the Cambanis–Huang–Simons
marginalization), then shape-preserving interpolation. At β = 1 this reproduces the
standard-normal margin to machine precision.

| Primitive | File | Notes |
|---|---|---|
| marginal transform (g₁, Q1, F1) | `peMarginals.m` | analytic generator; machine-precision at β=1 |
| copula log-density | `peCopulaLogPDF.m` | `c(u)=\|R\|^{-1/2} g_d(z'R⁻¹z)/∏ g₁(z_j²)`, `z=Q1(u)` |
| sampler | `peCopulaRandom.m` | stochastic representation, radius `T^{1/2β}`, `T~Gamma(d/2β,2)` |
| CDF | `peCopulaCDF.m` | Monte Carlo (no closed form); deterministic seed |
| conditional CDF (Rosenblatt) | `peConditionalCDF.m` | conditional generator `exp(−½(t+q)^β)` marginalized to 1-D |

**Estimation** (`estimateModel.m`): R by inversion of Kendall's τ
(`R̂ = sin(πτ̂/2)`, shape-independent; projected to the nearest correlation matrix
in d > 2 if needed), then β by 1-D maximum pseudo-likelihood holding R fixed.

## Supported tests

| Test | d | Primitive used | Notes |
|---|---|---|---|
| `gofPIOSTn`, `gofPIOSRn` | ≥ 2 | density | **recommended**; PIOS-Rn info matrix is β-aware |
| `gofKendallCvM`, `gofKendallKS` | ≥ 2 | sampler | needs no CDF or density |
| `gofRosenblattSnB/SnC/Gamma/Chisq` | ≥ 2 | conditional CDF | d > 2 steps use a (q,r) tabulation |
| `gofCvM`, `gofKS` | ≥ 2 | Monte-Carlo CDF | the review's `S_n`; slower |
| `gofKernel` | 2 | sampler | bivariate only |
| `gofWhite` | 2 | density | bivariate; info matrix β-aware |

## Accuracy and resolution

Validated to machine precision at β = 1 (log g₁ vs N(0,1): 5e-15; F1/Q1 vs
`normcdf`/`norminv`: ~1e-7/1e-6; copula density vs `copulapdf('Gaussian')`: ~7e-6)
and to ~1e-9 against an independent `integral()` reference for β ≠ 1. Resolution
is tunable via `peMarginals` name-value options (`GridPoints`, `QuadraturePoints`,
`MaxRadius`) and `peCopulaCDF`'s `MonteCarloSize`. The marginal grid is
sinh-spaced (fine near the origin, exponentially coarse in the tail), so the very
heavy tails of small β are covered at a fixed point budget: β is estimated on
`[0.2, 5]` with roughly constant per-fit cost across that range.

## References

- Gómez, Gómez-Villegas & Marín (1998), *A multivariate generalization of the
  power exponential family of distributions*, Comm. Statist. Theory Methods 27(3).
- Genest, Rémillard & Beaudoin (2009); Zhang, Okhrin, Zhou & Song (2016) — the
  goodness-of-fit methodology reused here.
- Derumigny & Fermanian (2022), `ElliptCopulas` — the generator toolkit ported in
  `+internal/+elliptical/` and used for cross-validation.
