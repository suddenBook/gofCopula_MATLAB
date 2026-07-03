# Serial dependence and autocorrelated data

The goodness-of-fit tests in this package calibrate their p-values with a
**parametric bootstrap that assumes independent observations**. When the rows of
the data are a time series — electrophysiology, finance, climate, sensor logs —
consecutive rows are autocorrelated, the effective sample size is far below the
row count, and the i.i.d. bootstrap under-estimates the variance of the test
statistic. The result is an **anti-conservative p-value: the test over-rejects**,
increasingly so as the autocorrelation grows.

`gofcopula.runTestSerial` is a first-class, family-agnostic remedy.

## 1. The problem

A copula goodness-of-fit statistic (Cramér–von Mises, Kolmogorov–Smirnov,
Rosenblatt, White, PIOS, …) is a functional of the empirical copula of the
pseudo-observations. Its sampling variance under the null is governed by the
*long-run* variance of the underlying rank process. For i.i.d. data that is the
ordinary variance; for serially dependent data it is inflated by roughly the
integrated autocorrelation time (IAT). The package's bootstrap reproduces the
i.i.d. variance, so on autocorrelated data the null distribution it builds is too
narrow and the observed statistic looks more extreme than it is.

Block-resampling the *fitted-copula draws* does not help — resampling an i.i.d.
series in blocks is still i.i.d. Phase-randomized surrogates do carry the
autocorrelation, but they are Gaussian by construction (a phase-randomized series
converges to a Gaussian process), so they only provide a *Gaussian*-copula null
and cannot calibrate a t, power-exponential, or Archimedean test.

## 2. The default method: decimate to near-independence

`runTestSerial` estimates how far apart the rows must be taken for the retained
rows to be approximately independent, keeps one row per that interval, and then
runs the ordinary i.i.d. bootstrap on the thinned sub-sample. Because it is pure
preprocessing in front of `gofcopula.runTest`, it works for **every** supported
copula family and any family added later.

The decorrelation length is the integrated autocorrelation time

```
tau = 1 + 2 * sum_{k>=1} rho(k)
```

computed by `gofcopula.internal.resampling.decorrelationLength`:

- **on the column ranks** (pseudo-observations), so it is invariant to the
  marginals, robust to heavy tails, and matched to the rank-based statistics;
- the autocovariance `rho(k)` is formed by FFT (Wiener–Khinchin) in base MATLAB —
  **no extra toolbox is required**;
- the sum uses Geyer's *initial positive sequence* (truncate at the first
  non-positive lag) to discard the noisy negative tail;
- a one-sided lag-1 significance gate (`rho(1) > 2/sqrt(n)`) keeps genuinely
  independent columns at `tau = 1` instead of over-thinning on autocorrelation
  noise;
- the per-column IATs are aggregated by their **maximum** (the most conservative
  choice), and the *same* interval is applied to every column so the
  cross-sectional dependence — the copula being tested — is preserved.

The interval is `round(tau)`, capped so the retained sub-sample keeps at least
`MinRetained` (default 50) rows on short series.

### Phase randomization: full power for the Gaussian test

`Method="phase"` keeps **all** rows. Instead of thinning, it calibrates the
statistic against **coherent phase-randomized surrogates** of the data. Each
surrogate keeps every column's power spectrum (its autocorrelation) and the
cross-spectrum between columns (the correlation matrix `R`, preserved exactly by
Parseval), but is given one shared random phase spectrum. Such a surrogate is a
draw from a stationary Gaussian process with the data's full second-order
structure — an exact Gaussian-copula null carrying the data's serial dependence —
so no rows are discarded and no power is lost.

Because phase randomization always yields a Gaussian process, this method is
valid **only for the `normal` copula**; it raises
`gofcopula:Serial:PhaseRequiresGaussian` for any other family, and requires
`Margins="ranks"` (the raw surrogates are re-ranked like the data). The
surrogates are built in base MATLAB (`fft`/`ifft`) and injected through
`runTest`'s `BootstrapSamples` hook.

## 3. Power trade-off

Decimation discards roughly a fraction `1 - 1/interval` of the rows. This trades
statistical **power for correct size**: the test becomes less able to detect a
true departure, but its false-positive rate returns to nominal. On long records
the cost is negligible (e.g. 800,000 rows with a decorrelation length of 50 still
leaves ~16,000 near-independent rows); on short records the `MinRetained` floor
keeps the sub-sample usable and, in the limit, collapses the interval to 1 (no
thinning). For the Gaussian copula, `Method="phase"` (above) avoids the trade-off
entirely by keeping all rows; a family-general full-power method (dependent
multiplier bootstrap) is heavier and statistic-specific — see §7.

## 4. Usage

```matlab
% X is an [n x d] matrix of time-series observations (rows = time).

% Family-agnostic (any copula), thinned to near-independence:
[result, serial] = gofcopula.runTestSerial("gofCvM", "clayton", X, M=999, Seed=1);

result.Tests.PValue     % serial-dependence-robust p-value
serial.thinInterval     % rows discarded between retained observations
serial.nThinned         % rows kept for the bootstrap (see the serial struct below)

% Gaussian test, full-power surrogate null (keeps all rows):
[result, serial] = gofcopula.runTestSerial("gofCvM", "normal", X, ...
    Method="phase", M=999, Seed=1);
```

`result` is an ordinary `gofcopula.GofResult`, identical in type to
`gofcopula.runTest`'s output, so it drops into the same reporting and plotting
code. Every `runTest` name-value option is accepted and forwarded unchanged
(`M`, `Seed`, `Processes`, `NumericMode`, `Dispersion`, `Margins`, …).

Additional options:

| Option | Default | Meaning |
|---|---|---|
| `ThinInterval` | `[]` (auto) | Fixed decimation interval (`"decimate"` only); overrides the automatic estimate. |
| `MaxLag` | `[]` (auto) | Maximum ACF lag used to estimate the interval (`"decimate"` only). |
| `Method` | `"decimate"` | `"decimate"` (default, all families), `"phase"` (full-power, `normal` only), or `"multiplier"` (reserved; raises `gofcopula:Serial:NotImplemented`). |

To inspect the decorrelation length without running a test:

```matlab
[interval, info] = gofcopula.internal.resampling.decorrelationLength(X);
```

## 5. Diagnostics (the `serial` struct)

| Field | Meaning |
|---|---|
| `method` | Correction method used (`"decimate"` or `"phase"`). |
| `nObserved` | Rows of `X` before thinning. |
| `nThinned` | Rows actually passed to the bootstrap (equals `nObserved` for `"phase"`). |
| `thinInterval` | Decimation interval applied. |
| `offset` | Starting row of the decimation (always 1). |
| `iatPerColumn` | Per-column integrated autocorrelation time (`NaN` if `ThinInterval` was supplied). |
| `maxLag` | ACF lag used to estimate the interval (`NaN` if supplied). |
| `keepIndices` | Retained row indices into `X`. |

## 6. Limitations

- Targets **positive, lag-1-dominant (AR-type) serial dependence**, the common
  case. Purely seasonal or higher-lag-only structure is not detected by the
  lag-1 gate.
- The interval is a single global value; strongly heterogeneous per-column
  autocorrelation is handled conservatively (the maximum), never by thinning
  columns differently (which would destroy the copula).
- On very short series the `MinRetained` floor forces `interval = 1`, i.e. no
  correction — there are simply not enough rows to both decorrelate and keep a
  usable sample.
- Decimation is a size-correction, not a power-preserving method (§3).

## 7. Future upgrades

`Method="multiplier"` is reserved for a **dependent (block) multiplier
bootstrap** of the empirical copula process (Bücher & Kojadinovic, 2016). It
would keep all rows for *any* family — the full-power, family-general
counterpart to `"phase"` — but it is a substantially larger, statistic-specific
addition (it does not reuse the parametric-bootstrap engine and applies mainly to
the empirical-copula CvM/KS statistics), so it currently raises
`gofcopula:Serial:NotImplemented`.

## 8. Worked example

`examples/serialDependenceExample.m` demonstrates the fix on synthetic data: on
an AR(1) Gaussian-copula null the i.i.d. bootstrap over-rejects (more so as the
autocorrelation grows) while `runTestSerial` holds ~5% size; `Method="phase"`
holds the same size while keeping every row; and on i.i.d. draws from several
families the correction is inert (interval ≈ 1, size unchanged).

## See also

`gofcopula.runTest`, `gofcopula.gofCvM`,
`gofcopula.internal.resampling.decorrelationLength`, `docs/UserGuide.md` §11.
