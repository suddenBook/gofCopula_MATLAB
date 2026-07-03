# gofCopula for MATLAB

[![CI](https://github.com/suddenBook/gofCopula_MATLAB/actions/workflows/ci.yml/badge.svg)](https://github.com/suddenBook/gofCopula_MATLAB/actions/workflows/ci.yml)

Pure MATLAB source code for goodness-of-fit testing of parametric copula models. This project is a MATLAB migration of the R package `gofCopula` 0.4-3; it is not an installable MATLAB toolbox and requires no packaging or project file.

The default `NumericMode="corrected"` uses a calibrated parametric bootstrap in which every replicate follows the same margins pipeline as the observed data. `NumericMode="rCompatible"` reproduces the behavior of R `gofCopula` 0.4-3 for cross-language comparison. See [Numerics](docs/Numerics.md) for the exact differences.

## Requirements

- MATLAB R2025b
- Statistics and Machine Learning Toolbox
- Optimization Toolbox
- Parallel Computing Toolbox (optional, for parallel bootstrap execution)

R and WolframScript are only needed to regenerate development reference data; they are not runtime dependencies.

## Use from source

Clone or download the repository, then add `src` to the MATLAB path:

```matlab
repositoryRoot = "/path/to/gofCopula";
addpath(fullfile(repositoryRoot,"src"));
```

Run a test using one of the included datasets:

```matlab
sample = load(fullfile(repositoryRoot,"data","IndexReturns2D.mat"));
result = gofcopula.gofCvM("normal",sample.IndexReturns2D,M=99,Seed=42);
disp(result)
```

Bootstrap counts such as `M=99` are suitable for examples, not final inference. Use at least 1000 replicates and assess Monte Carlo uncertainty for reported results.

## Repository layout

| Path | Contents |
|---|---|
| `src/` | MATLAB source; `+gofcopula/` is the primary namespaced API |
| `data/` | Converted `.mat` datasets from the R package |
| `examples/` | Runnable MATLAB examples |
| `docs/` | Usage, migration, capability, numerical, and reference documentation |
| `tests/` | MATLAB tests and frozen cross-language reference values |
| `tools/` | Optional scripts for regenerating data and checking analytic derivations |
| `buildfile.m` | Static analysis, tests, and coverage automation |

The `+gofcopula` directory is MATLAB's package namespace. Public functions are called with the `gofcopula.` prefix, which prevents collisions with functions from other toolboxes.

## Main API

| MATLAB API | Purpose |
|---|---|
| `gofcopula.gof` | Run selected tests over one or more copula families |
| `gofcopula.gofco` | Run tests using a `gofcopula.CopulaModel` |
| `gofcopula.gofCvM`, `gofcopula.gofKS` | Empirical-copula tests |
| `gofcopula.gofKendallCvM`, `gofcopula.gofKendallKS` | Kendall-process tests |
| `gofcopula.gofRosenblatt*` | Rosenblatt-transform tests |
| `gofcopula.gofArchm*` | Archimedean-transform tests |
| `gofcopula.gofKernel`, `gofcopula.gofWhite` | Kernel and information-matrix tests |
| `gofcopula.gofPIOSTn`, `gofcopula.gofPIOSRn` | Pseudo in-and-out-of-sample tests |
| `gofcopula.CopulaTestTable` | Query supported family/test dimensions |

New to the package? Start with the **[User Guide](docs/UserGuide.md)** — a beginner-oriented manual covering usage, every option, and parameter tuning.

See also [Getting Started](docs/GettingStarted.m), the [R-to-MATLAB migration guide](docs/MigrationGuide.md), the [capability matrix](docs/CapabilityMatrix.md), the [Power-Exponential copula guide](docs/PowerExponential.md), and the [dataset guide](docs/Datasets.md).

## Development

From MATLAB, run:

```matlab
buildtool
```

This runs static analysis, the test suite, and the coverage check. Generated reports are written to `results/` and are ignored by Git. The same pipeline runs on every push and pull request through GitHub Actions (`.github/workflows/ci.yml`).

The frozen fixture in `tests/reference/rOracle.json` lets the MATLAB tests run without R. Regenerating it requires the R packages listed in `tests/reference/README.md`.

## Origin and citation

This code is a MATLAB migration of:

> Simon Trimborn, Ostap Okhrin, and Martin Waltz, `gofCopula` 0.4-3, *Goodness-of-Fit Tests for Copulae*, distributed through CRAN under GPL version 3 or later.

Original package: <https://cran.r-project.org/package=gofCopula>

The statistical implementation also follows the papers below:

1. Genest, C., Remillard, B., and Beaudoin, D. (2009), “Goodness-of-fit tests for copulas: A review and a power study.” <https://doi.org/10.1016/j.insmatheco.2007.10.005>
2. Zhang, S., Okhrin, O., Zhou, Q. M., and Song, P. X.-K. (2016), “Goodness-of-fit test for specification of semiparametric copula dependence models.” <https://doi.org/10.1016/j.jeconom.2016.02.017>

See [References](docs/References.md) for the complete bibliography. This repository does not redistribute the original R source tree or copyrighted paper text.

## License

Distributed under the GNU General Public License, version 3 or later. See [LICENSE](LICENSE). The original authors retain copyright in the R implementation; contributors retain copyright in their MATLAB modifications.
