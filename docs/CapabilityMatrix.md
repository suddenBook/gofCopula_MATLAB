# Test and copula capability matrix

The authoritative machine-readable table is returned by `gofcopula.CopulaTestTable`. In the table below, `>=2` means any implemented dimension of at least two, `2` or `3` is the maximum supported dimension, and `-` means unsupported.

Family columns are abbreviated only in this document: Gaussian (`N`), Student t (`t`), Clayton (`C`), Gumbel (`G`), Frank (`F`), Joe (`J`), Ali-Mikhail-Haq (`AMH`), Galambos (`Gal`), Husler-Reiss (`HR`), Tawn (`Tawn`), t-EV (`TEV`), Farlie-Gumbel-Morgenstern (`FGM`), and Plackett (`Pl`).

| Tests | N | t | C | G | F | J | AMH | Gal | HR | Tawn | TEV | FGM | Pl |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `gofCvM`, `gofKS` | >=2 | >=2 | >=2 | >=2 | >=2 | >=2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 |
| `gofKendallCvM`, `gofKendallKS` | >=2 | >=2 | >=2 | >=2 | >=2 | >=2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 |
| `gofRosenblattSnB`, `SnC` | >=2 | >=2 | >=2 | >=2 | >=2 | >=2 | 2 | 2 | - | - | - | 2 | 2 |
| `gofRosenblattGamma`, `Chisq` | >=2 | >=2 | >=2 | >=2 | >=2 | >=2 | 2 | 2 | - | - | - | 2 | 2 |
| `gofKernel` | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 |
| `gofWhite` | 2 | 2 | 2 | 2 | 2 | 2 | - | - | - | - | - | - | - |
| `gofPIOSTn`, `gofPIOSRn` | 3 | 2 | 3 | 3 | 3 | 3 | 2 | 2 | - | - | - | 2 | 2 |
| `gofArchmSnB`, `SnC` | - | - | >=2 | >=2 | >=2 | >=2 | 2 | - | - | - | - | - | - |
| `gofArchmGamma`, `Chisq` | - | - | >=2 | >=2 | >=2 | >=2 | 2 | - | - | - | - | - | - |

The Power-Exponential copula (`powerexp`, PE) is an elliptical family that generalizes the Gaussian copula with a shape parameter beta (beta = 1 is Gaussian; stored in `DegreesOfFreedom`, read via `CopulaModel.Beta`). It supports every Gaussian-column test — `gofCvM`/`gofKS` and `gofRosenblatt*` at any dimension `>=2` (CvM/KS via a Monte-Carlo CDF; Rosenblatt via a numerical conditional CDF), `gofKendallCvM`/`gofKendallKS` and `gofPIOSTn`/`gofPIOSRn` at any `>=2`, and `gofKernel`/`gofWhite` in dimension two. Because PE is not closed under marginalization, beta is dimension-specific. See [PowerExponential.md](PowerExponential.md).

All listed families support bivariate empirical-CDF and custom-statistic bootstrapping. Rotations are supported only in dimension two and follow R's convention (90 degrees maps (u1,u2) to (1-u2,u1)). In `NumericMode="rCompatible"` the Student-t degrees of freedom are ceiled to an integer for `gofCvM`/`gofKS` and capped at 60 for `gofPIOSTn`, matching the R implementation's computational restrictions; `"corrected"` mode uses the fractional fitted df directly.

