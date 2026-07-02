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

All listed families support bivariate empirical-CDF and custom-statistic bootstrapping. Rotations are supported only in dimension two and follow R's convention (90 degrees maps (u1,u2) to (1-u2,u1)). In `NumericMode="rCompatible"` the Student-t degrees of freedom are ceiled to an integer for `gofCvM`/`gofKS` and capped at 60 for `gofPIOSTn`, matching the R implementation's computational restrictions; `"corrected"` mode uses the fractional fitted df directly.

