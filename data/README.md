# Bundled datasets

These MAT-files are exact binary64 conversions of the `.RData` files in gofCopula 0.4-3. R emitted each numeric array with `writeBin(..., size=8, endian="little")`; MATLAB read those bytes as `double`. No decimal intermediate or numerical transformation was used, so the stored values are bit-identical to R's. Reference values are verified by `tests/dataTest.m`.

| File | Main variable | Contents |
|---|---|---|
| `Banks.mat` | `Banks` | Annual Citigroup/Bank of America adjusted log returns, 2004--2012 |
| `CryptoCurrencies.mat` | `CryptoCurrencies` | Annual Bitcoin/Litecoin adjusted log returns, 2015--2018 |
| `IndexReturns2D.mat` | `IndexReturns2D` | 100-by-2 DAX/SMI log returns for 1998 |
| `IndexReturns3D.mat` | `IndexReturns3D` | 200-by-3 DAX/SMI/CAC log returns for 1998 |

`Banks` and `CryptoCurrencies` are scalar structures with fields `Years`, `VariableNames`, and `Data`. `Data{k}` corresponds to `Years(k)`. The index files contain a numeric matrix and a `VariableNames` string vector. Every file also contains a `Provenance` structure recording its source and conversion.

The source attributions are Yahoo Finance (banks), CoinMetrics (cryptocurrencies), and Erste Bank AG, Vienna (index returns).

To reproduce the conversion from the repository root, export the arrays
from R into a scratch folder of your choice, then rebuild the MAT-files:

```text
Rscript tools/exportRData.R <scratch-folder>
```

```matlab
addpath("tools")
importRData("<scratch-folder>","data")
```

The export step reads the datasets from an installed copy of R `gofCopula` 0.4-3; the original R package source is not stored in this repository.
