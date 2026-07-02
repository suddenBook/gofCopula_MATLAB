# Dataset reference

The four datasets from gofCopula 0.4-3 are stored under the repository's `data/` directory. Locate it without depending on the current folder:

```matlab
sourceRoot = fileparts(fileparts(which("gofcopula.CopulaModel")));
repositoryRoot = fileparts(sourceRoot);
s = load(fullfile(repositoryRoot,"data","IndexReturns2D.mat"));
```

## Banks

Volatility-adjusted Citigroup (`C`) and Bank of America (`BoA`) log returns, split into years 2004--2012. The source R object is a nine-element named list. The MATLAB scalar structure has fields:

- `Years`: 9-by-1 numeric year vector;
- `VariableNames`: `["C" "BoA"]`;
- `Data`: 9-by-1 cell vector of annual two-column matrices.

Source: Yahoo Finance.

## CryptoCurrencies

Volatility-adjusted Bitcoin and Litecoin log returns, split into years 2015--2018. Its structure has the same three fields as `Banks`.

Source: CoinMetrics.

## IndexReturns2D and IndexReturns3D

European stock-index log returns for 1998. `IndexReturns2D` is 100-by-2 with columns DAX and SMI. `IndexReturns3D` is 200-by-3 with columns DAX, SMI, and CAC. Each MAT-file includes its column labels in `VariableNames`.

Source: raw data supplied by Erste Bank AG, Vienna, Austria.

See [the data conversion notes](../data/README.md) for precision and provenance.
