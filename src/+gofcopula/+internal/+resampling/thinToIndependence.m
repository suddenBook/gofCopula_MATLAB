function [xThin, keepIdx] = thinToIndependence(x, interval, offset)
%THINTOINDEPENDENCE Decimate matrix rows by a fixed interval.
%   XTHIN = THINTOINDEPENDENCE(X, INTERVAL) keeps every INTERVAL-th row of X,
%   starting from the first row, so consecutive retained rows are INTERVAL
%   apart. Paired with gofcopula.internal.resampling.decorrelationLength this
%   produces an approximately independent sub-sample for the parametric
%   bootstrap. The same rows are kept for every column, so the
%   cross-sectional dependence (the copula) is preserved.
%
%   XTHIN = THINTOINDEPENDENCE(X, INTERVAL, OFFSET) starts from row OFFSET
%   instead of row 1 (default 1). INTERVAL = 1 returns X unchanged.
%
%   [XTHIN, KEEPIDX] = THINTOINDEPENDENCE(...) also returns the retained row
%   indices as a column vector.

arguments
    x {mustBeFloat,mustBeReal,mustBeFinite,mustBeMatrix}
    interval (1,1) {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger,mustBePositive}
    offset (1,1) {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger,mustBePositive} = 1
end
n = size(x, 1);
if offset > n
    error("gofcopula:Serial:Offset", ...
        "Offset (%d) exceeds the number of rows (%d).", offset, n);
end
keepIdx = (offset:interval:n).';
xThin = x(keepIdx, :);
end
