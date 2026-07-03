function [interval, info] = decorrelationLength(x, options)
%DECORRELATIONLENGTH Row-decimation interval for near-independent thinning.
%   INTERVAL = DECORRELATIONLENGTH(X) estimates how many rows apart the rows
%   of X must be taken for the retained rows to be approximately independent.
%   It is the tool behind gofcopula.runTestSerial: decimating X by INTERVAL
%   yields a sub-sample on which the i.i.d. parametric bootstrap is valid.
%
%   The decorrelation length is the integrated autocorrelation time (IAT)
%       tau = 1 + 2 * sum_{k>=1} rho(k),
%   computed on the COLUMN RANKS of X (pseudo-observations). Working on ranks
%   makes the estimate invariant to monotone marginal transforms, robust to
%   heavy tails, and matched to the goodness-of-fit statistics in this
%   package, which are all functions of the ranks. The autocovariance is
%   formed by FFT (Wiener-Khinchin) in base MATLAB, so no additional toolbox
%   is required. The sum uses Geyer's initial-positive-sequence truncation
%   (stop at the first non-positive lag) to tame the noisy negative tail, and
%   a one-sided lag-1 significance gate so that genuinely independent data
%   returns tau = 1 rather than being over-thinned by autocorrelation noise.
%
%   [INTERVAL, INFO] = DECORRELATIONLENGTH(...) also returns diagnostics:
%     iatPerColumn - 1-by-d IAT per column (before the retained-rows cap)
%     interval     - the returned integer thinning interval
%     nEffective   - retained rows floor((n-1)/interval)+1
%     maxLag       - maximum ACF lag used
%     minRetained  - the retained-rows floor applied
%     method       - "rank-acf-iat"
%
%   Name-value options:
%     MaxLag      - maximum ACF lag (default auto: min(n-1,max(20,ceil(10*log10(n)))))
%     MinRetained - floor on retained rows; caps INTERVAL so the sub-sample
%                   stays usable on short series (default 50)
%
%   X must be a real, finite [n x d] matrix with n >= 2. Columns are treated
%   independently and aggregated by the maximum IAT (the most conservative
%   choice); the same INTERVAL is applied to every column so the
%   cross-sectional copula is preserved.

arguments
    x {mustBeFloat,mustBeReal,mustBeFinite,mustBeMatrix}
    options.MaxLag {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger,mustBePositive} = []
    options.MinRetained (1,1) {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger,mustBePositive} = 50
end
[n, d] = size(x);
if n < 2
    error("gofcopula:Serial:TooFewRows", ...
        "DecorrelationLength needs at least two rows; got %d.", n);
end

% Column ranks (pseudo-observations) make the estimate marginal-free.
r = zeros(n, d);
for j = 1:d
    r(:, j) = tiedrank(x(:, j));
end

if isempty(options.MaxLag)
    maxLag = min(n - 1, max(20, ceil(10 * log10(n))));
else
    maxLag = min(double(options.MaxLag), n - 1);
end

% FFT (Wiener-Khinchin) autocovariance, zero-padded past 2n-1 to avoid
% circular wrap-around. Biased, lag-0-normalized ACF: the implicit 1/n
% cancels between the numerator and the lag-0 denominator.
xc = r - mean(r, 1);
nfft = 2^nextpow2(2 * n - 1);
xf = fft(xc, nfft);
acov = real(ifft(xf .* conj(xf)));   % nfft-by-d, lag 0 in row 1
c0 = acov(1, :);                      % 1-by-d lag-0 energy

% Integrated autocorrelation time per column, with the lag-1 gate.
tau = ones(1, d);
gate = 2 / sqrt(n);
for j = 1:d
    if c0(j) <= 0
        continue   % constant/degenerate column: leave tau(j) = 1
    end
    rho = acov(2:maxLag + 1, j) / c0(j);
    if rho(1) <= gate
        continue   % no detectable positive lag-1 dependence: tau(j) = 1
    end
    partial = 0;
    k = 1;
    while k <= maxLag && rho(k) > 0
        partial = partial + rho(k);
        k = k + 1;
    end
    tau(j) = 1 + 2 * partial;
end

% Aggregate (most conservative), map to an integer interval, and cap so the
% retained sub-sample keeps at least MinRetained rows when feasible.
mRaw = max(1, round(max(tau)));
if options.MinRetained <= 1
    mCap = n;
else
    mCap = max(1, floor((n - 1) / (options.MinRetained - 1)));
end
interval = min(mRaw, mCap);
nEff = floor((n - 1) / interval) + 1;

info = struct( ...
    "iatPerColumn", tau, ...
    "interval", interval, ...
    "nEffective", nEff, ...
    "maxLag", maxLag, ...
    "minRetained", options.MinRetained, ...
    "method", "rank-acf-iat");
end
