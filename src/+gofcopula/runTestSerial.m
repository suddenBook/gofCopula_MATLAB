function [result, serial] = runTestSerial(testName, copula, x, options)
%RUNTESTSERIAL Serial-dependence-robust parametric-bootstrap GoF test.
%   RESULT = RUNTESTSERIAL(TESTNAME, COPULA, X) runs the same goodness-of-fit
%   test as gofcopula.runTest, but corrects the parametric bootstrap so its
%   p-value stays valid on autocorrelated (time-series) data. The plain i.i.d.
%   bootstrap under-estimates the variance of the test statistic when the
%   observations are serially dependent and therefore over-rejects.
%
%   Two correction methods are available (name-value Method):
%     "decimate" (default) - thin the rows to near-independence (one row per
%         decorrelation time) and run the ordinary i.i.d. bootstrap on the
%         sub-sample. Family-agnostic; trades power for correct size.
%     "phase" - keep all rows and calibrate against coherent
%         phase-randomized surrogates (injected as BootstrapSamples). Exact
%         and full-power, but the surrogates are Gaussian by construction, so
%         this is valid ONLY for the "normal" copula and errors otherwise.
%     "multiplier" is reserved for a future dependent multiplier bootstrap and
%         currently errors.
%
%   [RESULT, SERIAL] = RUNTESTSERIAL(...) also returns a diagnostics struct:
%     method       - correction method used
%     nObserved    - rows of X
%     nThinned     - rows passed to the bootstrap (= nObserved for "phase")
%     thinInterval - decimation interval (1 for "phase")
%     offset       - starting row of the decimation (always 1)
%     iatPerColumn - per-column integrated autocorrelation time (NaN when not
%                    estimated, e.g. ThinInterval supplied or Method="phase")
%     maxLag       - ACF lag used for the interval estimate (NaN if not used)
%     keepIndices  - retained row indices into X
%
%   RESULT has the same type as gofcopula.runTest's output
%   (gofcopula.GofResult); RESULT.Tests carries the serial-robust p-value.
%
%   All gofcopula.runTest name-value options are accepted and forwarded.
%   Additional serial options:
%     ThinInterval - fixed decimation interval ("decimate"); empty (default)
%                    estimates it from the rank ACF.
%     MaxLag       - maximum ACF lag for the automatic interval estimate.
%     Method       - "decimate" (default), "phase", or "multiplier".
%
%   Example:
%     % Family-agnostic (any copula), thinned:
%     [r, s] = gofcopula.runTestSerial("gofCvM", "clayton", X, M=999, Seed=1);
%     % Gaussian test, full-power surrogate null:
%     [r, s] = gofcopula.runTestSerial("gofCvM", "normal", X, ...
%                  Method="phase", M=999, Seed=1);
%
%   See also gofcopula.runTest, gofcopula.gofCvM.

arguments
    testName {mustBeTextScalar}
    copula {mustBeTextScalar}
    x {mustBeFloat,mustBeReal,mustBeFinite,mustBeMatrix}
    options.Param {mustBeNumeric,mustBeReal,mustBeFinite} = 0.5
    options.ParamEst (1,1) logical = true
    options.DF (1,1) {mustBeNumeric,mustBeReal,mustBeFinite,mustBePositive} = 4
    options.DFEst (1,1) logical = true
    options.Margins = "ranks"
    options.Flip (1,1) {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger} = 0
    options.M (1,1) {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger,mustBeNonnegative} = 1000
    options.MJ (1,1) {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger,mustBePositive} = 100
    options.Dispersion {mustBeTextScalar} = "exchangeable"
    options.BlockSize (1,1) {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger,mustBePositive} = 1
    options.KernelScale (1,1) {mustBeNumeric,mustBeReal,mustBeFinite,mustBePositive} = 0.5
    options.IntegrationNodes (1,1) {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger,mustBePositive} = 12
    options.Lower {mustBeNumeric,mustBeReal} = []
    options.Upper {mustBeNumeric,mustBeReal} = []
    options.Seed {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger} = []
    options.Processes (1,1) {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger,mustBePositive} = 1
    options.NumericMode {mustBeTextScalar} = "corrected"
    options.CustomTest = []
    options.ModelSample {mustBeNumeric,mustBeReal,mustBeFinite} = []
    options.BootstrapSamples {mustBeNumeric,mustBeReal,mustBeFinite} = []
    options.ThinInterval {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger,mustBePositive} = []
    options.MaxLag {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger,mustBePositive} = []
    options.Method {mustBeTextScalar} = "decimate"
end

methodName = string(validatestring(options.Method, ["decimate", "phase", "multiplier"]));
switch methodName
    case "decimate"
        [result, serial] = decimateAndRun(testName, copula, x, options);
    case "phase"
        [result, serial] = phaseAndRun(testName, copula, x, options);
    otherwise
        error("gofcopula:Serial:NotImplemented", ...
            "Method '%s' is not implemented; use 'decimate' or 'phase'.", methodName);
end
end

function [result, serial] = decimateAndRun(testName, copula, x, options)
% Thin to near-independence, then run the ordinary i.i.d. bootstrap.
if isempty(options.ThinInterval)
    [interval, dInfo] = gofcopula.internal.resampling.decorrelationLength( ...
        x, MaxLag=options.MaxLag);
    iat = dInfo.iatPerColumn;
    maxLagUsed = dInfo.maxLag;
else
    interval = double(options.ThinInterval);
    iat = NaN(1, size(x, 2));
    maxLagUsed = NaN;
end

% Decimate the original rows; runTest applies its own Margins transform, so
% pass raw rows (not ranks) to avoid double-ranking.
[xThin, keepIdx] = gofcopula.internal.resampling.thinToIndependence(x, interval, 1);

% runTest requires BlockSize to divide the number of rows (PIOS tests).
if options.BlockSize > 1
    nThin = size(xThin, 1);
    usable = nThin - mod(nThin, options.BlockSize);
    xThin = xThin(1:usable, :);
    keepIdx = keepIdx(1:usable);
end

forward = rmfield(options, ["ThinInterval", "MaxLag", "Method"]);
args = namedargs2cell(forward);
result = gofcopula.runTest(testName, copula, xThin, args{:});

serial = struct( ...
    "method", "decimate", ...
    "nObserved", size(x, 1), ...
    "nThinned", size(xThin, 1), ...
    "thinInterval", interval, ...
    "offset", 1, ...
    "iatPerColumn", iat, ...
    "maxLag", maxLagUsed, ...
    "keepIndices", keepIdx);
end

function [result, serial] = phaseAndRun(testName, copula, x, options)
% Keep all rows; calibrate against coherent phase-randomized surrogates.
family = lower(string(copula));
if family ~= "normal" && family ~= "gaussian"
    error("gofcopula:Serial:PhaseRequiresGaussian", ...
        "Method 'phase' builds a Gaussian-copula null and is valid only for " + ...
        "the 'normal' copula; got '%s'. Use Method='decimate' for other families.", copula);
end
if ~(isscalar(string(options.Margins)) && string(options.Margins) == "ranks")
    error("gofcopula:Serial:PhaseRequiresRanks", ...
        "Method 'phase' generates raw surrogates that must be re-ranked; " + ...
        "use Margins='ranks' (the default).");
end
if ~isempty(options.BootstrapSamples)
    error("gofcopula:Serial:PhaseConflict", ...
        "Method 'phase' generates its own BootstrapSamples; do not also supply them.");
end

[n, d] = size(x);
forward = rmfield(options, ["ThinInterval", "MaxLag", "Method"]);
if options.M >= 1
    surrStream = seededStream(options.Seed);
    samples = zeros(n, d, options.M);
    for b = 1:options.M
        samples(:, :, b) = gofcopula.internal.resampling.phaseSurrogate(x, surrStream);
    end
    forward.BootstrapSamples = samples;
end
args = namedargs2cell(forward);
result = gofcopula.runTest(testName, copula, x, args{:});

serial = struct( ...
    "method", "phase", ...
    "nObserved", n, ...
    "nThinned", n, ...
    "thinInterval", 1, ...
    "offset", 1, ...
    "iatPerColumn", NaN(1, d), ...
    "maxLag", NaN, ...
    "keepIndices", (1:n).');
end

function stream = seededStream(seed)
% Deterministic surrogate stream when a seed is supplied, random otherwise.
if isempty(seed)
    seed = randi(intmax("int32"));
end
stream = RandStream("Threefry", "Seed", double(seed));
end
