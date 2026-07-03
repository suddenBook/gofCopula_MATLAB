function tk = peMarginals(dimension, beta, options)
%PEMARGINALS High-accuracy 1-D marginal transforms of a power-exponential copula.
%   TK = peMarginals(DIMENSION,BETA) returns a struct holding the univariate
%   marginal generator, its log interpolant, the quantile Q1 and distribution
%   F1 of one coordinate of a DIMENSION-variate standardized power-exponential
%   law with shape BETA.
%
%   A power-exponential margin is NOT itself power-exponential (the family is
%   not closed under marginalization), so the marginal generator is the
%   Cambanis-Huang-Simons reduction of the joint generator g_d(s)=k exp(-s^b/2):
%       g_1(t) = pi^((d-1)/2)/Gamma((d-1)/2) * int_0^inf u^((d-3)/2) g_d(t+u) du
%              = 2 k exp(logCd) int_0^inf w^(d-2) exp(-(t+w^2)^b / 2) dw   (u=w^2).
%   The generator is evaluated ANALYTICALLY inside the radial quadrature (no
%   interpolation of generator values). Both the radial (w) and the marginal (x)
%   grids are sinh-spaced -- fine near the origin, exponentially coarse in the
%   tail -- so a fixed point budget resolves any shape, including the very heavy
%   tails of small beta, at bounded cost. At BETA=1 this reproduces the
%   standard-normal margin, so the copula reduces to Gaussian.
%
%   Options (all default to high resolution; runtime is not a constraint):
%     SampleSize        - sizes the tail range to cover ~1/(4n) marginal
%                         quantiles (default 1000).
%     GridPoints        - sinh-spaced x-grid nodes for F1/Q1 (default 8000).
%     QuadraturePoints  - sinh-spaced radial nodes for g_1 (default 8000).
%     MaxRadius         - x-grid upper bound; default adapts to SampleSize,BETA.

arguments
    dimension (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(dimension,2)}
    beta (1,1) double {mustBeReal, mustBeFinite, mustBePositive}
    options.SampleSize (1,1) double {mustBeInteger, mustBePositive} = 1000
    options.GridPoints (1,1) double {mustBeInteger, mustBePositive} = 8000
    options.QuadraturePoints (1,1) double {mustBeInteger, mustBePositive} = 8000
    options.MaxRadius double = []
end

d = dimension;
b = beta;
logk = peLogConstant(d, b);
logCd = (d - 1)/2 * log(pi) - gammaln((d - 1)/2);   % log[ pi^((d-1)/2)/Gamma((d-1)/2) ]

% Tail range: cover the ~1/(4n) marginal quantile with margin. T = R^(2b) ~
% Gamma(d/2b, scale 2), so the squared radius is S = T^(1/b) and x ranges to
% sqrt(S). For small beta this is huge, which the sinh grid handles at fixed cost.
if isempty(options.MaxRadius)
    nEff = max(options.SampleSize, 100);
    tQuantile = 2 * gammaincinv(1 - 1/(4*nEff), d/(2*b), "lower");
    xmax = sqrt(tQuantile^(1/b)) * 1.3 + 2;
else
    xmax = options.MaxRadius;
end

xs = sinh(linspace(0, asinh(xmax), options.GridPoints)).';   % sinh-spaced, fine near 0
g1AtX = peGeneratorQuad(xs.^2, d, b, logk, logCd, options.QuadraturePoints, xmax);

% Drop the far tail where the density underflows so F1 stays strictly monotone.
keep = g1AtX > 0 & isfinite(g1AtX);
xs = xs(keep);
f1 = g1AtX(keep);

% Marginal CDF by cumulative quadrature on the (non-uniform) grid; symmetric.
Fpositive = 0.5 + cumtrapz(xs, f1);
xFull = [-flipud(xs(2:end)); xs];
Ffull = [1 - flipud(Fpositive(2:end)); Fpositive];
[xFull, Ffull] = makeStrictlyIncreasing(xFull, Ffull);

cdfGrid = griddedInterpolant(xFull, Ffull, "makima", "nearest");
quantileGrid = griddedInterpolant(Ffull, xFull, "makima", "nearest");
logG1Grid = griddedInterpolant(xs.^2, log(f1), "makima", "nearest");

tMax = xs(end)^2;
uLo = Ffull(1);
uHi = Ffull(end);
tk = struct( ...
    "Dimension", d, "Beta", b, "LogConstant", logk, "MaxSquaredRadius", tMax, ...
    "Quantile", @(u) quantileGrid(min(max(u, uLo), uHi)), ...
    "CDF", @(x) min(max(cdfGrid(x), 0), 1), ...
    "LogMarginalGenerator", @(t) logG1Grid(min(max(t, 0), tMax)));
end

function g1 = peGeneratorQuad(t, d, b, logk, logCd, m, xmax)
% g_1(t) = 2 k exp(logCd) * int_0^inf w^(d-2) exp(-(t+w^2)^b / 2) dw, on a
% sinh-spaced w-grid (non-uniform trapezoid). The generator is exact; the
% integrand matrix is formed in chunks of t to bound memory for heavy tails.
wMax = xmax + 6;                          % w is an auxiliary radius; integrand ~0 beyond
w = sinh(linspace(0, asinh(wMax), m)).';
dw = diff(w);
trapWeight = zeros(m, 1);
trapWeight(1) = dw(1)/2;
trapWeight(end) = dw(end)/2;
trapWeight(2:end-1) = (dw(1:end-1) + dw(2:end)) / 2;
radialWeight = (w.^(d - 2)) .* trapWeight;         % w^(d-2) dw
wSquared = (w.^2).';                               % 1-by-m
scale = 2 * exp(logk + logCd);
t = t(:);
g1 = zeros(numel(t), 1);
chunk = max(1, floor(2e7 / m));                    % bound the integrand matrix size
for a = 1:chunk:numel(t)
    z = min(a + chunk - 1, numel(t));
    integrand = exp(-0.5 * (t(a:z) + wSquared).^b);   % (<=chunk)-by-m, generator exact
    g1(a:z) = scale * (integrand * radialWeight);
end
end

function logk = peLogConstant(d, beta)
%PELOGCONSTANT log of k_{d,beta} = d Gamma(d/2) / (pi^(d/2) Gamma(1+d/2b) 2^(1+d/2b)).
logk = log(d) + gammaln(d/2) - (d/2) * log(pi) ...
    - gammaln(1 + d/(2*beta)) - (1 + d/(2*beta)) * log(2);
end

function [x, y] = makeStrictlyIncreasing(x, y)
% Ensure y is strictly increasing (griddedInterpolant sample requirement) by
% nudging any non-increasing entry up by a negligible amount.
step = 1e-14 * max(1, max(abs(y)));
for k = 2:numel(y)
    if y(k) <= y(k-1)
        y(k) = y(k-1) + step;
    end
end
end
