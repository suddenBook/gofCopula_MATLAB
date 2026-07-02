function V = frailty(family, n, theta, stream)
%FRAILTY Draw Marshall--Olkin frailty variables for Archimedean copulas.
%   V = FRAILTY(FAMILY,N,THETA,STREAM) returns N positive draws whose
%   Laplace transform is the family's completely monotone generator
%   inverse psi, so that U_j = psi(E_j/V) with E_j ~ Exp(1) i.i.d. has the
%   requested copula (Marshall & Olkin, 1988):
%     clayton (theta>0):   V ~ Gamma(1/theta, 1)
%     gumbel  (theta>1):   V ~ positive stable, LT exp(-t^(1/theta))
%     frank   (theta>0):   V ~ Logarithmic(1-exp(-theta))
%     joe     (theta>1):   V ~ Sibuya(1/theta)
%     amh     (0<theta<1): V ~ Geometric(1-theta) on {1,2,...}
%   The Sibuya sampler is a port of copula:::rSibuyaR (Hofert, 2011); the
%   logarithmic sampler is Kemp's (1981) LS algorithm. All randomness is
%   drawn from STREAM ([] uses the global stream).

arguments
    family (1,1) string
    n (1,1) {mustBeInteger, mustBePositive}
    theta (1,1) double {mustBeReal, mustBeFinite}
    stream = []
end

switch family
    case "clayton"
        mustBePositive(theta);
        V = gammaSample(stream, 1/theta, n);
    case "gumbel"
        assert(theta > 1, "gofcopula:frailty:Domain", "Gumbel frailty needs theta > 1.");
        V = stableSample(stream, 1/theta, n);
    case "frank"
        mustBePositive(theta);
        V = logarithmicSample(stream, -expm1(-theta), n);
    case "joe"
        assert(theta > 1, "gofcopula:frailty:Domain", "Joe frailty needs theta > 1.");
        V = sibuyaSample(stream, 1/theta, n);
    case "amh"
        assert(theta > 0 && theta < 1, "gofcopula:frailty:Domain", ...
            "AMH frailty needs theta in (0,1).");
        V = ceil(log(uniform(stream, n, 1)) ./ log(theta));
    otherwise
        error("gofcopula:frailty:UnsupportedFamily", ...
            "No frailty distribution for the %s copula.", family);
end
end

function g = gammaSample(stream, a, n)
% Marsaglia--Tsang (2000) with the standard a<1 boost, fed from STREAM.
boost = a < 1;
shape = a + boost;
d = shape - 1/3;
c = 1 / sqrt(9 * d);
g = zeros(n, 1);
need = true(n, 1);
while any(need)
    m = nnz(need);
    x = normal(stream, m, 1);
    v = (1 + c .* x).^3;
    u = uniform(stream, m, 1);
    ok = v > 0 & log(u) < 0.5 .* x.^2 + d - d .* v + d .* log(max(v, realmin));
    index = find(need);
    g(index(ok)) = d .* v(ok);
    need(index(ok)) = false;
end
if boost
    g = g .* uniform(stream, n, 1).^(1 / a);
end
end

function v = stableSample(stream, alpha, n)
% Kanter (1975): positive alpha-stable with Laplace transform exp(-t^alpha).
angle = pi .* uniform(stream, n, 1);
w = -log(uniform(stream, n, 1));
a = sin((1 - alpha) .* angle) .* sin(alpha .* angle).^(alpha ./ (1 - alpha)) ...
    ./ sin(angle).^(1 ./ (1 - alpha));
v = (a ./ w).^((1 - alpha) ./ alpha);
end

function v = logarithmicSample(stream, p, n)
% Kemp (1981) LS: P(V=k) = p^k / (-k*log(1-p)), k = 1,2,...
u = uniform(stream, n, 1);
v = ones(n, 1);
branch = u <= p;
if any(branch)
    ub = u(branch);
    q = -expm1(uniform(stream, nnz(branch), 1) .* log1p(-p)); % 1-(1-p)^W
    value = 2 .* ones(nnz(branch), 1);
    small = ub < q.^2;
    value(small) = floor(1 + log(ub(small)) ./ log(q(small)));
    value(ub > q) = 1;
    v(branch) = value;
end
end

function v = sibuyaSample(stream, alpha, n)
% Port of copula:::rSibuyaR (Hofert 2011): exact inversion using
% P(V <= x) = 1 - 1/(x*Beta(x, 1-alpha)) for real x >= 1.
if alpha == 1
    v = ones(n, 1);
    return
end
u = uniform(stream, n, 1);
v = ones(n, 1);
big = u > alpha;
if any(big)
    ub = u(big);
    inverse = ((1 - ub) .* gamma(1 - alpha)).^(-1 / alpha);
    floored = floor(inverse);
    logBeta = betaln(floored, 1 - alpha);
    passed = 1 - exp(-log(floored) - logBeta) < ub;
    value = floored;
    value(passed) = ceil(inverse(passed));
    v(big) = value;
end
end

function u = uniform(stream, varargin)
if isempty(stream), u = rand(varargin{:}); else, u = rand(stream, varargin{:}); end
end

function z = normal(stream, varargin)
if isempty(stream), z = randn(varargin{:}); else, z = randn(stream, varargin{:}); end
end
