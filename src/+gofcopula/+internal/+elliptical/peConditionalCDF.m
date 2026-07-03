function C = peConditionalCDF(prefix, Rfull, beta, ambientDim)
%PECONDITIONALCDF Rosenblatt conditional CDF for the power-exponential copula.
%   C = peConditionalCDF(PREFIX,RFULL,BETA,AMBIENTDIM) returns
%   P(U_j <= PREFIX(:,j) | U_1,...,U_{j-1}) for the PE copula, where PREFIX is
%   n-by-j, RFULL is the d-by-d correlation matrix and AMBIENTDIM = d.
%
%   Working in the elliptical latent space z = Q1(u) (Q1 the true d-dimensional
%   PE marginal quantile), the conditional law of z_j given z_{<j} is the first
%   coordinate of an elliptically contoured EC_{d-j+1} with generator
%   g_c(t) = exp(-(t+q)^beta/2), q = z_{<j}' R_{<j}^{-1} z_{<j} (Gomez et al.
%   1998, Prop 5.1(ii)). With standardized residual r = (z_j - mu)/sqrt(v),
%       C = int_{-inf}^{r} g_{c,1}(s^2) ds  /  int_{-inf}^{inf} g_{c,1}(s^2) ds,
%   where g_{c,1} marginalizes g_c from dimension p=d-j+1 to 1. The ratio is
%   tabulated on a (q,r) grid and interpolated per observation.

arguments
    prefix (:,:) double
    Rfull (:,:) double
    beta (1,1) double {mustBeReal, mustBeFinite, mustBePositive}
    ambientDim (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(ambientDim,2)}
end

[n, j] = size(prefix);
tk = gofcopula.internal.elliptical.peMarginals(ambientDim, beta, SampleSize=max(n,2000));
z = tk.Quantile(prefix);                       % latent coordinates, n-by-j

Rj = Rfull(1:j, 1:j);
Rprev = Rj(1:j-1, 1:j-1);
cross = Rj(j, 1:j-1);                           % 1-by-(j-1)
Lprev = chol(Rprev, "lower");
zprev = z(:, 1:j-1);

weights = cross / Rprev;                        % elliptical regression weights
mu = zprev * weights.';                         % conditional mean, n-by-1
condVar = max(Rj(j,j) - weights * cross.', realmin);
r = (z(:,j) - mu) / sqrt(condVar);              % standardized residual, n-by-1
q = sum((Lprev \ zprev.').^2, 1).';             % conditioning Mahalanobis, n-by-1
p = ambientDim - j + 1;                          % conditional elliptical dimension

C = conditionalRatio(r, q, beta, p);
C = min(max(C, 0), 1);
end

function H = conditionalRatio(r, q, beta, p)
% Tabulate H(r;q) = cumulative(g_{c,1}(r^2;q)) / total on a fine grid, then
% shape-preserving (makima) interpolate at (r,q). Z_j depends on q through z_1^2,
% so a fine, smooth q-interpolation is needed to avoid spurious dependence.
rMax = max(9, max(abs(r)) * 1.15 + 1);
rGrid = linspace(-rMax, rMax, 2401).';
qGrid = linspace(0, max(q) * 1.02 + 0.5, 600);
Htab = zeros(numel(rGrid), numel(qGrid));
for iq = 1:numel(qGrid)
    density = conditionalGenerator(rGrid.^2, qGrid(iq), beta, p);   % along rGrid
    cumulative = cumtrapz(rGrid, density);
    Htab(:, iq) = cumulative / cumulative(end);
end
interpolant = griddedInterpolant({rGrid, qGrid}, Htab, "makima", "nearest");
H = interpolant(r, min(max(q, 0), qGrid(end)));
end

function g = conditionalGenerator(s, q, beta, p)
% 1-D conditional margin generator g_{c,1}(s;q): reduce exp(-(t+q)^b/2) from
% dimension p to 1 at squared-residual s (unnormalized; the ratio cancels k).
if p == 1
    g = exp(-0.5 * (s + q).^beta);
else
    logCd = (p - 1)/2 * log(pi) - gammaln((p - 1)/2);
    m = 4000;
    wMax = sqrt(max(s) + q) + 12;
    w = linspace(0, wMax, m).';
    hw = w(2) - w(1);
    wt = w.^(p - 2); wt(1) = 0.5*wt(1); wt(end) = 0.5*wt(end); wt = wt * hw;
    w2 = (w.^2).';
    s = s(:);
    g = zeros(numel(s), 1);
    chunk = max(1, floor(2e7 / m));
    for a = 1:chunk:numel(s)
        z = min(a + chunk - 1, numel(s));
        g(a:z) = exp(-0.5 * (s(a:z) + w2 + q).^beta) * wt;
    end
    g = 2 * exp(logCd) * g;
end
end
