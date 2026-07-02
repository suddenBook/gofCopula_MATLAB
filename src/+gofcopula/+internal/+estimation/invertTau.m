function [theta, df] = invertTau(model, u)
%INVERTTAU Estimate copula parameters by exact inversion of Kendall's tau.
%   Fallback estimator mirroring R's fitCopula(method="itau"). Families
%   with closed-form tau(theta) use it directly; the remaining families
%   invert exact integral representations:
%     Archimedean:   tau = 1 + 4*int_0^1 phi(t)/phi'(t) dt
%     Extreme value: tau = int_0^1 t(1-t)*A''(t)/A(t) dt   (Ghoudi et al.)
%     Plackett:      tau = 1 - 4*int int dC/du * dC/dv du dv
%   No estimate is ever produced by simulation.

arguments
    model (1,1) gofcopula.CopulaModel
    u {mustBeFloat,mustBeReal,mustBeFinite,mustBeMatrix}
end

d = size(u,2);
df = model.DegreesOfFreedom;
family = model.Family;
tauMatrix = corr(u, "Type", "Kendall", "Rows", "complete");
pairTaus = tauMatrix(tril(true(d), -1));
tau = mean(pairTaus);

switch family
    case {"normal", "t"}
        if d > 2 && model.Dispersion == "unstructured"
            theta = sin(pi .* pairTaus(:).' ./ 2);
        else
            theta = sin(pi * tau / 2);
        end
    case "tev"
        theta = sin(pi * tau / 2);
    case "clayton"
        theta = 2 * tau / max(1 - tau, eps);
    case "gumbel"
        theta = 1 / max(1 - tau, eps);
    case "frank"
        theta = copulaparam("Frank", tau);
    case "fgm"
        % tau = 2*theta/9 exactly.
        theta = max(-1, min(1, 4.5 * tau));
    case "amh"
        theta = invertAMH(tau);
    case {"joe"}
        theta = invertArchimedean(family, tau);
    case "plackett"
        theta = invertPlackett(tau);
    case {"galambos", "huslerreiss", "tawn"}
        theta = invertExtremeValue(family, tau, df);
    otherwise
        error("gofcopula:Estimation:UnsupportedTauInversion", ...
            "No Kendall inversion is implemented for the %s copula.", family);
end
end

function theta = invertAMH(tau)
% Exact tau(theta) = 1 - 2*((1-theta)^2*log(1-theta) + theta)/(3*theta^2),
% attainable on [(5-8*log(2))/3, 1/3).
tauMin = (5 - 8*log(2)) / 3;
tauMax = 1/3;
clamped = min(max(tau, tauMin + 1e-6), tauMax - 1e-6);
if clamped ~= tau
    warning("gofcopula:Estimation:TauOutOfRange", ...
        "Kendall's tau %.4f clamped into the attainable AMH range.", tau);
end
theta = fzero(@(t) amhTau(t) - clamped, [-1 + 1e-9, 1 - 1e-9]);
end

function value = amhTau(theta)
if abs(theta) < 1e-8
    value = 2 * theta / 9; % series limit of the exact expression at zero
else
    value = 1 - 2 * ((1 - theta)^2 * log1p(-theta) + theta) / (3 * theta^2);
end
end

function theta = invertArchimedean(family, tau)
tau = min(max(tau, 1e-6), 0.995);
objective = @(theta) archimedeanTau(family, theta) - tau;
lo = 1 + 1e-6;
hi = 50;
while objective(hi) < 0 && hi < 1e6
    hi = hi * 10;
end
theta = fzero(objective, [lo, hi]);
end

function value = archimedeanTau(family, theta)
integrand = @(t) gofcopula.internal.copulas.archimedean("phi", family, t, theta) ...
    ./ (-exp(gofcopula.internal.copulas.archimedean("logphiprime", family, t, theta)));
value = 1 + 4 * integral(integrand, 0, 1, AbsTol=1e-10, RelTol=1e-8);
end

function theta = invertPlackett(tau)
tau = min(max(tau, -0.995), 0.995);
objective = @(logTheta) plackettTau(exp(logTheta)) - tau;
lo = log(1e-4);
hi = log(1e4);
while objective(lo) > 0 && lo > log(1e-8), lo = lo - log(10); end
while objective(hi) < 0 && hi < log(1e8), hi = hi + log(10); end
theta = exp(fzero(objective, [lo, hi]));
end

function value = plackettTau(theta)
k = 40;
[x, w] = gaussLegendre01(k);
[U, V] = ndgrid(x, x);
weights = w(:) * w(:).';
Cu = reshape(gofcopula.internal.copulas.conditionalCDF( ...
    "plackett", [U(:), V(:)], theta, 4, "unstructured"), k, k);
Cv = reshape(gofcopula.internal.copulas.conditionalCDF( ...
    "plackett", [V(:), U(:)], theta, 4, "unstructured"), k, k);
value = 1 - 4 * sum(weights .* Cu .* Cv, "all");
end

function theta = invertExtremeValue(family, tau, df)
switch family
    case "tawn"
        bracket = [1e-6, 1 - 1e-6];
    otherwise
        bracket = [1e-6, 50];
end
objective = @(theta) extremeValueTau(family, theta, df) - tau;
low = objective(bracket(1));
high = objective(bracket(2));
if low > 0
    warning("gofcopula:Estimation:TauOutOfRange", ...
        "Kendall's tau %.4f below the attainable %s range.", tau, family);
    theta = bracket(1);
elseif high < 0
    warning("gofcopula:Estimation:TauOutOfRange", ...
        "Kendall's tau %.4f above the attainable %s range.", tau, family);
    theta = bracket(2);
else
    theta = fzero(objective, bracket);
end
end

function value = extremeValueTau(family, theta, df)
value = integral(@(w) integrandAt(w), 1e-8, 1 - 1e-8, AbsTol=1e-9, RelTol=1e-7);
    function y = integrandAt(w)
        [A, ~, A2] = gofcopula.internal.copulas.evDependence(family, w, theta, df);
        y = w .* (1 - w) .* A2 ./ A;
    end
end

function [nodes, weights] = gaussLegendre01(k)
index = (1:k-1).';
beta = index ./ sqrt(4 * index.^2 - 1);
[V, D] = eig(diag(beta, 1) + diag(beta, -1));
[raw, order] = sort(diag(D));
V = V(:, order);
nodes = (raw + 1) / 2;
weights = (2 * V(1, :).^2).' / 2;
end
