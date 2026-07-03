function [fitted, method] = estimateModel(model, u, options)
%ESTIMATEMODEL Estimate copula parameters by pseudo-likelihood, then iTau.
%   The pseudo-likelihood step maximizes this implementation's log-density for
%   every family and dimension, mirroring R's fitCopula(method="mpl"):
%   scalar parameters use golden-section search on the family's parameter
%   domain intersected with the user's Lower/Upper bounds; Student and
%   t-EV degrees of freedom are profiled jointly when requested; only the
%   unstructured elliptical case delegates to copulafit (the unconstrained
%   maximum likelihood) when no bounds are given. On failure, estimation
%   falls back to exact inversion of Kendall's tau (R: method="itau").
%
%   Domain guards mirror gofCopula 0.4-3: a negative Clayton estimate, a
%   negative Frank estimate above dimension two, and an estimate exactly
%   at the family's parameter-domain boundary are errors.

arguments
    model (1,1) gofcopula.CopulaModel
    u {mustBeFloat,mustBeReal,mustBeFinite,mustBeMatrix}
    options.Lower {mustBeNumeric,mustBeReal} = []
    options.Upper {mustBeNumeric,mustBeReal} = []
end

wantDF = ismember(model.Family, ["t","tev","powerexp"]) && model.EstimateDegreesOfFreedom;
if ~model.EstimateTheta && ~wantDF
    fitted = model; method = "fixed"; return
end
u = min(max(u, eps), 1 - eps);

try
    [theta, df] = pseudoLikelihood(model, u, options.Lower, options.Upper);
    if any(~isfinite(theta)) || ~isfinite(df)
        error("gofcopula:Estimation:Nonfinite", "Nonfinite estimate.");
    end
    method = "mpl";
catch reason
    [theta, df] = gofcopula.internal.estimation.invertTau(model, u);
    method = "itau";
    warning("gofcopula:Estimation:Fallback", ...
        "Pseudo-likelihood estimation failed for the %s copula (%s) " + ...
        "Inversion of Kendall's tau was used.", model.Family, reason.message);
end
domainGuards(model.Family, size(u,2), theta);
fitted = model.withFit(theta, df);
end

function [theta, df] = pseudoLikelihood(model, u, lower, upper)
d = size(u,2);
family = model.Family;
df = model.DegreesOfFreedom;
wantTheta = model.EstimateTheta;
wantDF = ismember(family, ["t","tev"]) && model.EstimateDegreesOfFreedom;
dfBounds = [0.1, 1000];

if family == "powerexp"
    [theta, df] = powerExponential(model, u, lower, upper);
    return
end

if ismember(family, ["normal","t"]) && d > 2 && model.Dispersion == "unstructured"
    [theta, df] = unstructuredElliptical(model, u, lower, upper, wantTheta, wantDF);
    return
end

[lo, hi] = scalarBounds(family, d, lower, upper);
if ismember(family, ["t","tev"])
    if wantTheta && wantDF
        start = [middle(model.Theta(1), lo, hi), middle(df, dfBounds(1), dfBounds(2))];
        p = maximizeVector(@(p) negLog(model, u, p(1), p(2)), start, ...
            [lo, dfBounds(1)], [hi, dfBounds(2)]);
        theta = p(1); df = p(2);
    elseif wantTheta
        theta = maximizeScalar(@(t) negLog(model, u, t, df), lo, hi);
    else
        theta = model.Theta;
        df = maximizeScalar(@(nu) negLog(model, u, theta, nu), dfBounds(1), dfBounds(2));
    end
elseif family == "plackett"
    % Optimize on the log scale: the domain spans several decades.
    logTheta = maximizeScalar(@(lt) negLog(model, u, exp(lt), df), log(lo), log(hi));
    theta = exp(logTheta);
else
    theta = maximizeScalar(@(t) negLog(model, u, t, df), lo, hi);
end
end

function [theta, df] = unstructuredElliptical(model, u, lower, upper, wantTheta, wantDF)
d = size(u,2);
df = model.DegreesOfFreedom;
if ~wantTheta
    theta = model.Theta;
    df = maximizeScalar(@(nu) negLog(model, u, theta, nu), 0.1, 1000);
    return
end
if isempty(lower) && isempty(upper)
    if model.Family == "normal"
        rho = copulafit("Gaussian", u);
    else
        [rho, nu] = copulafit("t", u);
        if wantDF, df = nu; end
    end
    theta = rho(tril(true(d), -1)).';
else
    theta = boundedUnstructured(model, u, lower, upper);
    if wantDF
        df = maximizeScalar(@(nu) negLog(model, u, theta, nu), 0.1, 1000);
    end
end
end

function theta = boundedUnstructured(model, u, lower, upper)
d = size(u,2);
m = d * (d - 1) / 2;
lo = -0.999999 * ones(1, m);
hi = 0.999999 * ones(1, m);
if ~isempty(lower), lo = max(lo, lower(1)); end
if ~isempty(upper), hi = min(hi, upper(1)); end
tauMatrix = corr(u, "Type", "Kendall", "Rows", "complete");
start = sin(pi .* tauMatrix(tril(true(d), -1)).' ./ 2);
start = min(max(start, lo + 1e-6), hi - 1e-6);
theta = maximizeVector(@(p) negLog(model, u, p, model.DegreesOfFreedom), start, lo, hi);
end

function [theta, df] = powerExponential(model, u, lower, upper)
% Power-exponential copula: correlation R by inversion of Kendall's tau
% (elliptical, shape-independent), then shape beta by 1-D pseudo-likelihood
% with R held fixed. df carries beta. Mirrors the two-step estimator in the
% Power-Exponential validation pipeline.
d = size(u, 2);
if model.EstimateTheta
    tauMatrix = corr(u, "Type", "Kendall", "Rows", "complete");
    if d == 2
        rho = sin(pi * tauMatrix(2,1) / 2);
        if ~isempty(lower), rho = max(rho, lower(1)); end
        if ~isempty(upper), rho = min(rho, upper(1)); end
        theta = rho;
    elseif model.Dispersion == "unstructured"
        R = sin(pi .* tauMatrix ./ 2);
        R(1:d+1:end) = 1;
        [~, flag] = chol(R);
        if flag ~= 0
            R = nearestCorrelation(R);   % project the tau matrix to a valid correlation
        end
        theta = R(tril(true(d), -1)).';
    else
        % Exchangeable: one correlation from the mean pairwise Kendall tau.
        meanTau = mean(tauMatrix(tril(true(d), -1)));
        rho = sin(pi * meanTau / 2);
        if ~isempty(lower), rho = max(rho, lower(1)); end
        if ~isempty(upper), rho = min(rho, upper(1)); end
        theta = rho;
    end
else
    theta = model.Theta;
end
df = model.DegreesOfFreedom;
if model.EstimateDegreesOfFreedom
    betaBounds = [0.2, 5];
    df = maximizeShape(@(b) negLog(model, u, theta, b), betaBounds(1), betaBounds(2));
end
end

function R = nearestCorrelation(R)
% Nearest correlation matrix by clipping negative eigenvalues and rescaling to
% a unit diagonal (sufficient for a Kendall-tau starting value).
R = (R + R.') / 2;
[V, D] = eig(R);
D = max(diag(D), 1e-8);
R = V * diag(D) * V.';
s = sqrt(diag(R));
R = R ./ (s * s.');
R = (R + R.') / 2;
R(1:size(R,1)+1:end) = 1;
end

function b = maximizeShape(negative, lo, hi)
% Golden-section search for the shape parameter (tight tolerance; runtime is
% not a constraint and a precise optimum is preferred).
b = fminbnd(negative, lo, hi, optimset("TolX", 1e-5, "Display", "off"));
end

function value = negLog(model, u, theta, df)
try
    logDensity = gofcopula.copulaPDF(model.Family, u, theta, DF=df, ...
        Dispersion=model.Dispersion, Rotation=model.Rotation, Log=true);
    value = -sum(logDensity);
    if ~isfinite(value)
        value = 1e10;
    end
catch
    value = 1e10; % invalid parameter (e.g. non-positive-definite correlation)
end
end

function x = maximizeScalar(negative, lo, hi)
if ~(lo < hi)
    error("gofcopula:Estimation:EmptyBounds", ...
        "The lower bound must lie strictly below the upper bound.");
end
x = fminbnd(negative, lo, hi, optimset("TolX", 1e-12, "Display", "off"));
end

function x = maximizeVector(negative, start, lo, hi)
% Central differences and explicit scaling: several objectives (Student
% and t-EV profiles) are shallow in the df direction, and the extreme
% value densities carry small internal finite-difference noise, so the
% default forward-difference gradients can stall at the start point.
x = fmincon(negative, start, [], [], [], [], lo, hi, [], ...
    optimoptions("fmincon", Display="off", Algorithm="interior-point", ...
    OptimalityTolerance=1e-9, StepTolerance=1e-12, ...
    FiniteDifferenceType="central", TypicalX=max(abs(start), 0.25)));
end

function [lo, hi] = scalarBounds(family, d, lower, upper)
switch family
    case {"normal", "t"}
        if d == 2
            lo = -0.999999;
        else
            lo = -1 / (d - 1) + 1e-6; % exchangeable positive definiteness
        end
        hi = 0.999999;
    case "tev"
        lo = -0.999999; hi = 0.999999;
    case "clayton"
        if d == 2, lo = -0.999999; else, lo = 1e-8; end
        hi = 100;
    case {"gumbel", "joe"}
        lo = 1 + 1e-10; hi = 100;
    case "frank"
        if d == 2, lo = -100; else, lo = 1e-8; end
        hi = 100;
    case "amh"
        if d == 2, lo = -0.999999; else, lo = 1e-8; end
        hi = 0.999999;
    case {"galambos", "huslerreiss"}
        lo = 1e-8; hi = 100;
    case "tawn"
        lo = 1e-8; hi = 1 - 1e-8;
    case "fgm"
        lo = -1 + 1e-8; hi = 1 - 1e-8;
    case "plackett"
        lo = 1e-4; hi = 1e4;
    otherwise
        lo = -100; hi = 100;
end
if ~isempty(lower), lo = max(lo, lower(1)); end
if ~isempty(upper), hi = min(hi, upper(1)); end
end

function value = middle(value, lo, hi)
value = min(max(value, lo + sqrt(eps)), hi - sqrt(eps));
end

function domainGuards(family, d, theta)
% R parity (internal_param_est.R:156-169 and 312-321).
if family == "clayton" && any(theta < 0)
    error("gofcopula:Estimation:InvalidClayton", ...
        "The dependence parameter is negative; the Clayton copula is not " + ...
        "an appropriate model for this dataset.");
end
if family == "frank" && d > 2 && any(theta < 0)
    error("gofcopula:Estimation:InvalidFrank", ...
        "A negative Frank parameter is only possible in dimension two.");
end
[lo, hi] = intrinsicDomain(family, d);
if any(theta == lo) || (isfinite(hi) && any(theta == hi))
    error("gofcopula:Estimation:BoundaryEstimate", ...
        "The estimated %s parameter lies at its domain boundary, which " + ...
        "produces unstable goodness-of-fit results. Consider Lower/Upper " + ...
        "bounds or a different copula.", family);
end
end

function [lo, hi] = intrinsicDomain(family, d)
switch family
    case {"normal", "t", "tev", "powerexp"}
        lo = -1; hi = 1;
    case "clayton"
        if d == 2, lo = -1; else, lo = 0; end
        hi = Inf;
    case {"gumbel", "joe"}
        lo = 1; hi = Inf;
    case "frank"
        lo = -Inf; hi = Inf;
    case "amh"
        lo = -1; hi = 1;
    case {"galambos", "huslerreiss"}
        lo = 0; hi = Inf;
    case "tawn"
        lo = 0; hi = 1;
    case "fgm"
        lo = -1; hi = 1;
    case "plackett"
        lo = 0; hi = Inf;
    otherwise
        lo = -Inf; hi = Inf;
end
end
