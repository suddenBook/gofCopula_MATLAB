%% Serial-dependence-robust goodness-of-fit (synthetic demonstration)
% gofcopula.runTestSerial decimates autocorrelated data to near-independence
% before the parametric bootstrap, so the goodness-of-fit p-value keeps its
% nominal size. The plain i.i.d. bootstrap (gofcopula.runTest) over-rejects
% under serial dependence. This script demonstrates both facts on synthetic
% data; see docs/SerialDependence.md for the method.
%
% Requires the gofcopula package on the path; no external data is used.
% Increase T and M for tighter (publication-quality) size estimates.

alpha = 0.05;

%% Type-I error under AR(1) serial dependence (true Gaussian-copula null)
% Y has an exact Gaussian(R) cross-sectional copula, but every column is an
% AR(1) process, so the rows are serially dependent. Under this true null a
% correctly sized test should reject at about alpha. The i.i.d. bootstrap
% instead over-rejects (increasingly so as the autocorrelation phi grows),
% while the decimated test stays near alpha.
d = 2;
n = 800;
rho = 0.5;
R = rho * ones(d) + (1 - rho) * eye(d);
T = 60;
M = 99;
phis = [0.3 0.6 0.9];

fprintf("AR(1) Gaussian-copula null (n=%d, rho=%.2f, T=%d, M=%d, alpha=%.2f)\n", ...
    n, rho, T, M, alpha);
fprintf("%6s %12s %12s %14s\n", "phi", "iid size", "serial size", "mean interval");
for phi = phis
    rejIID = 0;
    rejSerial = 0;
    intervals = zeros(T, 1);
    for t = 1:T
        Y = arGaussianNull(n, d, phi, R, 70000 + t);
        iidResult = gofcopula.runTest("gofCvM", "normal", Y, M=M, Seed=80000 + t);
        [serialResult, info] = gofcopula.runTestSerial("gofCvM", "normal", Y, ...
            M=M, Seed=80000 + t);
        rejIID = rejIID + double(iidResult.Tests.PValue < alpha);
        rejSerial = rejSerial + double(serialResult.Tests.PValue < alpha);
        intervals(t) = info.thinInterval;
    end
    fprintf("%6.2f %12.3f %12.3f %14.1f\n", ...
        phi, rejIID / T, rejSerial / T, mean(intervals));
end
fprintf("Expected: iid size well above %.2f and rising with phi; serial size ~ %.2f.\n", ...
    alpha, alpha);

%% Phase method: full-power Gaussian null (keeps all rows)
% For the normal copula, Method="phase" calibrates against coherent
% phase-randomized surrogates and keeps every row, so it controls size without
% the power loss of decimation.
phiPhase = 0.9;
rejPhase = 0;
rejDecimate = 0;
nThinned = n;
for t = 1:T
    Y = arGaussianNull(n, d, phiPhase, R, 60000 + t);
    phaseResult = gofcopula.runTestSerial("gofCvM", "normal", Y, ...
        Method="phase", M=M, Seed=82000 + t);
    [decimateResult, decimateInfo] = gofcopula.runTestSerial("gofCvM", "normal", Y, ...
        M=M, Seed=82000 + t);
    rejPhase = rejPhase + double(phaseResult.Tests.PValue < alpha);
    rejDecimate = rejDecimate + double(decimateResult.Tests.PValue < alpha);
    nThinned = decimateInfo.nThinned;
end
fprintf("\nPhase vs decimate at phi=%.1f (n=%d, T=%d, M=%d)\n", phiPhase, n, T, M);
fprintf("  phase    (all %4d rows) size: %.3f\n", n, rejPhase / T);
fprintf("  decimate (~%4d rows)    size: %.3f\n", nThinned, rejDecimate / T);
fprintf("Both should be ~ %.2f; phase keeps every row (more power to detect real departures).\n", ...
    alpha);

%% The correction is inert on i.i.d. data (size preserved for every family)
% On independent draws the estimated interval is about 1, so runTestSerial
% and runTest agree and both keep nominal size -- for any copula family.
families = ["normal" "clayton" "frank" "gumbel"];
thetas = [0.5 2.0 4.0 1.8];
Ti = 60;
ni = 800;

fprintf("\nI.i.d. copula draws, no serial dependence (n=%d, T=%d, M=%d)\n", ni, Ti, M);
fprintf("%10s %12s %12s %14s\n", "family", "iid size", "serial size", "mean interval");
for k = 1:numel(families)
    family = families(k);
    theta = thetas(k);
    rejIID = 0;
    rejSerial = 0;
    intervals = zeros(Ti, 1);
    for t = 1:Ti
        U = gofcopula.copulaRandom(family, ni, theta, ...
            Stream=RandStream("Threefry", "Seed", 90000 + t));
        iidResult = gofcopula.runTest("gofCvM", family, U, ...
            M=M, Margins="none", Seed=91000 + t);
        [serialResult, info] = gofcopula.runTestSerial("gofCvM", family, U, ...
            M=M, Margins="none", Seed=91000 + t);
        rejIID = rejIID + double(iidResult.Tests.PValue < alpha);
        rejSerial = rejSerial + double(serialResult.Tests.PValue < alpha);
        intervals(t) = info.thinInterval;
    end
    fprintf("%10s %12.3f %12.3f %14.2f\n", ...
        family, rejIID / Ti, rejSerial / Ti, mean(intervals));
end

function Y = arGaussianNull(n, d, phi, R, seed)
%ARGAUSSIANNULL Gaussian(R) copula with AR(1) serial dependence per column.
%   Each column is a unit-variance AR(1) process with ACF phi^|h|; the
%   spatial mix by chol(R) makes every row ~ N(0,R), i.e. an exact
%   Gaussian(R) cross-sectional copula. The whole vector is a VAR(1) process,
%   so the copula is the null while the rows are serially dependent.
stream = RandStream("Threefry", "Seed", seed);
burn = 200;
e = randn(stream, n + burn, d);
z = filter(sqrt(1 - phi^2), [1 -phi], e);
z = z(burn + 1:end, :);
Y = z * chol(R);
end
