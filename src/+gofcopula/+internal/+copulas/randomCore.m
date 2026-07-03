function U = randomCore(family, n, dimension, theta, df, dispersion, rotation, stream)
%RANDOMCORE Draw copula observations.
%   Sampling strategy per family (fast paths validated against the
%   bisection reference randomCoreBisection by differential tests):
%     normal/t          - correlated Gaussian / Student scores (vectorized)
%     clayton, frank    - d=2: closed-form conditional inversion (any
%                         admissible theta); d>=3: Marshall--Olkin frailty
%     gumbel, joe, amh  - Marshall--Olkin frailty (amh theta<0 bisects)
%     fgm               - closed-form conditional inversion
%     EV families,
%     plackett, others  - conditional inversion by vectorized bisection
%   Near-independence parameters short-circuit to uniform draws. All
%   randomness comes from STREAM ([] uses the global stream).

if family == "powerexp"
    % Elliptical power-exponential copula; df carries the shape beta.
    R = gofcopula.internal.copulas.correlationMatrix(theta, dimension, dispersion);
    U = gofcopula.internal.elliptical.peCopulaRandom(n, R, df, stream);
elseif ismember(family, ["normal", "t"])
    R = gofcopula.internal.copulas.correlationMatrix(theta, dimension, dispersion);
    L = chol(R, "lower");
    z = normalRandom(stream, n, dimension) * L.';
    if family == "normal"
        U = normcdf(z);
    else
        chiSquare = 2 .* gammaincinv(uniformRandom(stream,n,1), df/2, "lower");
        U = tcdf(z ./ sqrt(chiSquare ./ df), df);
    end
elseif isIndependence(family, theta)
    U = uniformRandom(stream, n, dimension);
elseif family == "clayton" && dimension == 2
    U = conditionalClayton(stream, n, theta);
elseif family == "frank" && dimension == 2
    U = conditionalFrank(stream, n, theta);
elseif family == "fgm"
    U = conditionalFGM(stream, n, theta);
elseif ismember(family, ["clayton", "gumbel", "frank", "joe"]) || ...
        (family == "amh" && theta > 0)
    U = marshallOlkin(family, stream, n, dimension, theta);
else
    U = bisectionSample(family, stream, n, dimension, theta, df, dispersion);
end

if rotation ~= 0
    switch rotation
        case 90
            U(:,1) = 1-U(:,1);
        case 180
            U = 1-U;
        case 270
            U(:,2) = 1-U(:,2);
    end
end
U = min(max(U, 0), 1);
end

function tf = isIndependence(family, theta)
switch family
    case {"clayton", "frank", "amh", "fgm", "galambos", "huslerreiss", "tawn"}
        tf = abs(theta) < sqrt(eps);
    case {"gumbel", "joe"}
        tf = abs(theta - 1) < sqrt(eps);
    case "plackett"
        tf = abs(theta - 1) < sqrt(eps);
    otherwise
        tf = false;
end
end

function U = marshallOlkin(family, stream, n, dimension, theta)
% U_j = psi(E_j / V) with the standard-scale generator inverse psi.
V = gofcopula.internal.copulas.frailty(family, n, theta, stream);
E = -log(uniformRandom(stream, n, dimension));
T = E ./ V;
switch family
    case "clayton"
        U = (1 + T).^(-1/theta);
    case "gumbel"
        U = exp(-T.^(1/theta));
    case "frank"
        U = -log1p(exp(-T) .* expm1(-theta)) ./ theta;
    case "joe"
        U = 1 - (-expm1(-T)).^(1/theta);
    case "amh"
        U = (1 - theta) ./ (exp(T) - theta);
end
end

function U = conditionalClayton(stream, n, theta)
% v = [1 + u^-theta*(w^(-theta/(1+theta)) - 1)]^(-1/theta), exact for
% theta in (-1,0) and theta > 0.
if theta < -0.999
    U = bisectionSample("clayton", stream, n, 2, theta, 4, "unstructured");
    return
end
u = uniformRandom(stream, n, 1);
w = uniformRandom(stream, n, 1);
base = 1 + u.^(-theta) .* (w.^(-theta/(1+theta)) - 1);
v = max(base, realmin).^(-1/theta);
U = [u, min(max(v, 0), 1)];
end

function U = conditionalFrank(stream, n, theta)
% v = -(1/theta)*log1p(w*(exp(-theta)-1) / ((1-w)*exp(-theta*u) + w)).
u = uniformRandom(stream, n, 1);
w = uniformRandom(stream, n, 1);
v = -log1p(w .* expm1(-theta) ./ ((1 - w) .* exp(-theta .* u) + w)) ./ theta;
U = [u, min(max(v, 0), 1)];
end

function U = conditionalFGM(stream, n, theta)
% h(v|u) = A*v^2 + (1-A)*v with A = theta*(2u-1); invert the quadratic.
u = uniformRandom(stream, n, 1);
w = uniformRandom(stream, n, 1);
A = theta .* (2 .* u - 1);
v = w;
solvable = abs(A) >= 1e-9;
As = A(solvable);
v(solvable) = ((As - 1) + sqrt((1 - As).^2 + 4 .* As .* w(solvable))) ./ (2 .* As);
U = [u, min(max(v, 0), 1)];
end

function U = bisectionSample(family, stream, n, dimension, theta, df, dispersion)
% Conditional inversion, bisecting all rows simultaneously: 54 vectorized
% conditionalCDF calls per dimension instead of 54*n scalar calls.
targets = uniformRandom(stream, n, dimension);
U = zeros(n, dimension);
U(:,1) = targets(:,1);
for j = 2:dimension
    low = zeros(n, 1);
    high = ones(n, 1);
    for iteration = 1:54
        middle = (low + high) / 2;
        value = gofcopula.internal.copulas.conditionalCDF( ...
            family, [U(:,1:j-1), middle], theta, df, dispersion);
        below = value < targets(:,j);
        low(below) = middle(below);
        high(~below) = middle(~below);
    end
    U(:,j) = (low + high) / 2;
end
end

function u = uniformRandom(stream, varargin)
if isempty(stream), u = rand(varargin{:}); else, u = rand(stream, varargin{:}); end
end

function z = normalRandom(stream, varargin)
if isempty(stream), z = randn(varargin{:}); else, z = randn(stream, varargin{:}); end
end
