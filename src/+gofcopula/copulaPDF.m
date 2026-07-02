function density = copulaPDF(family, U, theta, options)
%COPULAPDF Evaluate a copula density or log-density.
%   F = gofcopula.copulaPDF(FAMILY,U,THETA) returns the density at rows of U.
%   Use Log=true to return the numerically stable log-density.

arguments
    family {mustBeTextScalar}
    U {mustBeFloat, mustBeReal, mustBeFinite}
    theta {mustBeFloat, mustBeReal, mustBeFinite}
    options.DF (1,1) {mustBeFloat, mustBeReal, mustBeFinite, mustBePositive} = 4
    options.Dispersion (1,1) string = "unstructured"
    options.Rotation (1,1) {mustBeNumeric, mustBeReal, mustBeMember(options.Rotation,[0 90 180 270])} = 0
    options.Log (1,1) logical = false
end
if ~ismatrix(U) || isempty(U)
    error("gofcopula:copula:InvalidData", "U must be a nonempty numeric matrix.");
end
if any(U(:) < 0 | U(:) > 1)
    error("gofcopula:copula:OutsideUnitCube", "Every entry of U must lie in [0,1].");
end
family = gofcopula.internal.copulas.normalizeFamily(family);
dimension = size(U,2);
gofcopula.internal.copulas.validateModel( ...
    family, dimension, theta, options.DF, options.Dispersion, options.Rotation);

X = U;
switch options.Rotation
    case 90
        X(:,1) = 1-X(:,1);
    case 180
        X = 1-X;
    case 270
        X(:,2) = 1-X(:,2);
end
logDensity = gofcopula.internal.copulas.baseLogPDF( ...
    family, X, theta, options.DF, options.Dispersion);
if options.Log, density = logDensity; else, density = exp(logDensity); end
end
