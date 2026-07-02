function C = copulaCDF(family, U, theta, options)
%COPULACDF Evaluate a copula cumulative distribution function.
%   C = gofcopula.copulaCDF(FAMILY,U,THETA) evaluates the named family at
%   rows of U. Elliptical THETA may be a correlation matrix or parameters
%   of the selected Dispersion structure.
%
%   Name-value options:
%     DF          - Student/t-EV degrees of freedom (default 4)
%     Dispersion  - "unstructured", "exchangeable", "ar1", or "toeplitz"
%     Rotation    - bivariate rotation: 0, 90, 180, or 270

arguments
    family {mustBeTextScalar}
    U {mustBeFloat, mustBeReal, mustBeFinite}
    theta {mustBeFloat, mustBeReal, mustBeFinite}
    options.DF (1,1) {mustBeFloat, mustBeReal, mustBeFinite, mustBePositive} = 4
    options.Dispersion (1,1) string = "unstructured"
    options.Rotation (1,1) {mustBeNumeric, mustBeReal, mustBeMember(options.Rotation,[0 90 180 270])} = 0
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

if options.Rotation == 0
    C = gofcopula.internal.copulas.baseCDF( ...
        family, U, theta, options.DF, options.Dispersion);
else
    C = gofcopula.internal.copulas.rotatedCDF( ...
        family, U, theta, options.DF, options.Dispersion, options.Rotation);
end
end
