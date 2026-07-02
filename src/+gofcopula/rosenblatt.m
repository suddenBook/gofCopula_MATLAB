function Z = rosenblatt(family, U, theta, options)
%ROSENBLATT Transform copula observations to independent uniforms.
%   Z(:,1)=U(:,1); subsequent columns are sequential conditional CDFs.

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
Z = gofcopula.internal.copulas.rosenblattCore( ...
    family, U, theta, options.DF, options.Dispersion, options.Rotation);
end
