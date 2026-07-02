function U = copulaRandom(family, n, theta, options)
%COPULARANDOM Generate observations from a supported copula.
%   U = gofcopula.copulaRandom(FAMILY,N,THETA,Dimension=D) returns N-by-D
%   observations. Stream may be a RandStream for isolated reproducibility.

arguments
    family {mustBeTextScalar}
    n (1,1) {mustBeNumeric, mustBeReal, mustBeFinite, mustBeInteger, mustBeNonnegative}
    theta {mustBeFloat, mustBeReal, mustBeFinite}
    options.Dimension (1,1) {mustBeNumeric, mustBeReal, mustBeFinite, mustBeInteger, mustBePositive} = 2
    options.DF (1,1) {mustBeFloat, mustBeReal, mustBeFinite, mustBePositive} = 4
    options.Dispersion (1,1) string = "unstructured"
    options.Rotation (1,1) {mustBeNumeric, mustBeReal, mustBeMember(options.Rotation,[0 90 180 270])} = 0
    options.Stream = []
end
if ~isempty(options.Stream) && ~isa(options.Stream, "RandStream")
    error("gofcopula:copula:InvalidStream", "Stream must be a RandStream object or empty.");
end
family = gofcopula.internal.copulas.normalizeFamily(family);
gofcopula.internal.copulas.validateModel( ...
    family, options.Dimension, theta, options.DF, options.Dispersion, options.Rotation);
U = gofcopula.internal.copulas.randomCore( ...
    family, n, options.Dimension, theta, options.DF, options.Dispersion, ...
    options.Rotation, options.Stream);
end
