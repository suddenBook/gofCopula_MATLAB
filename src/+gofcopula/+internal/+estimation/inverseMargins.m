function x = inverseMargins(u, family, parameters)
%INVERSEMARGINS Quantile transform of one column through a fitted margin.
%   X = INVERSEMARGINS(U, FAMILY, PARAMETERS) maps copula-scale values U in
%   (0,1) back to the data scale using the margin parameters produced by
%   transformMargins for the same column. Used by the corrected-mode
%   bootstrap to rebuild data-scale replicates before margins are refitted.

arguments
    u (:,1) double {mustBeReal}
    family (1,1) string
    parameters
end

switch lower(family)
    case {"ranks", "none"}
        x = u;
    case "cauchy"
        x = parameters.Location + parameters.Scale .* tan(pi .* (u - 0.5));
    case "chisq"
        x = chi2inv(u, parameters.DegreesOfFreedom);
    case "f"
        x = finv(u, parameters.DF1, parameters.DF2);
    otherwise
        x = icdf(parameters.Distribution, u);
end
end
