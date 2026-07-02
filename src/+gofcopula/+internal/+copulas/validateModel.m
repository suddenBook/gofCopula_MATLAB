function validateModel(family, dimension, theta, df, dispersion, rotation)
%VALIDATEMODEL Validate family-specific dimensions and parameter domains.

if dimension < 2 || dimension ~= fix(dimension)
    error("gofcopula:copula:InvalidDimension", ...
        "Copula dimension must be an integer greater than or equal to two.");
end
if ~ismember(rotation, [0 90 180 270])
    error("gofcopula:copula:InvalidRotation", ...
        "Rotation must be 0, 90, 180, or 270 degrees.");
end
if rotation ~= 0 && dimension ~= 2
    error("gofcopula:copula:RotationDimension", ...
        "Rotations are defined only for bivariate copulas.");
end
if ~isfinite(df) || df <= 0
    error("gofcopula:copula:InvalidDF", ...
        "Degrees of freedom must be positive and finite.");
end

bivariateOnly = ["galambos", "huslerreiss", "tawn", "tev", "fgm", "plackett"];
if dimension > 2 && ismember(family, bivariateOnly)
    error("gofcopula:copula:UnsupportedDimension", ...
        "The %s copula is implemented only in dimension two.", family);
end

if ismember(family, ["normal", "t"])
    gofcopula.internal.copulas.correlationMatrix(theta, dimension, dispersion);
    return
end
if ~isscalar(theta)
    error("gofcopula:copula:InvalidParameter", ...
        "The %s copula requires one scalar dependence parameter.", family);
end

switch family
    case "clayton"
        valid = theta >= (dimension == 2) * -1;
        description = "theta >= -1 in dimension two and theta >= 0 otherwise";
    case "gumbel"
        valid = theta >= 1;
        description = "theta >= 1";
    case "frank"
        valid = dimension == 2 || theta >= 0;
        description = "finite theta (and theta >= 0 above dimension two)";
    case "joe"
        valid = theta >= 1;
        description = "theta >= 1";
    case "amh"
        valid = theta >= -1 && theta < 1 && (dimension == 2 || theta >= 0);
        description = "-1 <= theta < 1 (and theta >= 0 above dimension two)";
    case {"galambos", "huslerreiss"}
        valid = theta >= 0;
        description = "theta >= 0";
    case "tawn"
        valid = theta >= 0 && theta <= 1;
        description = "0 <= theta <= 1";
    case "tev"
        valid = theta > -1 && theta < 1;
        description = "-1 < theta < 1";
    case "fgm"
        valid = abs(theta) <= 1;
        description = "-1 <= theta <= 1";
    case "plackett"
        valid = theta > 0;
        description = "theta > 0";
    otherwise
        valid = false;
        description = "a valid parameter";
end
if ~valid
    error("gofcopula:copula:InvalidParameter", ...
        "The %s copula requires %s.", family, description);
end
end
