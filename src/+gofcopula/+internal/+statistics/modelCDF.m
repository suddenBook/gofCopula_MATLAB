function C = modelCDF(model, U)
%MODELCDF Evaluate the public copula CDF from a CopulaModel specification.

[family,theta,df,dispersion,rotation] = modelParts(model);
C = gofcopula.copulaCDF(family, U, theta, DF=df, ...
    Dispersion=dispersion, Rotation=rotation);
end

function [family,theta,df,dispersion,rotation] = modelParts(model)
if isa(model, "gofcopula.CopulaModel")
    family = model.Family;
    theta = model.Theta;
    df = model.DegreesOfFreedom;
    dispersion = model.Dispersion;
    rotation = model.Rotation;
elseif isstruct(model)
    family = string(model.Family);
    theta = model.Theta;
    df = fieldOr(model, "DegreesOfFreedom", 4);
    dispersion = string(fieldOr(model, "Dispersion", "exchangeable"));
    rotation = fieldOr(model, "Rotation", 0);
else
    error("gofcopula:statistics:InvalidModel", ...
        "Model must be a gofcopula.CopulaModel or a compatible structure.");
end
end

function value = fieldOr(s, name, default)
if isfield(s,name)
    value = s.(name);
else
    value = default;
end
end
