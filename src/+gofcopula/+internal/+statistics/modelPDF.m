function f = modelPDF(model, U, logFlag)
%MODELPDF Evaluate the public copula density from a CopulaModel specification.

arguments
    model
    U (:,:) double
    logFlag (1,1) logical = false
end

if isa(model, "gofcopula.CopulaModel")
    family = model.Family; theta = model.Theta;
    df = model.DegreesOfFreedom; dispersion = model.Dispersion;
    rotation = model.Rotation;
elseif isstruct(model)
    family = string(model.Family); theta = model.Theta;
    df = fieldOr(model,"DegreesOfFreedom",4);
    dispersion = string(fieldOr(model,"Dispersion","exchangeable"));
    rotation = fieldOr(model,"Rotation",0);
else
    error("gofcopula:statistics:InvalidModel", ...
        "Model must be a gofcopula.CopulaModel or a compatible structure.");
end
f = gofcopula.copulaPDF(family, U, theta, DF=df, ...
    Dispersion=dispersion, Rotation=rotation, Log=logFlag);
end

function value = fieldOr(s,name,default)
if isfield(s,name), value=s.(name); else, value=default; end
end
