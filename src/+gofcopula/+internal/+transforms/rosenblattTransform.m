function Z = rosenblattTransform(U, model)
%ROSENBLATTTRANSFORM Conditional probability transform for a copula sample.
% Under a correctly specified continuous copula, rows of Z are independent
% observations from the product uniform distribution.

arguments
    U (:,:) double {mustBeReal,mustBeFinite}
    model
end

if isa(model,"gofcopula.CopulaModel")
    family=model.Family; theta=model.Theta; df=model.DegreesOfFreedom;
    dispersion=model.Dispersion; rotation=model.Rotation;
elseif isstruct(model)
    family=string(model.Family); theta=model.Theta;
    df=fieldOr(model,"DegreesOfFreedom",4);
    dispersion=string(fieldOr(model,"Dispersion","exchangeable"));
    rotation=fieldOr(model,"Rotation",0);
else
    error("gofcopula:transforms:InvalidModel", "Invalid copula model.");
end

Z = gofcopula.rosenblatt(family, U, theta, DF=df, ...
    Dispersion=dispersion, Rotation=rotation);
Z = min(max(Z,0),1); % remove harmless roundoff outside the unit cube
end

function value=fieldOr(s,name,default)
if isfield(s,name), value=s.(name); else, value=default; end
end
