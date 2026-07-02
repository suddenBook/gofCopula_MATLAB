function U = modelRandom(model, n, stream, dimension)
%MODELRANDOM Draw from the public copula random-number interface.

arguments
    model
    n (1,1) double {mustBeInteger,mustBePositive}
    stream = []
    dimension (1,1) double {mustBeInteger,mustBeGreaterThan(dimension,1)} = 2
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
    error("gofcopula:statistics:InvalidModel", "Invalid copula model.");
end
U = gofcopula.copulaRandom(family, n, theta, Dimension=dimension, DF=df, ...
    Dispersion=dispersion, Rotation=rotation, Stream=stream);
end

function value=fieldOr(s,name,default)
if isfield(s,name), value=s.(name); else, value=default; end
end
