function fitted = modelWithTheta(model, theta)
%MODELWITHTHETA Copy a model while replacing its copula parameter vector.

if isa(model,"gofcopula.CopulaModel")
    fitted = model.withFit(theta,model.DegreesOfFreedom);
elseif isstruct(model)
    fitted = model;
    fitted.Theta = theta;
else
    error("gofcopula:statistics:InvalidModel", "Invalid copula model.");
end
end
