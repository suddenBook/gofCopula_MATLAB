function fitted = modelWithTheta(model, theta)
%MODELWITHTHETA Copy a model while replacing its copula parameter vector.

if isa(model,"gofcopula.CopulaModel")
    if model.Family == "powerexp"
        % Trailing element is the shape beta (see numericalInformation).
        fitted = model.withFit(theta(1:end-1),theta(end));
    else
        fitted = model.withFit(theta,model.DegreesOfFreedom);
    end
elseif isstruct(model)
    fitted = model;
    fitted.Theta = theta;
else
    error("gofcopula:statistics:InvalidModel", "Invalid copula model.");
end
end
