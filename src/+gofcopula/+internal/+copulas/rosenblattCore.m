function Z = rosenblattCore(family, U, theta, df, dispersion, rotation)
%ROSENBLATTCORE Compute a sequential Rosenblatt transform.

[n, dimension] = size(U);
Z = zeros(n, dimension, "like", U);
Z(:,1) = U(:,1);

if rotation ~= 0
    baseU = U;
    switch rotation
        case 90
            baseU(:,1) = 1-U(:,1);
            value = stepValue(family, baseU, theta, df, dispersion, dimension);
        case 180
            baseU = 1-U;
            value = 1-stepValue(family, baseU, theta, df, dispersion, dimension);
        case 270
            baseU(:,2) = 1-U(:,2);
            value = 1-stepValue(family, baseU, theta, df, dispersion, dimension);
    end
    Z(:,2) = min(max(value, 0), 1);
    return
end

for j = 2:dimension
    Z(:,j) = stepValue(family, U(:,1:j), theta, df, dispersion, dimension);
end
end

function value = stepValue(family, prefix, theta, df, dispersion, dimension)
% Conditional CDF P(U_j <= prefix(:,j) | prefix(:,1:j-1)) for one Rosenblatt step.
% Power-exponential conditionals depend on the ambient dimension (its margins are
% not power-exponential), so they are handled here where DIMENSION is known.
if family == "powerexp"
    R = gofcopula.internal.copulas.correlationMatrix(theta, dimension, dispersion);
    value = gofcopula.internal.elliptical.peConditionalCDF(prefix, R, df, dimension);
    return
end
conditionalTheta = theta;
if ismember(family, ["normal", "t"])
    fullR = gofcopula.internal.copulas.correlationMatrix(theta, dimension, dispersion);
    conditionalTheta = fullR(1:size(prefix,2), 1:size(prefix,2));
end
value = gofcopula.internal.copulas.conditionalCDF( ...
    family, prefix, conditionalTheta, df, dispersion);
end
