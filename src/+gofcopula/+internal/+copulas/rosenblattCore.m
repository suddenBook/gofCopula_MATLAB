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
            value = gofcopula.internal.copulas.conditionalCDF( ...
                family, baseU, theta, df, dispersion);
        case 180
            baseU = 1-U;
            value = 1-gofcopula.internal.copulas.conditionalCDF( ...
                family, baseU, theta, df, dispersion);
        case 270
            baseU(:,2) = 1-U(:,2);
            value = 1-gofcopula.internal.copulas.conditionalCDF( ...
                family, baseU, theta, df, dispersion);
    end
    Z(:,2) = min(max(value, 0), 1);
    return
end

for j = 2:dimension
    conditionalTheta = theta;
    if ismember(family, ["normal", "t"])
        fullR = gofcopula.internal.copulas.correlationMatrix(theta, dimension, dispersion);
        conditionalTheta = fullR(1:j,1:j);
    end
    Z(:,j) = gofcopula.internal.copulas.conditionalCDF( ...
        family, U(:,1:j), conditionalTheta, df, dispersion);
end
end
