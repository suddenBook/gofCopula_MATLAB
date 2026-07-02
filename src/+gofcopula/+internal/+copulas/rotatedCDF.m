function C = rotatedCDF(family, U, theta, df, dispersion, rotation)
%ROTATEDCDF Evaluate a possibly rotated bivariate copula CDF.

base = @(X) gofcopula.internal.copulas.baseCDF(family, X, theta, df, dispersion);
u = U(:,1); v = U(:,2);
switch rotation
    case 0
        C = base(U);
    case 90
        C = v - base([1-u, v]);
    case 180
        C = u + v - 1 + base([1-u, 1-v]);
    case 270
        C = u - base([u, 1-v]);
end
C = min(max(C, 0), min(U, [], 2));
end
