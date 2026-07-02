function statistic = empiricalStatistics(method, U, model)
%EMPIRICALSTATISTICS Empirical-copula KS and Cramer--von Mises statistics.

n = size(U,1);
Cn = gofcopula.internal.statistics.empiricalCopula(U, U);
Ctheta = gofcopula.internal.statistics.modelCDF(model, U);
residual = Cn - Ctheta;

switch lower(string(method))
    case "ks"
        statistic = sqrt(n) * max(abs(residual));
    case "cvm"
        % n * integral (Cn-Ctheta)^2 dCn is the unscaled sample sum.
        statistic = sum(residual.^2);
    otherwise
        error("gofcopula:statistics:InternalMethod", ...
            "Unknown empirical-copula statistic '%s'.", method);
end
end
