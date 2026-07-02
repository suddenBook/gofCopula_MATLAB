function statistic = kendallStatistics(method, U, referenceValues)
%KENDALLSTATISTICS CvM/KS statistics for the empirical Kendall process.
% REFERENCEVALUES contains empirical-copula values from a large model sample,
% as in gofCopula 0.4-3. It is not the sample's raw copula observations.

n = size(U,1);
Cn = gofcopula.internal.statistics.empiricalCopula(U, U);
grid = (1:n).' ./ n;
Kn = gofcopula.internal.statistics.empiricalCopula(Cn, grid);
Kmodel = gofcopula.internal.statistics.empiricalCopula(referenceValues(:), grid);

switch lower(string(method))
    case "cvm"
        increments = Kmodel(2:n) - Kmodel(1:n-1);
        squaredIncrements = Kmodel(2:n).^2 - Kmodel(1:n-1).^2;
        statistic = n/3 + n*sum(Kn(1:n-1).^2 .* increments) - ...
            n*sum(Kn(1:n-1) .* squaredIncrements);
    case "ks"
        left = abs(Kn(1:n-1) - Kmodel(1:n-1));
        right = abs(Kn(1:n-1) - Kmodel(2:n));
        statistic = sqrt(n) * max([left; right]);
    otherwise
        error("gofcopula:statistics:InternalMethod", ...
            "Unknown Kendall statistic '%s'.", method);
end
end
