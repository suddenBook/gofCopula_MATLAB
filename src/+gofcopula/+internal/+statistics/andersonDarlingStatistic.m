function statistic = andersonDarlingStatistic(method, U)
%ANDERSONDARLINGSTATISTIC AD statistic after chi-square or gamma mapping.

[n,d] = size(U);
switch lower(string(method))
    case "chisq"
        z = -sqrt(2) .* erfcinv(2 .* U);
        mapped = gammainc(sum(z.^2,2) ./ 2, d/2, "lower");
    case "gamma"
        mapped = gammainc(sum(-log(U),2), d, "lower");
    otherwise
        error("gofcopula:statistics:InternalMethod", ...
            "Unknown Anderson--Darling mapping '%s'.", method);
end

% Inputs generated from pseudo-observations are interior. Clipping protects
% against log(0) if callers deliberately supply an exact boundary value.
mapped = sort(min(max(mapped, realmin("double")), 1-eps("double")));
j = (1:n).';
statistic = -n - sum(((2*j-1)./n) .* ...
    (log(mapped) + log1p(-mapped(n+1-j))));
end
