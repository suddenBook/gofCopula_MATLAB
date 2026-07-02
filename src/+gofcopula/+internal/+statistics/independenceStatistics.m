function statistic = independenceStatistics(method, U)
%INDEPENDENCESTATISTICS Genest--Remillard--Beaudoin Sn(B) and Sn(C).
% U is assumed to contain independent-uniform observations, normally after
% a Rosenblatt or Hering--Hofert transform.

[n,d] = size(U);
switch lower(string(method))
    case "snb"
        % Closed form of n int(Dn-prod(u))^2 du. This algebra is identical
        % to copula::gofTstat(method="SnB"). Log-sum-exp avoids losing all
        % contribution when products underflow in moderate/high dimensions.
        sum1 = exp(logSumExp(sum(log1p(-U.^2),2)));
        pairTerm = 0;
        for i = 1:n
            logs = sum(log1p(-max(U,U(i,:))),2);
            pairTerm = pairTerm + exp(logSumExp(logs));
        end
        statistic = n / 3^d - sum1 / 2^(d-1) + pairTerm / n;
    case "snc"
        Dn = gofcopula.internal.statistics.empiricalCopula(U, U);
        statistic = sum((Dn - prod(U,2)).^2);
    otherwise
        error("gofcopula:statistics:InternalMethod", ...
            "Unknown independence statistic '%s'.", method);
end
end

function value = logSumExp(x)
maximum=max(x);
if isinf(maximum) && maximum < 0
    value=-Inf;
else
    value=maximum+log(sum(exp(x-maximum)));
end
end
