function h = conditionalCDF(family, prefix, theta, df, dispersion)
%CONDITIONALCDF Conditional CDF of the last component given its predecessors.
%   PREFIX is n-by-j and the result is P(U_j <= prefix(:,j) | U_1,...,U_{j-1}).

[~, j] = size(prefix);
candidate = prefix(:,j);
h = zeros(size(candidate), "like", prefix);
h(candidate >= 1) = 1;
active = candidate > 0 & candidate < 1;
if ~any(active), return, end
X = min(max(prefix(active,:), realmin(class(prefix))), 1-eps(class(prefix)));

switch family
    case "normal"
        R = gofcopula.internal.copulas.correlationMatrix(theta, j, dispersion);
        z = norminv(X);
        previous = 1:j-1;
        weights = R(j,previous) / R(previous,previous);
        conditionalMean = z(:,previous) * weights.';
        conditionalSD = sqrt(max(realmin, 1 - weights * R(previous,j)));
        value = normcdf((z(:,j) - conditionalMean) ./ conditionalSD);
    case "t"
        R = gofcopula.internal.copulas.correlationMatrix(theta, j, dispersion);
        z = tinv(X, df);
        previous = 1:j-1;
        Rpp = R(previous,previous);
        weights = R(j,previous) / Rpp;
        conditionalMean = z(:,previous) * weights.';
        q = sum((z(:,previous) / chol(Rpp)) .^ 2, 2);
        residualVariance = max(realmin, 1 - weights * R(previous,j));
        conditionalScale = sqrt((df + q) ./ (df + j - 1) .* residualVariance);
        value = tcdf((z(:,j)-conditionalMean)./conditionalScale, df+j-1);
    case {"clayton", "gumbel", "frank", "joe", "amh"}
        phi = gofcopula.internal.copulas.archimedean("phi", family, X, theta);
        numerator = gofcopula.internal.copulas.archimedean( ...
            "logpsideriv", family, sum(phi,2), theta, j-1);
        denominator = gofcopula.internal.copulas.archimedean( ...
            "logpsideriv", family, sum(phi(:,1:j-1),2), theta, j-1);
        value = exp(numerator - denominator);
    otherwise
        if j ~= 2
            error("gofcopula:copula:UnsupportedDimension", ...
                "The %s conditional distribution is bivariate only.", family);
        end
        u = X(:,1); v = X(:,2);
        switch family
            case "fgm"
                value = v + theta .* (2*u-1) .* (v-1) .* v;
            case "plackett"
                eta = theta - 1;
                root = sqrt(max(0, (1 + eta.*(u+v)).^2 ...
                    - 4.*theta.*eta.*u.*v));
                value = 0.5 .* (1 + (-1 - eta.*u + (2+eta).*v) ./ root);
            case {"galambos", "huslerreiss", "tawn", "tev"}
                s = log(u) + log(v);
                w = log(v) ./ s;
                [A, A1] = gofcopula.internal.copulas.evDependence(family,w,theta,df);
                C = exp(s .* A);
                value = C .* (A - w.*A1) ./ u;
        end
end
h(active) = min(max(value, 0), 1);
end
