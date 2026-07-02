function C = baseCDF(family, U, theta, df, dispersion)
%BASECDF Evaluate an unrotated copula distribution function.

[n, dimension] = size(U);
C = zeros(n, 1, "like", U);
zeroRows = any(U == 0, 2);
active = ~zeroRows;
if dimension == 2
    firstMargin = active & U(:,2) == 1;
    secondMargin = active & U(:,1) == 1;
    C(firstMargin) = U(firstMargin,1);
    C(secondMargin) = U(secondMargin,2);
    active = active & ~firstMargin & ~secondMargin;
end
if ~any(active), return, end
X = U(active, :);

switch family
    case "normal"
        R = gofcopula.internal.copulas.correlationMatrix(theta, dimension, dispersion);
        C(active) = mvncdf(norminv(X), zeros(1, dimension), R);
    case "t"
        R = gofcopula.internal.copulas.correlationMatrix(theta, dimension, dispersion);
        C(active) = mvtcdf(tinv(X, df), R, df);
    case {"clayton", "gumbel", "frank", "joe", "amh"}
        phi = gofcopula.internal.copulas.archimedean("phi", family, X, theta);
        C(active) = gofcopula.internal.copulas.archimedean( ...
            "psi", family, sum(phi, 2), theta);
    case "fgm"
        u = X(:, 1); v = X(:, 2);
        C(active) = u .* v .* (1 + theta .* (1-u) .* (1-v));
    case "plackett"
        u = X(:, 1); v = X(:, 2);
        if theta == 1
            C(active) = u .* v;
        else
            eta = theta - 1;
            b = 1 + eta .* (u + v);
            discriminant = max(0, b.^2 - 4 .* theta .* eta .* u .* v);
            % Rationalization avoids cancellation when theta is near one.
            C(active) = (2 .* theta .* u .* v) ./ (b + sqrt(discriminant));
        end
    case {"galambos", "huslerreiss", "tawn", "tev"}
        u = X(:, 1); v = X(:, 2);
        loguv = log(u) + log(v);
        w = log(v) ./ loguv;
        A = gofcopula.internal.copulas.evDependence(family, w, theta, df);
        C(active) = exp(loguv .* A);
end
C = min(max(C, 0), min(U, [], 2));
end
