function logDensity = baseLogPDF(family, U, theta, df, dispersion)
%BASELOGPDF Evaluate the log-density of an unrotated copula.

[~, dimension] = size(U);
% Formulas are evaluated by their interior limits at exact cube boundaries.
X = min(max(U, realmin(class(U))), 1 - eps(class(U)));

switch family
    case "normal"
        R = gofcopula.internal.copulas.correlationMatrix(theta, dimension, dispersion);
        z = norminv(X);
        L = chol(R, "lower");
        solved = L \ z.';
        quadratic = sum(solved.^2, 1).' - sum(z.^2, 2);
        logDensity = -sum(log(diag(L))) - 0.5 .* quadratic;
    case "t"
        R = gofcopula.internal.copulas.correlationMatrix(theta, dimension, dispersion);
        z = tinv(X, df);
        L = chol(R, "lower");
        solved = L \ z.';
        q = sum(solved.^2, 1).';
        logJoint = gammaln((df + dimension)/2) - gammaln(df/2) ...
            - sum(log(diag(L))) - dimension/2 * log(df*pi) ...
            - (df + dimension)/2 .* log1p(q ./ df);
        logMarginal = gammaln((df + 1)/2) - gammaln(df/2) ...
            - 0.5*log(df*pi) - (df + 1)/2 .* log1p(z.^2 ./ df);
        logDensity = logJoint - sum(logMarginal, 2);
    case "powerexp"
        % Elliptical power-exponential copula; df carries the shape beta.
        R = gofcopula.internal.copulas.correlationMatrix(theta, dimension, dispersion);
        logDensity = gofcopula.internal.elliptical.peCopulaLogPDF(X, R, df);
    case {"clayton", "gumbel", "frank", "joe", "amh"}
        phi = gofcopula.internal.copulas.archimedean("phi", family, X, theta);
        logPhiPrime = gofcopula.internal.copulas.archimedean( ...
            "logphiprime", family, X, theta);
        logPsiDerivative = gofcopula.internal.copulas.archimedean( ...
            "logpsideriv", family, sum(phi, 2), theta, dimension);
        logDensity = logPsiDerivative + sum(logPhiPrime, 2);
    case "fgm"
        density = 1 + theta .* (1 - 2*X(:,1)) .* (1 - 2*X(:,2));
        logDensity = log(density);
    case "plackett"
        u = X(:,1); v = X(:,2);
        if theta == 1
            logDensity = zeros(size(u));
        else
            eta = theta - 1;
            delta = (1 + eta .* (u+v)).^2 - 4 .* theta .* eta .* u .* v;
            logDensity = log(theta) + log1p(eta .* (u+v-2*u.*v)) ...
                - 1.5 .* log(delta);
        end
    case {"galambos", "huslerreiss", "tawn", "tev"}
        u = X(:,1); v = X(:,2);
        s = log(u) + log(v);
        w = log(v) ./ s;
        [A, A1, A2] = gofcopula.internal.copulas.evDependence(family, w, theta, df);
        factor = (A - w.*A1) .* (A + (1-w).*A1) ...
            - w.*(1-w).*A2 ./ s;
        logDensity = s .* A - log(u) - log(v) + log(max(factor, realmin));
end
end
