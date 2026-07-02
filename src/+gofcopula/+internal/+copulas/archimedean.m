function value = archimedean(action, family, x, theta, order)
%ARCHIMEDEAN Stable primitives for supported Archimedean copulas.
%   PHI maps uniforms to the generator scale. LOGPHIPRIME returns the log
%   absolute derivative. PSI evaluates the inverse generator. LOGPSIDERIV
%   returns log(abs(psi^(order)(x))).

arguments
    action (1,1) string
    family (1,1) string
    x {mustBeFloat, mustBeReal}
    theta (1,1) {mustBeFloat, mustBeReal, mustBeFinite}
    order (1,1) {mustBeInteger, mustBeNonnegative, mustBeReal} = 0
end

switch action
    case "phi"
        u = x;
        switch family
            case "clayton"
                if theta == 0, value = -log(u); else, value = expm1(-theta .* log(u)) ./ theta; end
            case "gumbel"
                value = (-log(u)) .^ theta;
            case "frank"
                if theta == 0, value = -log(u); else, value = -log(expm1(-theta .* u) ./ expm1(-theta)); end
            case "joe"
                value = -log1p(-(1 - u) .^ theta);
            case "amh"
                value = log1p(-theta .* (1 - u)) - log(u);
        end
    case "logphiprime"
        u = x;
        switch family
            case "clayton"
                value = (-theta - 1) .* log(u);
            case "gumbel"
                value = log(theta) + (theta - 1) .* log(-log(u)) - log(u);
            case "frank"
                if theta == 0, value = -log(u); else, value = log(abs(theta)) - theta .* u - log(abs(expm1(-theta .* u))); end
            case "joe"
                value = log(theta) + (theta - 1) .* log1p(-u) ...
                    - log1p(-(1 - u) .^ theta);
            case "amh"
                value = log1p(-theta) - log(u) - log1p(-theta .* (1 - u));
        end
    case "psi"
        t = x;
        switch family
            case "clayton"
                if theta == 0
                    value = exp(-t);
                else
                    base = max(1 + theta .* t, 0);
                    value = base .^ (-1/theta);
                end
            case "gumbel"
                value = exp(-(t .^ (1 / theta)));
            case "frank"
                if theta == 0, value = exp(-t); else, value = -log1p(expm1(-theta) .* exp(-t)) ./ theta; end
            case "joe"
                value = 1 - (1 - exp(-t)) .^ (1 / theta);
            case "amh"
                value = (1 - theta) ./ (exp(t) - theta);
        end
    case "logpsideriv"
        value = logPsiDerivative(family, x, theta, order);
    otherwise
        error("gofcopula:copula:InternalAction", "Unknown Archimedean action.");
end
end

function out = logPsiDerivative(family, t, theta, order)
if order == 0
    out = log(archimedean("psi", family, t, theta));
    return
end
switch family
    case "clayton"
        if theta == 0
            out = -t;
        else
            factors = 1 + (0:order-1) .* theta;
            base = 1 + theta .* t;
            out = -inf(size(t));
            inside = base > 0;
            out(inside) = sum(log(abs(factors))) ...
                - (order + 1/theta) .* log(base(inside));
        end
    case "gumbel"
        alpha = 1 / theta;
        polynomial = 1;
        for n = 0:order-1
            derivative = (1:numel(polynomial)-1) .* polynomial(2:end);
            term1 = -n .* polynomial;
            term2 = -alpha .* [0 polynomial];
            term3 = alpha .* [0 derivative];
            lengthOut = max([numel(term1), numel(term2), numel(term3)]);
            polynomial = pad(term1, lengthOut) + pad(term2, lengthOut) + pad(term3, lengthOut);
        end
        z = t .^ alpha;
        p = polyval(fliplr(polynomial), z);
        out = -z - order .* log(t) + log(abs(p));
    case "frank"
        if theta == 0
            out = -t;
        else
            q = -expm1(-theta) .* exp(-t);
            li = negativePolylog(q, order - 1);
            out = log(abs(li ./ theta));
        end
    case "amh"
        if theta == 0
            out = -t;
        else
            q = theta .* exp(-t);
            li = negativePolylog(q, order);
            out = log(abs((1 - theta) .* li ./ theta));
        end
    case "joe"
        if theta == 1
            out = -t;
        else
            out = joeDerivativeLog(t, 1/theta, order);
        end
end
end

function y = pad(x, n)
y = zeros(1, n, "like", x);
y(1:numel(x)) = x;
end

function value = negativePolylog(q, m)
% Li_{-m}(q), using its Eulerian-polynomial rational representation.
if m == 0
    value = q ./ (1 - q);
    return
end
eulerian = 1;
for n = 2:m
    next = zeros(1, n);
    for k = 0:n-1
        left = 0; right = 0;
        if k + 1 <= numel(eulerian), left = (k + 1) * eulerian(k + 1); end
        if k >= 1, right = (n - k) * eulerian(k); end
        next(k + 1) = left + right;
    end
    eulerian = next;
end
value = q .* polyval(fliplr(eulerian), q) ./ (1 - q) .^ (m + 1);
end

function out = joeDerivativeLog(t, alpha, order)
% Positive binomial series for (-1)^order psi^(order), 0 < alpha < 1.
out = zeros(size(t));
for i = 1:numel(t)
    q = exp(-t(i));
    coefficient = alpha;
    term = coefficient * q;
    total = term;
    k = 1;
    while k < 200000
        ratio = ((k - alpha) / (k + 1)) * q * ((k + 1) / k) ^ order;
        term = term * ratio;
        totalNew = total + term;
        k = k + 1;
        if term <= 8 * eps(totalNew) * totalNew, break, end
        total = totalNew;
    end
    out(i) = log(totalNew);
end
end
