function [A, A1, A2] = evDependence(family, w, theta, df)
%EVDEPENDENCE Pickands function and its first two derivatives (analytic).
%   Galambos:      A = 1 - S^(-1/theta), S = w^-theta + (1-w)^-theta
%   Husler-Reiss:  A = w*Phi(z+) + (1-w)*Phi(z-),
%                  z+- = 1/theta +- (theta/2)*log(w/(1-w)); the identity
%                  w*phi(z+) = (1-w)*phi(z-) collapses A' to Phi(z+)-Phi(z-)
%   Tawn:          A = 1 - theta*w*(1-w)
%   t-EV:          A = w*G(x) + (1-w)*G(y), G = t CDF with df+1,
%                  x=(r-theta)*s, y=(1/r-theta)*s, r=(w/(1-w))^(1/df),
%                  s=sqrt((df+1)/(1-theta^2)); derivatives by product rule.
%   All formulas are validated against five-point finite differences in
%   the test suite and by WolframScript derivations in tools/wolfram.

A = evaluateA(family, w, theta, df);
if nargout < 2
    return
end

switch family
    case "tawn"
        A1 = theta .* (2 .* w - 1);
        A2 = 2 .* theta .* ones(size(w));
    case "galambos"
        [A1, A2] = galambosDerivatives(w, theta);
    case "huslerreiss"
        [A1, A2] = huslerReissDerivatives(w, theta);
    case "tev"
        [A1, A2] = tevDerivatives(w, theta, df);
    otherwise
        error("gofcopula:copula:InternalFamily", ...
            "Family '%s' is not an extreme-value family.", family);
end
end

function A = evaluateA(family, w, theta, df)
switch family
    case "galambos"
        if theta == 0
            A = ones(size(w));
        else
            logTerms = [-theta .* log(w(:)), -theta .* log1p(-w(:))];
            maximum = max(logTerms, [], 2);
            logSum = maximum + log(sum(exp(logTerms - maximum), 2));
            A = reshape(1 - exp(-logSum ./ theta), size(w));
        end
    case "huslerreiss"
        if theta == 0
            A = max(w, 1 - w);
        else
            z = 0.5 .* theta .* log(w ./ (1 - w));
            A = w .* normcdf(1/theta + z) ...
                + (1 - w) .* normcdf(1/theta - z);
        end
    case "tawn"
        A = 1 + theta .* w .* (w - 1);
    case "tev"
        ratio = (w ./ (1 - w)) .^ (1 / df);
        scale = sqrt((df + 1) / (1 - theta^2));
        x = (ratio - theta) .* scale;
        y = (1 ./ ratio - theta) .* scale;
        A = w .* tcdf(x, df + 1) + (1 - w) .* tcdf(y, df + 1);
    otherwise
        error("gofcopula:copula:InternalFamily", ...
            "Family '%s' is not an extreme-value family.", family);
end
end

function [A1, A2] = galambosDerivatives(w, theta)
if theta == 0
    A1 = zeros(size(w));
    A2 = zeros(size(w));
    return
end
% Work in logs: w^-theta overflows for large theta.
logW = log(w);
logV = log1p(-w);
la = -theta .* logW;
lb = -theta .* logV;
maximum = max(la, lb);
logS = maximum + log(exp(la - maximum) + exp(lb - maximum));
% First derivative: A1 = -S^(-1/theta-1) * (w^(-theta-1) - (1-w)^(-theta-1))
t1 = la - logW;
t2 = lb - logV;
m1 = max(t1, t2);
difference = exp(t1 - m1) - exp(t2 - m1);
A1 = -exp((-1/theta - 1) .* logS + m1) .* difference;
% Second derivative:
% A2 = (1+theta)*[ (w^(-t-2)+(1-w)^(-t-2))*S^(-1/t-1)
%                  - (w^(-t-1)-(1-w)^(-t-1))^2 * S^(-1/t-2) ]
s1 = la - 2 .* logW;
s2 = lb - 2 .* logV;
m2 = max(s1, s2);
sumTerm = exp((-1/theta - 1) .* logS + m2) .* (exp(s1 - m2) + exp(s2 - m2));
squareTerm = exp((-1/theta - 2) .* logS + 2 .* m1) .* difference.^2;
A2 = (1 + theta) .* (sumTerm - squareTerm);
end

function [A1, A2] = huslerReissDerivatives(w, theta)
if theta == 0
    A1 = sign(w - 0.5);
    A2 = zeros(size(w));
    return
end
z = 0.5 .* theta .* log(w ./ (1 - w));
zPlus = 1/theta + z;
zMinus = 1/theta - z;
A1 = normcdf(zPlus) - normcdf(zMinus);
A2 = (theta ./ (2 .* w .* (1 - w))) .* (normpdf(zPlus) + normpdf(zMinus));
end

function [A1, A2] = tevDerivatives(w, theta, df)
nu = df;
r = (w ./ (1 - w)).^(1 / nu);
s = sqrt((nu + 1) ./ (1 - theta.^2));
x = (r - theta) .* s;
y = (1 ./ r - theta) .* s;
g = @(t) tpdf(t, nu + 1);
gPrime = @(t) -(nu + 2) .* t ./ ((nu + 1) + t.^2) .* g(t);

k = 1 ./ (nu .* w .* (1 - w));
rPrime = r .* k;
xPrime = s .* rPrime;
yPrime = -s .* rPrime ./ r.^2;

P = w .* r .* g(x);
Q = (1 - w) .* g(y) ./ r;
A1 = tcdf(x, nu + 1) - tcdf(y, nu + 1) + s .* k .* (P - Q);

kPrime = -k .* (1 - 2 .* w) ./ (w .* (1 - w));
PPrime = r .* g(x) + w .* rPrime .* g(x) + w .* r .* gPrime(x) .* xPrime;
QPrime = -g(y) ./ r + (1 - w) .* gPrime(y) .* yPrime ./ r ...
    - (1 - w) .* g(y) .* rPrime ./ r.^2;
A2 = g(x) .* xPrime - g(y) .* yPrime ...
    + s .* (kPrime .* (P - Q) + k .* (PPrime - QPrime));
end
