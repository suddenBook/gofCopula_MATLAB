function t = inverseGenerator(family, u, theta)
%INVERSEGENERATOR Inverse Archimedean generator psi^{-1}(u).
% Stable formulas are used at independence limits where practical.

family = lower(string(family));
switch family
    case "clayton"
        if abs(theta) < sqrt(eps)
            t = -log(u);
        else
            t = expm1(-theta .* log(u));
        end
    case "gumbel"
        t = (-log(u)).^theta;
    case "frank"
        if abs(theta) < sqrt(eps)
            t = -log(u);
        else
            t = -log(expm1(-theta.*u) ./ expm1(-theta));
        end
    case "joe"
        t = -log1p(-(1-u).^theta);
    case "amh"
        t = log(theta + (1-theta)./u);
    otherwise
        error("gofcopula:transforms:UnsupportedArchimedeanFamily", ...
            "The Hering--Hofert transform is not implemented for '%s'.", family);
end
end
