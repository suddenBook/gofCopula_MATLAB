function family = normalizeFamily(family)
%NORMALIZEFAMILY Return the canonical lower-case copula family name.

family = lower(string(family));
family = replace(family, ["-", "_", " "], "");

switch family
    case {"gaussian", "gauss", "normal"}
        family = "normal";
    case {"student", "studentt", "tcopula", "t"}
        family = "t";
    case {"clayton", "gumbel", "frank", "joe", "amh", ...
            "galambos", "tawn", "tev", "fgm", "plackett"}
        % Already canonical.
    case {"huslerreiss", "hueslerreiss"}
        family = "huslerreiss";
    case {"powerexp", "pe", "powerexponential"}
        family = "powerexp";
    otherwise
        error("gofcopula:copula:UnknownFamily", ...
            "Unknown copula family '%s'.", family);
end
end
