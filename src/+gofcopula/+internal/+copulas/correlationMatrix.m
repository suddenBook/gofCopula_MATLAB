function R = correlationMatrix(theta, dimension, dispersion)
%CORRELATIONMATRIX Construct and validate an elliptical correlation matrix.

dispersion = lower(string(dispersion));
dispersion = replace(dispersion, ["-", "_", " "], "");
if ismatrix(theta) && size(theta, 1) == dimension && size(theta, 2) == dimension
    R = theta;
else
    switch dispersion
        case {"exchangeable", "ex", "compound", "compoundsymmetry"}
            if ~isscalar(theta)
                error("gofcopula:copula:InvalidParameter", ...
                    "Exchangeable dispersion requires a scalar correlation.");
            end
            R = theta + (1 - theta) * eye(dimension);
        case {"ar1", "autoregressive"}
            if ~isscalar(theta)
                error("gofcopula:copula:InvalidParameter", ...
                    "AR(1) dispersion requires a scalar correlation.");
            end
            R = toeplitz(theta .^ (0:dimension-1));
        case {"toeplitz", "toep"}
            if numel(theta) ~= dimension - 1
                error("gofcopula:copula:InvalidParameter", ...
                    "Toeplitz dispersion requires dimension-1 correlations.");
            end
            R = toeplitz([1; theta(:)]);
        case {"unstructured", "un"}
            count = dimension * (dimension - 1) / 2;
            if isscalar(theta) && dimension == 2
                R = [1 theta; theta 1];
            elseif numel(theta) == count
                R = eye(dimension);
                lowerMask = tril(true(dimension), -1);
                R(lowerMask) = theta(:);
                R = R + tril(R, -1).';
            else
                error("gofcopula:copula:InvalidParameter", ...
                    "Unstructured dispersion requires d(d-1)/2 correlations or a d-by-d matrix.");
            end
        otherwise
            error("gofcopula:copula:InvalidDispersion", ...
                "Unsupported dispersion structure '%s'.", dispersion);
    end
end

tol = 100 * eps(class(R));
if any(abs(diag(R) - 1) > tol) || norm(R - R.', "fro") > tol * dimension
    error("gofcopula:copula:InvalidCorrelation", ...
        "The dispersion matrix must be symmetric with a unit diagonal.");
end
[~, flag] = chol((R + R.') / 2);
if flag ~= 0
    error("gofcopula:copula:InvalidCorrelation", ...
        "The dispersion matrix must be positive definite.");
end
R = (R + R.') / 2;
end
