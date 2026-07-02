function U = randomCoreBisection(family, n, dimension, theta, df, dispersion, rotation, stream)
%RANDOMCOREBISECTION Reference sampler by per-row conditional bisection.
%   This is the original (slow, exact to 2^-54) conditional-inversion
%   sampler. It is retained as the differential-testing oracle for the
%   fast samplers in randomCore; production code should not call it.

if ismember(family, ["normal", "t"])
    R = gofcopula.internal.copulas.correlationMatrix(theta, dimension, dispersion);
    L = chol(R, "lower");
    z = normalRandom(stream, n, dimension) * L.';
    if family == "normal"
        U = normcdf(z);
    else
        chiSquare = 2 .* gammaincinv(uniformRandom(stream,n,1), df/2, "lower");
        U = tcdf(z ./ sqrt(chiSquare ./ df), df);
    end
else
    targets = uniformRandom(stream, n, dimension);
    U = zeros(n, dimension);
    U(:,1) = targets(:,1);
    for j = 2:dimension
        for i = 1:n
            prefix = U(i,1:j);
            target = targets(i,j);
            low = 0;
            high = 1;
            for iteration = 1:54
                middle = (low + high) / 2;
                prefix(j) = middle;
                value = gofcopula.internal.copulas.conditionalCDF( ...
                    family, prefix, theta, df, dispersion);
                if value < target
                    low = middle;
                else
                    high = middle;
                end
            end
            U(i,j) = (low + high) / 2;
        end
    end
end

if rotation ~= 0
    switch rotation
        case 90
            U(:,1) = 1-U(:,1);
        case 180
            U = 1-U;
        case 270
            U(:,2) = 1-U(:,2);
    end
end
U = min(max(U, 0), 1);
end

function u = uniformRandom(stream, varargin)
if isempty(stream), u = rand(varargin{:}); else, u = rand(stream, varargin{:}); end
end

function z = normalRandom(stream, varargin)
if isempty(stream), z = randn(varargin{:}); else, z = randn(stream, varargin{:}); end
end
