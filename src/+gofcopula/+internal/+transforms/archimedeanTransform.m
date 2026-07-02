function Z = archimedeanTransform(U, model)
%ARCHIMEDEANTRANSFORM Hering--Hofert transform used by Archm tests.
% For j=1,...,d-1, Z_j=(sum_{k<=j} psi^-1(U_k) /
% sum_{k<=j+1} psi^-1(U_k))^j. The last coordinate is the empirical Kendall
% transform of C(U), divided by n+1 exactly as in gofCopula 0.4-3.

arguments
    U (:,:) double {mustBeReal,mustBeFinite}
    model
end

if isa(model,"gofcopula.CopulaModel")
    family=model.Family; theta=model.Theta;
elseif isstruct(model)
    family=string(model.Family); theta=model.Theta;
else
    error("gofcopula:transforms:InvalidModel", "Invalid copula model.");
end
if ~isscalar(theta)
    error("gofcopula:transforms:ScalarParameterRequired", ...
        "The Archimedean transform requires one scalar copula parameter.");
end

[n,d] = size(U);
T = gofcopula.internal.transforms.inverseGenerator(family,U,theta);
cumulative = cumsum(T,2);
Z = zeros(n,d);
for j = 1:d-1
    denominator = cumulative(:,j+1);
    ratio = cumulative(:,j) ./ denominator;
    ratio(denominator == 0) = 1;
    Z(:,j) = ratio.^j;
end
C = gofcopula.internal.statistics.modelCDF(model,U);
Z(:,d) = gofcopula.internal.statistics.empiricalCopula(C,C,1);
Z = min(max(Z,0),1);
end
