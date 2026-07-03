function logDensity = peCopulaLogPDF(U, R, beta, options)
%PECOPULALOGPDF Log-density of the power-exponential (PE) copula.
%   L = peCopulaLogPDF(U,R,BETA) evaluates log c(u) at the rows of U for the
%   PE copula with correlation matrix R and shape BETA. With
%     z = Q1(u)  (the true PE marginal quantile) and g_d, g1 the joint and
%     marginal generators,
%       c(u) = |R|^(-1/2) g_d(z' R^{-1} z) / prod_j g1(z_j^2).
%   At BETA = 1 this equals the Gaussian copula density. Pass precomputed
%   marginal transforms through the Transforms option to avoid rebuilding the
%   generator grid (e.g. within an optimizer or bootstrap).

arguments
    U (:,:) double
    R (:,:) double
    beta (1,1) double {mustBeReal, mustBeFinite, mustBePositive}
    options.Transforms = []
end

[n, d] = size(U);
if isempty(options.Transforms)
    tk = gofcopula.internal.elliptical.peMarginals(d, beta, SampleSize=n);
else
    tk = options.Transforms;
end

X = min(max(U, realmin(class(U))), 1 - eps(class(U)));
Z = zeros(n, d);
for j = 1:d
    Z(:,j) = tk.Quantile(X(:,j));
end

L = chol(R, "lower");
solved = L \ Z.';
s = sum(solved.^2, 1).';                      % z' R^{-1} z per row
logJoint = -sum(log(diag(L))) + tk.LogConstant - 0.5 * s.^beta;
logMarginal = sum(tk.LogMarginalGenerator(Z.^2), 2);
logDensity = logJoint - logMarginal;
end
