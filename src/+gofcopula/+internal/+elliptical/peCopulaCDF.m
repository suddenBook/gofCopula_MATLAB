function C = peCopulaCDF(U, R, beta, options)
%PECOPULACDF Monte-Carlo distribution function of the power-exponential copula.
%   C = peCopulaCDF(U,R,BETA) returns C(u)=P(V<=u) at the rows of U, where V is a
%   PE copula with correlation R and shape BETA. The PE copula has no closed-form
%   CDF, so C is estimated from a large PE-copula reference sample. The reference
%   is drawn from a fixed-seed stream, so C is DETERMINISTIC given (R,BETA) and
%   the goodness-of-fit statistic is reproducible. At BETA=1 it matches the
%   Gaussian copula CDF.
%
%   Options:
%     MonteCarloSize - reference sample size (default 2e5; larger = more accurate).
%     Seed           - reference-stream seed (fixed so the CDF is reproducible).

arguments
    U (:,:) double
    R (:,:) double
    beta (1,1) double {mustBeReal, mustBeFinite, mustBePositive}
    options.MonteCarloSize (1,1) double {mustBeInteger, mustBePositive} = 500000
    options.Seed (1,1) double = 20240101
end

stream = RandStream("Threefry", Seed=options.Seed);
reference = gofcopula.internal.elliptical.peCopulaRandom( ...
    options.MonteCarloSize, R, beta, stream);

nq = size(U, 1);
C = zeros(nq, 1);
for i = 1:nq                                  % empirical multivariate CDF at each query
    C(i) = mean(all(reference <= U(i,:), 2));
end
end
