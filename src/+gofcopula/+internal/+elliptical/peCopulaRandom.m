function U = peCopulaRandom(n, R, beta, stream, options)
%PECOPULARANDOM Draw observations from the power-exponential (PE) copula.
%   U = peCopulaRandom(N,R,BETA,STREAM) returns an N-by-d matrix of PE-copula
%   pseudo-observations. Uses the elliptical stochastic representation
%     X = radius * (uniform-on-sphere) * chol(R)',   radius = T^(1/2beta),
%     T ~ Gamma(shape = d/2beta, scale = 2)   (Gomez et al. 1998),
%   then maps each margin through the true PE marginal CDF F1. All randomness
%   is drawn from STREAM ([] uses the global stream), mirroring randomCore.

arguments
    n (1,1) double {mustBeInteger, mustBePositive}
    R (:,:) double
    beta (1,1) double {mustBeReal, mustBeFinite, mustBePositive}
    stream = []
    options.Transforms = []
end

d = size(R, 1);
if isempty(options.Transforms)
    tk = gofcopula.internal.elliptical.peMarginals(d, beta, SampleSize=n);
else
    tk = options.Transforms;
end

directions = normalRandom(stream, n, d);
sphere = directions ./ vecnorm(directions, 2, 2);        % uniform on the unit sphere
% T = radius^(2beta) ~ Gamma(d/2beta, scale 2), sampled by inverse-CDF for stream control.
T = 2 .* gammaincinv(uniformRandom(stream, n, 1), d/(2*beta), "lower");
radius = T.^(1/(2*beta));
L = chol(R, "lower");
X = radius .* (sphere * L.');                             % PE(0,R,beta) latent draw
U = tk.CDF(X);                                            % true PE marginal CDF -> uniform
U = min(max(U, 0), 1);
end

function u = uniformRandom(stream, varargin)
if isempty(stream), u = rand(varargin{:}); else, u = rand(stream, varargin{:}); end
end

function z = normalRandom(stream, varargin)
if isempty(stream), z = randn(varargin{:}); else, z = randn(stream, varargin{:}); end
end
