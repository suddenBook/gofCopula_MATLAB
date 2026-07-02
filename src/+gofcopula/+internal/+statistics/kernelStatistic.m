function statistic = kernelStatistic(U, model, options)
%KERNELSTATISTIC Integrated squared difference of biweight kernel densities.
% This reproduces gofCopula's implemented bivariate statistic. Despite the R
% front end accepting larger dimensions, its .ksmooth2 kernel uses two only;
% this implementation rejects d~=2 instead of silently discarding columns.

[n,d] = size(U);
if d ~= 2
    error("gofcopula:statistics:KernelDimension", ...
        "The kernel statistic is defined here for bivariate data only.");
end
k = fieldAny(options,["NodesIntegration","IntegrationNodes"],12);
MJ = fieldAny(options,["MonteCarloSize","MJ"],10000);
delta = fieldAny(options,["BandwidthScale","KernelScale"],1);
stream = fieldOr(options,"Stream",[]);
if isfield(options,"ModelSample")
    modelSample = options.ModelSample;
else
    modelSample = gofcopula.internal.statistics.modelRandom(model,MJ,stream,d);
end
if size(modelSample,2) ~= d
    error("gofcopula:statistics:DimensionMismatch", ...
        "The model sample must have the same dimension as U.");
end

L = chol(cov(U),"lower");
h = 2.6073 * n^(-1/6) * diag(L).' * delta;
if any(~isfinite(h) | h <= 0)
    error("gofcopula:statistics:DegenerateBandwidth", ...
        "Kernel bandwidth is nonpositive; input columns must vary.");
end
[x,w] = gaussLegendre01(k);
statistic = 0;
for i = 1:k
    for j = 1:k
        point = [x(i),x(j)];
        difference = biweight(point,U,h) - biweight(point,modelSample,h);
        statistic = statistic + w(i)*w(j)*difference^2;
    end
end
end

function y = biweight(point,sample,h)
q = (sample-point)./h;
inside = all(abs(q)<=1,2);
kernel = prod((1-q(inside,:).^2).^2,2);
y = (15/16)^2 * sum(kernel) / (size(sample,1)*prod(h));
end

function [nodes,weights] = gaussLegendre01(k)
if ~(isscalar(k) && isfinite(k) && k == fix(k) && k >= 1)
    error("gofcopula:statistics:InvalidQuadratureOrder", ...
        "NodesIntegration must be a positive integer.");
end
if k == 1
    nodes=0.5; weights=1; return
end
index=(1:k-1).'; beta=index./sqrt(4*index.^2-1);
[V,D]=eig(diag(beta,1)+diag(beta,-1));
[raw,order]=sort(diag(D)); V=V(:,order);
nodes=(raw+1)/2; weights=(2*V(1,:).^2).'/2;
end

function value=fieldOr(s,name,default)
if isfield(s,name), value=s.(name); else, value=default; end
end

function value=fieldAny(s,names,default)
value=default;
for name=names
    if isfield(s,name), value=s.(name); return, end
end
end
