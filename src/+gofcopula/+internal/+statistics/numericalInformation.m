function [scores,hessians] = numericalInformation(U, model, stepScale)
%NUMERICALINFORMATION Per-observation score and Hessian by central differences.
% Intended for PIOS Rn and the White approximation when analytic derivatives
% are unavailable. STEPScale defaults to eps^(1/4), balancing cancellation
% and second-derivative truncation error in double precision.

if nargin < 3 || isempty(stepScale), stepScale=eps^(1/4); end
theta = parameterVector(model);
p = numel(theta); n = size(U,1);
if p == 0
    error("gofcopula:statistics:NoParameters", ...
        "The information-matrix statistic requires fitted parameters.");
end
h = stepScale*(1+abs(theta));
base = logDensity(theta);
scores = zeros(n,p); hessians = zeros(p,p,n);
plus = cell(p,1); minus = cell(p,1);
for a=1:p
    tp=theta; tm=theta; tp(a)=tp(a)+h(a); tm(a)=tm(a)-h(a);
    plus{a}=logDensity(tp); minus{a}=logDensity(tm);
    scores(:,a)=(plus{a}-minus{a})/(2*h(a));
    hessians(a,a,:)=reshape((plus{a}-2*base+minus{a})/h(a)^2,1,1,n);
end
for a=1:p
    for b=a+1:p
        tpp=theta; tpm=theta; tmp=theta; tmm=theta;
        tpp([a b])=tpp([a b])+[h(a) h(b)];
        tpm(a)=tpm(a)+h(a); tpm(b)=tpm(b)-h(b);
        tmp(a)=tmp(a)-h(a); tmp(b)=tmp(b)+h(b);
        tmm([a b])=tmm([a b])-[h(a) h(b)];
        mixed=(logDensity(tpp)-logDensity(tpm)-logDensity(tmp)+logDensity(tmm)) ...
            /(4*h(a)*h(b));
        hessians(a,b,:)=reshape(mixed,1,1,n);
        hessians(b,a,:)=reshape(mixed,1,1,n);
    end
end

    function value=logDensity(candidate)
        candidateModel=gofcopula.internal.statistics.modelWithTheta(model,candidate);
        value=gofcopula.internal.statistics.modelPDF(candidateModel,U,true);
        value=value(:);
        if any(~isfinite(value))
            error("gofcopula:statistics:DerivativeOutsideDomain", ...
                "A finite-difference perturbation produced a nonfinite log density.");
        end
    end
end

function theta=parameterVector(model)
if isa(model,"gofcopula.CopulaModel")
    theta=double(model.Theta(:).');
elseif isstruct(model) && isfield(model,"Theta")
    theta=double(model.Theta(:).');
else
    theta=[];
end
end
