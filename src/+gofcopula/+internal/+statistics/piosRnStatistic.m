function statistic = piosRnStatistic(U,model,options)
%PIOSRNSTATISTIC Information-matrix ratio statistic of Zhang et al. (2016).
% Analytic score/Hessian handles may be supplied as options.ScoreFunction and
% options.HessianFunction. Otherwise stable central differences are used.

if isfield(options,"ScoreFunction") && isfield(options,"HessianFunction")
    scores=options.ScoreFunction(U,model);
    hessians=options.HessianFunction(U,model);
else
    step=fieldOr(options,"DerivativeStep",eps^(1/4));
    [scores,hessians]=gofcopula.internal.statistics.numericalInformation(U,model,step);
end
p=size(scores,2); H=sum(hessians,3); J=scores.'*scores;
statistic=trace(-(H\J))-p;
end

function value=fieldOr(s,name,default)
if isfield(s,name), value=s.(name); else, value=default; end
end
