function statistic = piosTnStatistic(U,model,options)
%PIOSTNSTATISTIC Leave-block-out pseudo-likelihood statistic.
% The R-compatible centering subtracts the number of fitted dependence
% parameters. Use options.Estimator(trainingU,model) to override estimation.

n=size(U,1);
m=fieldAny(options,["BlockLength","BlockSize"],1);
if ~(isscalar(m) && isfinite(m) && m==fix(m) && m>=1 && mod(n,m)==0)
    error("gofcopula:statistics:InvalidBlockLength", ...
        "BlockLength must be a positive divisor of the sample size.");
end
if isfield(options,"Estimator")
    estimator=options.Estimator;
else
    estimator=@defaultEstimator;
end
if ~isa(estimator,"function_handle")
    error("gofcopula:statistics:InvalidEstimator", ...
        "Estimator must be a function handle accepting (trainingU,model).");
end

fullLog=gofcopula.internal.statistics.modelPDF(model,U,true);
leaveLog=zeros(n,1); B=n/m;
for b=1:B
    held=(b-1)*m+(1:m);
    keep=true(n,1); keep(held)=false;
    fitted=estimator(U(keep,:),model);
    if iscell(fitted), fitted=fitted{1}; end
    leaveLog(held)=gofcopula.internal.statistics.modelPDF(fitted,U(held,:),true);
end
p=fieldOr(options,"ParameterCount",parameterCount(model));
statistic=sum(fullLog(:)-leaveLog)-p;
end

function fitted=defaultEstimator(trainingU,model)
if isa(model,"gofcopula.CopulaModel")
    % A fitted model is immutable and has EstimateTheta=false. Construct an
    % estimable copy so every omitted block really is refitted. As in the R
    % implementation, Student degrees of freedom remain fixed for Tn.
    source=gofcopula.CopulaModel(model.Family,Theta=model.Theta, ...
        DegreesOfFreedom=model.DegreesOfFreedom,EstimateTheta=true, ...
        EstimateDegreesOfFreedom=false,Dispersion=model.Dispersion, ...
        Rotation=model.Rotation);
else
    source=model;
end
[fitted,~]=gofcopula.internal.estimation.estimateModel(source,trainingU);
end

function p=parameterCount(model)
if isa(model,"gofcopula.CopulaModel")
    p=numel(model.Theta);
elseif isstruct(model) && isfield(model,"Theta")
    p=numel(model.Theta);
else
    p=1;
end
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
