function statistic = whiteStatistic(U,model,options)
%WHITESTATISTIC White information-matrix test with influence correction.
% The covariance uses d_i-D*Hbar^(-1)*s_i, where d_i is the half-vectorized
% information-matrix discrepancy, D is the derivative of mean(d_i) with
% respect to the fitted parameters, Hbar is the mean Hessian, and s_i is the
% score. Derivatives are numerical unless analytic handles/matrices are given.

numericMode=string(fieldOr(options,"NumericMode","corrected"));
if strcmpi(numericMode,"rCompatible")
    statistic=rCompatibleStatistic(U,model);
    return
end

n=size(U,1);
[scores,hessians]=informationAt(model);
[discrepancies,mask]=discrepancyRows(scores,hessians);
p=size(scores,2); q=size(discrepancies,2);
dbar=mean(discrepancies,1).';
Hbar=mean(hessians,3);

if isfield(options,"DerivativeMatrix")
    derivative=options.DerivativeMatrix;
    if ~isequal(size(derivative),[q,p]) || any(~isfinite(derivative),"all")
        error("gofcopula:statistics:InvalidDerivativeMatrix", ...
            "DerivativeMatrix must be a finite q-by-p matrix.");
    end
else
    theta=parameterVector(model);
    outerScale=fieldOr(options,"OuterDerivativeStep",eps^(1/4));
    h=outerScale*(1+abs(theta));
    derivative=zeros(q,p);
    for a=1:p
        plus=theta; minus=theta;
        plus(a)=plus(a)+h(a); minus(a)=minus(a)-h(a);
        [sp,hp]=informationAt(modelWithParameters(model,plus));
        [sm,hm]=informationAt(modelWithParameters(model,minus));
        dp=meanDiscrepancy(sp,hp,mask);
        dm=meanDiscrepancy(sm,hm,mask);
        derivative(:,a)=(dp-dm)/(2*h(a));
    end
end

% Solve rather than explicitly invert Hbar. Each column corresponds to an
% observation's influence on the fitted parameter.
influence=(Hbar\scores.').';
adjusted=discrepancies-influence*derivative.';
V=cov(adjusted,1);
statistic=n*(dbar.'*pinv(V)*dbar);

    function [s,hessian]=informationAt(candidateModel)
        if isfield(options,"ScoreFunction") && isfield(options,"HessianFunction")
            s=options.ScoreFunction(U,candidateModel);
            hessian=options.HessianFunction(U,candidateModel);
        else
            step=fieldOr(options,"DerivativeStep",eps^(1/4));
            if isStudent(candidateModel)
                [s,hessian]=studentInformation(U,candidateModel,step);
            else
                [s,hessian]=gofcopula.internal.statistics.numericalInformation( ...
                    U,candidateModel,step);
            end
        end
        if size(s,1) ~= n || size(hessian,3) ~= n || ...
                size(hessian,1) ~= size(s,2) || size(hessian,2) ~= size(s,2)
            error("gofcopula:statistics:InvalidInformationDerivatives", ...
                "Score and Hessian outputs have inconsistent dimensions.");
        end
    end
end

function [rows,mask]=discrepancyRows(scores,hessians)
n=size(scores,1); p=size(scores,2); mask=tril(true(p));
rows=zeros(n,p*(p+1)/2);
for i=1:n
    discrepancy=hessians(:,:,i)+scores(i,:).'*scores(i,:);
    rows(i,:)=discrepancy(mask).';
end
end

function meanValue=meanDiscrepancy(scores,hessians,mask)
n=size(scores,1); total=zeros(nnz(mask),1);
for i=1:n
    discrepancy=hessians(:,:,i)+scores(i,:).'*scores(i,:);
    total=total+discrepancy(mask);
end
meanValue=total/n;
end

function theta=parameterVector(model)
if isa(model,"gofcopula.CopulaModel")
    theta=double(model.Theta(:).');
    if ismember(model.Family,["t","powerexp"]), theta=[theta,model.DegreesOfFreedom]; end
elseif isstruct(model) && isfield(model,"Theta")
    theta=double(model.Theta(:).');
    if isfield(model,'Family') && ismember(string(model.Family),["t","powerexp"])
        theta=[theta,model.DegreesOfFreedom];
    end
else
    error("gofcopula:statistics:MissingWhiteDerivative", ...
        "Provide a parameterized model or options.DerivativeMatrix.");
end
end

function statistic=rCompatibleStatistic(U,model)
% Reproduce VineCopula::BiCopGofTest(method="white",B=0).
n=size(U,1);
if isStudent(model)
    [scores,hessians]=studentInformation(U,model,1e-4);
    discrepancies=discrepancyRows(scores,hessians);
    dbar=mean(discrepancies,1).';
    V=(discrepancies.'*discrepancies)/n;
    statistic=n*(dbar.'*pinv(V)*dbar);
    return
end

% VineCopula re-estimates every non-t family even when par is supplied.
if isa(model,"gofcopula.CopulaModel")
    source=gofcopula.CopulaModel(model.Family,Theta=model.Theta, ...
        DegreesOfFreedom=model.DegreesOfFreedom,EstimateTheta=true, ...
        EstimateDegreesOfFreedom=false,Dispersion=model.Dispersion, ...
        Rotation=model.Rotation);
    fitted=gofcopula.internal.estimation.estimateModel(source,U);
else
    fitted=model;
end
theta=parameterVector(fitted); theta=theta(1); h=1e-4;
density=@(value)gofcopula.internal.statistics.modelPDF( ...
    modelWithParameters(fitted,value),U,false);
b=density(theta);
first=(density(theta+h)-density(theta-h))/(2*h);
second=(density(theta+h)-2*b+density(theta-h))/h^2;
d=second./b;
scoreAt=@(value)(density(value+h)-density(value-h))./(2*h)./density(value);
gradD=(mean(scoreAt(theta+h))-mean(scoreAt(theta-h)))/(2*h);
H=mean(-first./b.^2+d);
adjusted=d-(gradD/H)*(first./b);
V=mean(adjusted.^2);
statistic=n*mean(d)^2/V;
end

function [scores,hessians]=studentInformation(U,model,stepScale)
theta=parameterVector(model); p=numel(theta); n=size(U,1);
h=stepScale*(1+abs(theta)); base=logDensity(theta);
scores=zeros(n,p); hessians=zeros(p,p,n);
for a=1:p
    plus=theta; minus=theta; plus(a)=plus(a)+h(a); minus(a)=minus(a)-h(a);
    fp=logDensity(plus); fm=logDensity(minus);
    scores(:,a)=(fp-fm)/(2*h(a));
    hessians(a,a,:)=reshape((fp-2*base+fm)/h(a)^2,1,1,n);
end
for a=1:p
    for b=a+1:p
        pp=theta; pm=theta; mp=theta; mm=theta;
        pp([a,b])=pp([a,b])+[h(a),h(b)];
        pm(a)=pm(a)+h(a); pm(b)=pm(b)-h(b);
        mp(a)=mp(a)-h(a); mp(b)=mp(b)+h(b);
        mm([a,b])=mm([a,b])-[h(a),h(b)];
        mixed=(logDensity(pp)-logDensity(pm)-logDensity(mp)+logDensity(mm)) ...
            /(4*h(a)*h(b));
        hessians(a,b,:)=reshape(mixed,1,1,n);
        hessians(b,a,:)=reshape(mixed,1,1,n);
    end
end
    function values=logDensity(parameters)
        candidate=modelWithParameters(model,parameters);
        values=gofcopula.internal.statistics.modelPDF(candidate,U,true);
    end
end

function candidate=modelWithParameters(model,parameters)
if isa(model,"gofcopula.CopulaModel")
    df=model.DegreesOfFreedom; theta=parameters;
    if ismember(model.Family,["t","powerexp"]), theta=parameters(1:end-1); df=parameters(end); end
    candidate=gofcopula.CopulaModel(model.Family,Theta=theta, ...
        DegreesOfFreedom=df,EstimateTheta=false, ...
        EstimateDegreesOfFreedom=false,Dispersion=model.Dispersion, ...
        Rotation=model.Rotation);
else
    candidate=model; candidate.Theta=parameters;
end
end

function tf=isStudent(model)
tf=(isa(model,"gofcopula.CopulaModel") && model.Family == "t") || ...
    (isstruct(model) && isfield(model,'Family') && string(model.Family) == "t");
end

function value=fieldOr(s,name,default)
if isfield(s,name), value=s.(name); else, value=default; end
end
