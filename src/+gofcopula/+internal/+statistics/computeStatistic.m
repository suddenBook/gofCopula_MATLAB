function statistic = computeStatistic(testName,U,model,options)
%COMPUTESTATISTIC Dispatch a raw goodness-of-fit statistic.
%   T = COMPUTESTATISTIC(NAME,U,MODEL,OPTIONS) computes one statistic only;
%   bootstrap calibration belongs to the bootstrap subsystem. U must be an
%   n-by-d matrix in [0,1]. OPTIONS is a scalar structure.
%
% Supported NAME aliases include CvM/Sn, KS, KendallCvM/SnK,
% KendallKS/TnK, RosenblattSnB/SnB, RosenblattSnC/SnC,
% RosenblattGamma/AnGamma, RosenblattChisq/AnChisq, the four corresponding
% Archm variants, PIOSRn/Rn, PIOSTn/Tn, Kernel, and White.
%
% For transform-based names, U is transformed internally. Set
% OPTIONS.InputTransformed=true when passing an already transformed sample.

if nargin < 4 || isempty(options), options=struct(); end
if ~isstruct(options) || ~isscalar(options)
    error("gofcopula:statistics:InvalidOptions", ...
        "Options must be a scalar structure.");
end
validateSample(U);
if isa(testName,"function_handle")
    statistic=testName(U,model);
    validateStatistic(statistic);
    return
end
if ~(ischar(testName) || (isstring(testName) && isscalar(testName)))
    error("gofcopula:statistics:InvalidTestName", ...
        "TestName must be text or a function handle.");
end
name=lower(regexprep(string(testName),'[^a-zA-Z0-9]',''));
if startsWith(name,"gof"), name=extractAfter(name,3); end
inputTransformed=fieldOr(options,"InputTransformed",false);

switch name
    case {"cvm","sn"}
        statistic=gofcopula.internal.statistics.empiricalStatistics("cvm",U,model);
    case "ks"
        statistic=gofcopula.internal.statistics.empiricalStatistics("ks",U,model);
    case {"kendallcvm","snk"}
        reference=referenceCopulaValues(U,model,options);
        statistic=gofcopula.internal.statistics.kendallStatistics("cvm",U,reference);
    case {"kendallks","tnk"}
        reference=referenceCopulaValues(U,model,options);
        statistic=gofcopula.internal.statistics.kendallStatistics("ks",U,reference);
    case {"rosenblattsnb","snb"}
        Z=ordinaryTransform(U,model,inputTransformed);
        statistic=gofcopula.internal.statistics.independenceStatistics("snb",Z);
    case {"rosenblattsnc","snc"}
        Z=ordinaryTransform(U,model,inputTransformed);
        statistic=gofcopula.internal.statistics.independenceStatistics("snc",Z);
    case {"rosenblattgamma","angamma","gamma"}
        Z=ordinaryTransform(U,model,inputTransformed);
        statistic=gofcopula.internal.statistics.andersonDarlingStatistic("gamma",Z);
    case {"rosenblattchisq","anchisq","chisq"}
        Z=ordinaryTransform(U,model,inputTransformed);
        statistic=gofcopula.internal.statistics.andersonDarlingStatistic("chisq",Z);
    case {"archmsnb"}
        Z=archmTransform(U,model,inputTransformed);
        statistic=gofcopula.internal.statistics.independenceStatistics("snb",Z);
    case {"archmsnc"}
        Z=archmTransform(U,model,inputTransformed);
        statistic=gofcopula.internal.statistics.independenceStatistics("snc",Z);
    case {"archmgamma","archmangamma"}
        Z=archmTransform(U,model,inputTransformed);
        statistic=gofcopula.internal.statistics.andersonDarlingStatistic("gamma",Z);
    case {"archmchisq","archmanchisq"}
        Z=archmTransform(U,model,inputTransformed);
        statistic=gofcopula.internal.statistics.andersonDarlingStatistic("chisq",Z);
    case {"piosrn","rn"}
        [model,options]=piosNormalCompatibility(U,model,options,false);
        statistic=gofcopula.internal.statistics.piosRnStatistic(U,model,options);
    case {"piostn","tn"}
        [model,options]=piosNormalCompatibility(U,model,options,true);
        statistic=gofcopula.internal.statistics.piosTnStatistic(U,model,options);
    case "kernel"
        statistic=gofcopula.internal.statistics.kernelStatistic(U,model,options);
    case "white"
        statistic=gofcopula.internal.statistics.whiteStatistic(U,model,options);
    case "customtest"
        custom=fieldOr(options,"CustomTest",[]);
        if ~isa(custom,"function_handle")
            error("gofcopula:statistics:InvalidCustomTest", ...
                "Options.CustomTest must be a function handle.");
        end
        statistic=custom(U,model);
    otherwise
        error("gofcopula:statistics:UnsupportedTest", ...
            "Unsupported statistic '%s'.",testName);
end
validateStatistic(statistic);
end

function validateSample(U)
if ~(isnumeric(U) && ismatrix(U) && isreal(U) && all(isfinite(U),"all") && ...
        size(U,1)>0 && size(U,2)>1)
    error("gofcopula:statistics:InvalidSample", ...
        "U must be a nonempty finite real numeric matrix with at least two columns.");
end
if any(U<0 | U>1,"all")
    error("gofcopula:statistics:OutsideUnitCube", ...
        "All entries of U must lie in [0,1].");
end
end

function validateStatistic(value)
if ~(isnumeric(value) && isreal(value) && isscalar(value) && isfinite(value))
    error("gofcopula:statistics:InvalidStatistic", ...
        "The computed statistic must be a finite real scalar.");
end
end

function Z=ordinaryTransform(U,model,inputTransformed)
if inputTransformed
    Z=U;
else
    Z=gofcopula.internal.transforms.rosenblattTransform(U,model);
end
end

function Z=archmTransform(U,model,inputTransformed)
if inputTransformed
    Z=U;
else
    Z=gofcopula.internal.transforms.archimedeanTransform(U,model);
end
end

function reference=referenceCopulaValues(U,model,options)
d=size(U,2);
if isfield(options,"ReferenceCopulaValues")
    supplied=options.ReferenceCopulaValues;
    if ismatrix(supplied) && size(supplied,2)==d && d>1
        reference=gofcopula.internal.statistics.empiricalCopula(supplied,supplied);
    else
        reference=supplied(:);
    end
else
    count=fieldOr(options,"ReferenceSampleSize",10000);
    stream=fieldOr(options,"Stream",[]);
    sample=gofcopula.internal.statistics.modelRandom(model,count,stream,d);
    reference=gofcopula.internal.statistics.empiricalCopula(sample,sample);
end
if ~(isnumeric(reference) && isreal(reference) && all(isfinite(reference)) && ...
        ~isempty(reference) && all(reference>=0 & reference<=1))
    error("gofcopula:statistics:InvalidKendallReference", ...
        "Reference copula values must be a nonempty finite vector in [0,1].");
end
end

function [model,options]=piosNormalCompatibility(U,model,options,isTn)
%PIOSNORMALCOMPATIBILITY R parameterization of the normal-family PIOS tests.
% In NumericMode "rCompatible", R treats the Gaussian copula's PIOS
% information matrices over ALL d(d-1)/2 pairwise correlations regardless
% of the fitted dispersion (original R internal_Tstats.R:189-210), and the Tn
% leave-block refits use the Pearson correlation of the normal scores
% (original R internal_PIOS_normal.R:13-23). "corrected" mode keeps the fitted
% dispersion's own parameterization.
mode=string(fieldOr(options,"NumericMode","corrected"));
if mode ~= "rCompatible" || ~isa(model,"gofcopula.CopulaModel") || ...
        model.Family ~= "normal"
    return
end
d=size(U,2);
R=gofcopula.internal.copulas.correlationMatrix(model.Theta,d,model.Dispersion);
model=gofcopula.CopulaModel("normal",Theta=R(tril(true(d),-1)).', ...
    DegreesOfFreedom=model.DegreesOfFreedom,EstimateTheta=false, ...
    EstimateDegreesOfFreedom=false,Dispersion="unstructured", ...
    Rotation=model.Rotation);
if isTn
    options.ParameterCount=d*(d-1)/2;
    options.Estimator=@pearsonNormalRefit;
end
end

function fitted=pearsonNormalRefit(trainingU,~)
z=norminv(min(max(trainingU,eps),1-eps));
R=corr(z);
d=size(trainingU,2);
fitted=gofcopula.CopulaModel("normal",Theta=R(tril(true(d),-1)).', ...
    EstimateTheta=false,EstimateDegreesOfFreedom=false, ...
    Dispersion="unstructured");
end

function value=fieldOr(s,name,default)
if isfield(s,name), value=s.(name); else, value=default; end
end
