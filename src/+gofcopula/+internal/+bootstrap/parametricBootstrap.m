function [pValue, observed] = parametricBootstrap(testName,u,sourceModel,fitted,M,statOptions,options)
%PARAMETRICBOOTSTRAP Composite-null bootstrap with replicate refitting.
%   Two numeric regimes (see doc/Numerics.md):
%     "corrected"   - every replicate passes through the SAME margins
%                     pipeline as the observed data (Genest & Remillard,
%                     2008), refitting honors the user's estimation
%                     request, and p = (count+1)/(M+1).
%     "rCompatible" - reproduces gofCopula 0.4-3 (R): replicate statistics
%                     are computed on the raw copula draws (never
%                     re-ranked), parameters are ALWAYS refitted
%                     (original R internal_bootstrap.R:136), Student df follows the
%                     R ceil/cap rules, and p = count/M.
arguments
    testName (1,1) string
    u {mustBeFloat,mustBeReal,mustBeFinite,mustBeMatrix}
    sourceModel (1,1) gofcopula.CopulaModel
    fitted (1,1) gofcopula.CopulaModel
    M (1,1) {mustBeNumeric,mustBeInteger,mustBeNonnegative}
    statOptions (1,1) struct
    options.Seed {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger} = []
    options.Processes (1,1) {mustBeNumeric,mustBeInteger,mustBePositive} = 1
    options.Lower {mustBeNumeric,mustBeReal} = []
    options.Upper {mustBeNumeric,mustBeReal} = []
    options.NumericMode {mustBeTextScalar} = "corrected"
    options.BootstrapSamples {mustBeNumeric,mustBeReal,mustBeFinite} = []
    options.Margins string = "none"
    options.MarginParameters cell = {}
end
mode = string(validatestring(options.NumericMode,{'rCompatible','corrected'}));
seeds = makeSeeds(options.Seed,M+1);
if ismember(testName,["gofKendallCvM","gofKendallKS"])
    referenceStream = RandStream("Threefry","Seed",seeds(1));
    referenceSample = gofcopula.internal.statistics.modelRandom( ...
        fitted,10000,referenceStream,size(u,2));
    statOptions.ReferenceCopulaValues = ...
        gofcopula.internal.statistics.empiricalCopula(referenceSample,referenceSample);
end
observedOptions = statOptions;
observedOptions.Stream = RandStream("Threefry","Seed",seeds(1));
observed = compute(testName,u,fitted,observedOptions);
if M == 0
    pValue = NaN; return
end
boot = zeros(M,1); n=size(u,1); d=size(u,2);
if ~isempty(options.BootstrapSamples) && ...
        ~isequal(size(options.BootstrapSamples),[n,d,M])
    error("gofcopula:Bootstrap:SampleSize", ...
        "BootstrapSamples must be an n-by-d-by-M numeric array.");
end
% Replicate estimation policy: R refits unconditionally; corrected mode
% refits only the quantities the user asked to estimate. Refits start at
% the fitted parameters, as fitCopula does in R.
wantTheta = mode == "rCompatible" || sourceModel.EstimateTheta;
wantDF = fitted.Family == "t" && sourceModel.EstimateDegreesOfFreedom;
if wantTheta || wantDF
    replicateSource = gofcopula.CopulaModel(fitted.Family,Theta=fitted.Theta, ...
        DegreesOfFreedom=fitted.DegreesOfFreedom,EstimateTheta=wantTheta, ...
        EstimateDegreesOfFreedom=wantDF,Dispersion=fitted.Dispersion, ...
        Rotation=fitted.Rotation);
else
    replicateSource = [];
end
margins = options.Margins; marginParameters = options.MarginParameters;
useParallel = options.Processes > 1 && license("test","Distrib_Computing_Toolbox");
if useParallel
    parfor (b = 1:M, options.Processes)
        boot(b) = oneReplicate(testName,n,d,replicateSource,fitted, ...
            statOptions,seeds(b+1),options.Lower,options.Upper, ...
            mode,margins,marginParameters, ...
            selectSample(options.BootstrapSamples,b));
    end
else
    for b = 1:M
        boot(b) = oneReplicate(testName,n,d,replicateSource,fitted, ...
            statOptions,seeds(b+1),options.Lower,options.Upper, ...
            mode,margins,marginParameters, ...
            selectSample(options.BootstrapSamples,b));
    end
end
count = sum(abs(boot) >= abs(observed));
if mode == "corrected"
    pValue = (count+1)/(M+1);
else
    pValue = count/M;
end
end

function statistic = oneReplicate(testName,n,d,replicateSource,fitted,statOptions,seed,lower,upper,mode,margins,marginParameters,simulated)
if ~isempty(simulated)
    % Injected samples are deterministic; a failure must surface directly.
    statistic = replicateAttempt(testName,n,d,replicateSource,fitted, ...
        statOptions,seed,lower,upper,mode,margins,marginParameters,1,simulated);
    return
end
% As in R (internal_bootstrap.R:112-153), a replicate whose estimation or
% statistic fails is redrawn deterministically from the next substream.
lastError = [];
for attempt = 1:100
    try
        statistic = replicateAttempt(testName,n,d,replicateSource,fitted, ...
            statOptions,seed,lower,upper,mode,margins,marginParameters,attempt,[]);
        return
    catch lastError %#ok<NASGU>
    end
end
rethrow(lastError);
end

function statistic = replicateAttempt(testName,n,d,replicateSource,fitted,statOptions,seed,lower,upper,mode,margins,marginParameters,attempt,simulated)
stream = RandStream("Threefry","Seed",seed);
stream.Substream = attempt;
if isempty(simulated)
    simulated = gofcopula.copulaRandom(fitted.Family,n,fitted.Theta, ...
        Dimension=d,DF=fitted.DegreesOfFreedom,Dispersion=fitted.Dispersion, ...
        Rotation=fitted.Rotation,Stream=stream);
end
if mode == "corrected"
    simulated = gofcopula.internal.bootstrap.replicatePipeline( ...
        simulated,margins,marginParameters);
end
if isempty(replicateSource)
    bootstrapModel = fitted;
else
    bootstrapModel = gofcopula.internal.estimation.estimateModel( ...
        replicateSource,simulated,Lower=lower,Upper=upper);
end
bootstrapModel = gofcopula.internal.bootstrap.adjustDegreesOfFreedom( ...
    testName,bootstrapModel,mode);
statOptions.Stream = stream;
statistic = compute(testName,simulated,bootstrapModel,statOptions);
end

function sample=selectSample(samples,index)
if isempty(samples), sample=[]; else, sample=samples(:,:,index); end
end

function statistic = compute(testName,u,model,options)
statistic = gofcopula.internal.statistics.computeStatistic(testName,u,model,options);
if ~isscalar(statistic) || ~isreal(statistic) || ~isfinite(statistic)
    error("gofcopula:Statistic:Invalid", ...
        "A test statistic must return one finite real scalar.");
end
end

function seeds=makeSeeds(seed,count)
if isempty(seed)
    stream=RandStream.getGlobalStream;
    seeds=randi(stream,intmax("int32"),count,1);
elseif isscalar(seed)
    stream=RandStream("Threefry","Seed",double(seed));
    seeds=randi(stream,intmax("int32"),count,1);
elseif numel(seed)==count
    seeds=double(seed(:));
else
    error("gofcopula:Seed:WrongCount", ...
        "Seed must be scalar or contain M+1 integers.");
end
end
