function result = runTest(testName, copula, x, options)
%RUNTEST Run one parametric-bootstrap copula goodness-of-fit test.
arguments
    testName {mustBeTextScalar}
    copula {mustBeTextScalar}
    x {mustBeFloat,mustBeReal,mustBeFinite,mustBeMatrix}
    options.Param {mustBeNumeric,mustBeReal,mustBeFinite} = 0.5
    options.ParamEst (1,1) logical = true
    options.DF (1,1) {mustBeNumeric,mustBeReal,mustBeFinite,mustBePositive} = 4
    options.DFEst (1,1) logical = true
    options.Margins = "ranks"
    options.Flip (1,1) {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger} = 0
    options.M (1,1) {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger,mustBeNonnegative} = 1000
    options.MJ (1,1) {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger,mustBePositive} = 100
    options.Dispersion {mustBeTextScalar} = "exchangeable"
    options.BlockSize (1,1) {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger,mustBePositive} = 1
    options.KernelScale (1,1) {mustBeNumeric,mustBeReal,mustBeFinite,mustBePositive} = 0.5
    options.IntegrationNodes (1,1) {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger,mustBePositive} = 12
    options.Lower {mustBeNumeric,mustBeReal} = []
    options.Upper {mustBeNumeric,mustBeReal} = []
    options.Seed {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger} = []
    options.Processes (1,1) {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger,mustBePositive} = 1
    options.NumericMode {mustBeTextScalar} = "corrected"
    options.CustomTest = []
    options.ModelSample {mustBeNumeric,mustBeReal,mustBeFinite} = []
    options.BootstrapSamples {mustBeNumeric,mustBeReal,mustBeFinite} = []
end
testName = string(testName); family = lower(string(copula));
if family == "gaussian", family = "normal"; end
if size(x,2) < 2
    error("gofcopula:Input:Dimension", "At least two data columns are required.");
end
maxD = gofcopula.internal.utilities.maxDimension(testName,family);
if maxD == 0 || size(x,2) > maxD
    error("gofcopula:UnsupportedCombination", ...
        "%s is not supported for %s in dimension %d.",testName,family,size(x,2));
end
mode = validatestring(options.NumericMode,{'rCompatible','corrected'});
mustBeMember(options.Flip,[0,90,180,270]);
if options.BlockSize > size(x,1) || mod(size(x,1),options.BlockSize) ~= 0
    error("gofcopula:PIOS:BlockSize", ...
        "BlockSize must divide the number of observations.");
end
model = gofcopula.CopulaModel(family,Theta=options.Param, ...
    DegreesOfFreedom=options.DF,EstimateTheta=options.ParamEst, ...
    EstimateDegreesOfFreedom=options.DFEst, ...
    Dispersion=options.Dispersion,Rotation=0);
margins = string(options.Margins);
[u,marginParameters] = gofcopula.internal.estimation.transformMargins(x,margins);
u = gofcopula.internal.estimation.rotateData(u,options.Flip);
[fitted,fitMethod] = gofcopula.internal.estimation.estimateModel( ...
    model,u,Lower=options.Lower,Upper=options.Upper);
% R modifies the fitted Student df in place for some tests; the adjusted
% model is reported, simulated from, and evaluated (rCompatible only).
fitted = gofcopula.internal.bootstrap.adjustDegreesOfFreedom(testName,fitted,string(mode));
statOptions = struct(MJ=options.MJ,BlockSize=options.BlockSize, ...
    KernelScale=options.KernelScale,IntegrationNodes=options.IntegrationNodes, ...
    NumericMode=string(mode),CustomTest=options.CustomTest, ...
    ModelSample=options.ModelSample);
if isempty(options.ModelSample)
    statOptions=rmfield(statOptions,"ModelSample");
end
[pValue,statistic] = gofcopula.internal.bootstrap.parametricBootstrap( ...
    testName,u,model,fitted,options.M,statOptions,Seed=options.Seed, ...
    Processes=options.Processes,Lower=options.Lower,Upper=options.Upper, ...
    NumericMode=mode,BootstrapSamples=options.BootstrapSamples, ...
    Margins=margins,MarginParameters=marginParameters);
tests = table(testName,pValue,statistic, ...
    VariableNames=["Name","PValue","Statistic"]);
result = gofcopula.GofResult(Method=sprintf( ...
    "Parametric bootstrap %s test (%s fit)",testName,fitMethod), ...
    Copula=family,Margins=margins,MarginParameters=marginParameters, ...
    Theta=fitted.Theta,DegreesOfFreedom=fitted.DegreesOfFreedom, ...
    Rotation=options.Flip,Tests=tests,NumericMode=mode);
end
