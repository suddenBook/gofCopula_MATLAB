function results = gof(x, options)
%GOF Run compatible tests across one or more copula families.
%   M may be a scalar (shared by all tests) or a vector with one bootstrap
%   count per requested test, as in R. Failing tests degrade to NaN rows;
%   unsupported test/family combinations are skipped and summarized in one
%   warning.
arguments
    x {mustBeFloat,mustBeReal,mustBeFinite,mustBeMatrix}
    options.Priority {mustBeTextScalar} = "copula"
    options.Copulas string = strings(0,1)
    options.Tests string = strings(0,1)
    options.CustomTests cell = {}
    options.Param {mustBeNumeric,mustBeReal,mustBeFinite} = 0.5
    options.ParamEst (1,1) logical = true
    options.DF (1,1) {mustBeNumeric,mustBeReal,mustBeFinite,mustBePositive} = 4
    options.DFEst (1,1) logical = true
    options.Margins = "ranks"
    options.Flip {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger} = 0
    options.M {mustBeNumeric,mustBeInteger,mustBeNonnegative,mustBeVector} = 1000
    options.MJ (1,1) {mustBeNumeric,mustBeInteger,mustBePositive} = 100
    options.Dispersion {mustBeTextScalar} = "exchangeable"
    options.BlockSize (1,1) {mustBeNumeric,mustBeInteger,mustBePositive} = 1
    options.KernelScale (1,1) {mustBeNumeric,mustBePositive} = 0.5
    options.IntegrationNodes (1,1) {mustBeNumeric,mustBeInteger,mustBePositive} = 12
    options.Lower {mustBeNumeric,mustBeReal} = []
    options.Upper {mustBeNumeric,mustBeReal} = []
    options.Seed {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger} = []
    options.Processes (1,1) {mustBeNumeric,mustBeInteger,mustBePositive} = 1
    options.NumericMode {mustBeTextScalar} = "corrected"
end
priority = validatestring(options.Priority,{'copula','tests'});
allFamilies = gofcopula.internal.utilities.families();
allTests = string(gofcopula.CopulaTestTable().Properties.RowNames);
families = lower(options.Copulas(:));
gaussian = families == "gaussian";
if any(gaussian)
    warning("gofcopula:gof:GaussianAlias", ...
        "The pre-0.1-3 name 'gaussian' was replaced with 'normal'.");
    families(gaussian) = "normal";
end
tests = options.Tests(:);
if isempty(families) && isempty(tests)
    if strcmp(priority,"copula")
        families = allFamilies(:);
    else
        tests = allTests;
    end
end
if isempty(families)
    families = allFamilies(arrayfun(@(f)all(arrayfun(@(t) ...
        gofcopula.internal.utilities.maxDimension(t,f)>=size(x,2),tests)),allFamilies));
end
if isempty(tests)
    tests = allTests(arrayfun(@(t)all(arrayfun(@(f) ...
        gofcopula.internal.utilities.maxDimension(t,f)>=size(x,2),families)),allTests));
end
if isempty(families) || isempty(tests)
    error("gofcopula:gof:NoCombination", ...
        "No requested test/copula combination supports dimension %d.",size(x,2));
end
mustBeMember(families,allFamilies);
bootstrapCounts = resolveBootstrapCounts(options.M, numel(tests));
results = gofcopula.GofResult.empty(0,numel(families));
skipped = strings(0,1);
for f = 1:numel(families)
    rotation = options.Flip(min(f,numel(options.Flip)));
    rows = table(Size=[0,3],VariableTypes=["string","double","double"], ...
        VariableNames=["Name","PValue","Statistic"]);
    base = [];
    for t = 1:numel(tests)
        if gofcopula.internal.utilities.maxDimension(tests(t),families(f)) < size(x,2)
            skipped(end+1) = tests(t) + "/" + families(f); %#ok<AGROW>
            continue
        end
        try
            r = gofcopula.runTest(tests(t),families(f),x,Param=options.Param, ...
                ParamEst=options.ParamEst,DF=options.DF,DFEst=options.DFEst, ...
                Margins=options.Margins,Flip=rotation,M=bootstrapCounts(t), ...
                MJ=options.MJ, ...
                Dispersion=options.Dispersion,BlockSize=options.BlockSize, ...
                KernelScale=options.KernelScale,IntegrationNodes=options.IntegrationNodes, ...
                Lower=options.Lower,Upper=options.Upper,Seed=options.Seed, ...
                Processes=options.Processes,NumericMode=options.NumericMode);
        catch reason
            % Degrade like R's gofHybrid: record an NA row and continue.
            warning("gofcopula:gof:TestFailed", ...
                "%s failed for the %s copula: %s", ...
                tests(t), families(f), reason.message);
            rows = [rows; {tests(t), NaN, NaN}]; %#ok<AGROW>
            continue
        end
        rows = [rows;r.Tests]; %#ok<AGROW>
        base = r;
    end
    for c = 1:numel(options.CustomTests)
        fn = options.CustomTests{c};
        try
            r = gofcopula.runTest("gofCustomTest",families(f),x, ...
                Param=options.Param,ParamEst=options.ParamEst,DF=options.DF, ...
                DFEst=options.DFEst,Margins=options.Margins,Flip=rotation, ...
                M=bootstrapCounts(end),MJ=options.MJ, ...
                Dispersion=options.Dispersion,BlockSize=options.BlockSize, ...
                KernelScale=options.KernelScale,IntegrationNodes=options.IntegrationNodes, ...
                Lower=options.Lower,Upper=options.Upper, ...
                Seed=options.Seed,Processes=options.Processes, ...
                NumericMode=options.NumericMode,CustomTest=fn);
        catch reason
            warning("gofcopula:gof:TestFailed", ...
                "Custom test %s failed for the %s copula: %s", ...
                func2str(fn), families(f), reason.message);
            rows = [rows; {string(func2str(fn)), NaN, NaN}]; %#ok<AGROW>
            continue
        end
        r.Tests.Name = string(func2str(fn)); rows = [rows;r.Tests]; %#ok<AGROW>
        base = r;
    end
    if isempty(base), continue; end
    rows = [rows;gofcopula.internal.utilities.hybridRows(rows)];
    results(end+1) = gofcopula.GofResult(Method="Combined copula GOF analysis", ... %#ok<AGROW>
        Copula=base.Copula,Margins=base.Margins, ...
        MarginParameters=base.MarginParameters,Theta=base.Theta, ...
        DegreesOfFreedom=base.DegreesOfFreedom,Rotation=rotation,Tests=rows, ...
        NumericMode=base.NumericMode);
end
if ~isempty(skipped)
    warning("gofcopula:gof:SkippedCombinations", ...
        "Skipped unsupported combinations in dimension %d: %s", ...
        size(x,2), strjoin(skipped, ", "));
end
end

function counts = resolveBootstrapCounts(M, testCount)
if isscalar(M)
    counts = repmat(double(M), testCount, 1);
elseif numel(M) == testCount
    counts = double(M(:));
else
    error("gofcopula:gof:BootstrapCountMismatch", ...
        "M must be scalar or contain one entry per requested test.");
end
end
