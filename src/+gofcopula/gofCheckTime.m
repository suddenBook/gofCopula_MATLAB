function estimate = gofCheckTime(copula, x, options)
%GOFCHECKTIME Estimate the wall-clock time of a bootstrap run.
%   Times each requested test at two small bootstrap sizes and
%   extrapolates linearly in M (R fits the same line over an M grid).
%   Kernel-test cost is measured at the requested MJ directly. Returns a
%   duration.
arguments
    copula {mustBeTextScalar}
    x {mustBeFloat,mustBeReal,mustBeFinite,mustBeMatrix}
    options.Tests string = "gofCvM"
    options.CustomTests cell = {}
    options.M (1,1) {mustBeNumeric,mustBeInteger,mustBePositive} = 1000
    options.MJ (1,1) {mustBeNumeric,mustBeInteger,mustBePositive} = 100
    options.Param {mustBeNumeric,mustBeReal,mustBeFinite} = 0.5
    options.ParamEst (1,1) logical = true
    options.DF (1,1) {mustBeNumeric,mustBeReal,mustBeFinite,mustBePositive} = 4
    options.DFEst (1,1) logical = true
    options.Margins = "ranks"
    options.Flip (1,1) {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger} = 0
    options.Dispersion {mustBeTextScalar} = "exchangeable"
    options.BlockSize (1,1) {mustBeNumeric,mustBeInteger,mustBePositive} = 1
    options.KernelScale (1,1) {mustBeNumeric,mustBePositive} = 0.5
    options.IntegrationNodes (1,1) {mustBeNumeric,mustBeInteger,mustBePositive} = 12
    options.Lower {mustBeNumeric,mustBeReal} = []
    options.Upper {mustBeNumeric,mustBeReal} = []
    options.Seed {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger} = 1
    options.Processes (1,1) {mustBeNumeric,mustBeInteger,mustBePositive} = 1
    options.NumericMode {mustBeTextScalar} = "corrected"
end
probes = [2, 6];
total = 0;
for test = options.Tests(:).'
    times = arrayfun(@(M) timeRun(test, copula, x, options, M, []), probes);
    total = total + extrapolate(times, probes, options.M);
end
for c = 1:numel(options.CustomTests)
    times = arrayfun(@(M) timeRun("gofCustomTest", copula, x, options, M, ...
        options.CustomTests{c}), probes);
    total = total + extrapolate(times, probes, options.M);
end
estimate = seconds(total);
end

function elapsed = timeRun(test, copula, x, options, M, custom)
started = tic;
gofcopula.runTest(test, copula, x, Param=options.Param, ...
    ParamEst=options.ParamEst, DF=options.DF, DFEst=options.DFEst, ...
    Margins=options.Margins, Flip=options.Flip, M=M, MJ=options.MJ, ...
    Dispersion=options.Dispersion, BlockSize=options.BlockSize, ...
    KernelScale=options.KernelScale, IntegrationNodes=options.IntegrationNodes, ...
    Lower=options.Lower, Upper=options.Upper, Seed=options.Seed, ...
    Processes=options.Processes, NumericMode=options.NumericMode, ...
    CustomTest=custom);
elapsed = toc(started);
end

function projected = extrapolate(times, probes, M)
perReplicate = max((times(2) - times(1)) / (probes(2) - probes(1)), 0);
fixed = max(times(1) - probes(1) * perReplicate, 0);
projected = fixed + perReplicate * double(M);
end
