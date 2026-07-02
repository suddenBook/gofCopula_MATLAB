function result = gofco(model, x, options)
%GOFCO Run tests using a CopulaModel specification.
arguments
    model (1,1) gofcopula.CopulaModel
    x {mustBeFloat,mustBeReal,mustBeFinite,mustBeMatrix}
    options.Tests string = ["gofPIOSRn","gofKernel"]
    options.CustomTests cell = {}
    options.Margins = "ranks"
    options.M {mustBeNumeric,mustBeInteger,mustBeNonnegative} = 1000
    options.MJ (1,1) {mustBeNumeric,mustBeInteger,mustBePositive} = 100
    options.BlockSize (1,1) {mustBeNumeric,mustBeInteger,mustBePositive} = 1
    options.KernelScale (1,1) {mustBeNumeric,mustBePositive} = 0.5
    options.IntegrationNodes (1,1) {mustBeNumeric,mustBeInteger,mustBePositive} = 12
    options.Lower {mustBeNumeric,mustBeReal} = []
    options.Upper {mustBeNumeric,mustBeReal} = []
    options.Seed {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger} = []
    options.Processes (1,1) {mustBeNumeric,mustBeInteger,mustBePositive} = 1
    options.NumericMode {mustBeTextScalar} = "corrected"
end
result = gofcopula.gof(x,Copulas=model.Family,Tests=options.Tests, ...
    CustomTests=options.CustomTests, ...
    Param=model.Theta,ParamEst=model.EstimateTheta,DF=model.DegreesOfFreedom, ...
    DFEst=model.EstimateDegreesOfFreedom,Margins=options.Margins, ...
    Flip=model.Rotation,M=options.M,MJ=options.MJ,Dispersion=model.Dispersion, ...
    BlockSize=options.BlockSize,KernelScale=options.KernelScale, ...
    IntegrationNodes=options.IntegrationNodes, ...
    Lower=options.Lower,Upper=options.Upper,Seed=options.Seed, ...
    Processes=options.Processes,NumericMode=options.NumericMode);
end
