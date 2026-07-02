function result = gofPIOSTn(copula, x, options)
%GOFPIOSTN Copula goodness-of-fit test with parametric bootstrap.
arguments
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
end
args = namedargs2cell(options);
result = gofcopula.runTest("gofPIOSTn",copula,x,args{:});
end
