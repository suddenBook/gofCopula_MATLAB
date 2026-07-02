%% gofCopula for MATLAB
% Test whether multivariate observations are consistent with a parametric
% copula family. Rows are observations and columns are marginal variables.
% The source requires MATLAB R2025b, Statistics and Machine Learning Toolbox,
% and Optimization Toolbox. Add the repository's src folder to the MATLAB path.

%% Load a bundled dataset
% Locate repository data independently of the current working folder.

sourceRoot = fileparts(fileparts(which("gofcopula.CopulaModel")));
repositoryRoot = fileparts(sourceRoot);
sample = load(fullfile(repositoryRoot,"data","IndexReturns2D.mat"));
x = sample.IndexReturns2D;

figure
scatter(x(:,1),x(:,2),18,"filled")
xlabel(sample.VariableNames(1))
ylabel(sample.VariableNames(2))
title("European index log returns")
grid on

%% Inspect supported tests
% Family/test support depends on dimension. Query it before a large bootstrap.

available = gofcopula.gofTest4Copula("normal",size(x,2));
disp(available)

%% Run a goodness-of-fit test
% The default NumericMode="corrected" gives calibrated bootstrap p-values;
% pass NumericMode="rCompatible" only to reproduce R gofCopula 0.4-3 output
% (see Numerics.md for why the two differ). The small bootstrap count keeps
% this guide quick; use at least 1000 for inference.

result = gofcopula.gofCvM("normal",x,M=19,Seed=42);
disp(result)

%% Specify a copula model explicitly
% A model object is useful when parameters, estimation, dispersion, or rotation
% need to be controlled explicitly.

model = gofcopula.CopulaModel("clayton",Theta=1,EstimateTheta=true);
modelResult = gofcopula.gofco(model,x, ...
    Tests="gofKendallKS",M=19,Seed=42);
disp(modelResult)

%% Next steps
% See MigrationGuide.md for R-to-MATLAB option names, CapabilityMatrix.md for
% supported dimensions, and Numerics.md before comparing numerical results.
