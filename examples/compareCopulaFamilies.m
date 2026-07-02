%% Compare copula families with gof
% Apply one statistic to Gaussian and Clayton candidates. A fixed seed makes
% the bootstrap samples reproducible.

%% Load data

sourceRoot = fileparts(fileparts(which("gofcopula.CopulaModel")));
repositoryRoot = fileparts(sourceRoot);
sample = load(fullfile(repositoryRoot,"data","IndexReturns2D.mat"));
x = sample.IndexReturns2D;

%% Fit and test candidates

results = gofcopula.gof(x,Copulas=["normal" "clayton"], ...
    Tests="gofCvM",M=19,Seed=2718);
disp(results)
