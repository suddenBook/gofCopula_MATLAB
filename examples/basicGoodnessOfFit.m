%% Test a Gaussian copula with gofCvM
% Load the bundled bivariate index returns and run a reproducible bootstrap
% Cramer-von Mises test.

%% Load data

sourceRoot = fileparts(fileparts(which("gofcopula.CopulaModel")));
repositoryRoot = fileparts(sourceRoot);
sample = load(fullfile(repositoryRoot,"data","IndexReturns2D.mat"));
x = sample.IndexReturns2D;

%% Run the test
% Use a larger M for final statistical inference.

result = gofcopula.gofCvM("normal",x,M=19,Seed=2026);
disp(result.Tests)

%% Visualize the observations

figure
scatter(x(:,1),x(:,2),18,"filled")
xlabel(sample.VariableNames(1))
ylabel(sample.VariableNames(2))
title("Data tested against a Gaussian copula")
grid on
