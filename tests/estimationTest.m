classdef estimationTest < matlab.unittest.TestCase
    %ESTIMATIONTEST Estimation-engine properties beyond the R oracle.

    methods (TestClassSetup)
        function addToolbox(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(root, "src")));
        end
    end

    methods (Test)
        function trivariateArchimedeanUsesPseudoLikelihood(testCase)
            U = gofcopula.copulaRandom("clayton", 250, 2.0, Dimension=3, ...
                Stream=RandStream("Threefry","Seed",99));
            [fitted, method] = gofcopula.internal.estimation.estimateModel( ...
                gofcopula.CopulaModel("clayton"), U);
            testCase.verifyEqual(method, "mpl");
            testCase.verifyEqual(fitted.Theta, 2.0, AbsTol=0.5);
        end

        function boundsAreHonoredForArchimedeanFamilies(testCase)
            U = gofcopula.copulaRandom("clayton", 300, 2.5, ...
                Stream=RandStream("Threefry","Seed",100));
            [fitted, method] = gofcopula.internal.estimation.estimateModel( ...
                gofcopula.CopulaModel("clayton"), U, Upper=1.5);
            testCase.verifyEqual(method, "mpl");
            testCase.verifyLessThanOrEqual(fitted.Theta, 1.5 + 1e-9);
        end

        function negativeClaytonEstimateErrorsLikeR(testCase)
            stream = RandStream("Threefry","Seed",101);
            u1 = rand(stream, 150, 1);
            u2 = min(max(1 - u1 + 0.05 .* randn(stream, 150, 1), 0.01), 0.99);
            testCase.verifyError(@() gofcopula.internal.estimation.estimateModel( ...
                gofcopula.CopulaModel("clayton"), [u1, u2]), ...
                "gofcopula:Estimation:InvalidClayton");
        end

        function tawnEstimateStaysInValidDomain(testCase)
            U = gofcopula.copulaRandom("tawn", 300, 0.6, ...
                Stream=RandStream("Threefry","Seed",102));
            fitted = gofcopula.internal.estimation.estimateModel( ...
                gofcopula.CopulaModel("tawn"), U);
            testCase.verifyGreaterThanOrEqual(fitted.Theta, 0);
            testCase.verifyLessThanOrEqual(fitted.Theta, 1);
        end

        function tauInversionsRoundTrip(testCase)
            % invertTau inverts the sample's empirical Kendall tau; the
            % round trip pushes the returned theta back through the exact
            % forward mapping tau(theta) and must recover that empirical
            % tau (up to root-finder tolerance).
            weak = gofcopula.copulaRandom("normal", 400, 0.25, ...
                Stream=RandStream("Threefry","Seed",301)); % tau in AMH range
            strong = gofcopula.copulaRandom("normal", 400, 0.7, ...
                Stream=RandStream("Threefry","Seed",302));
            forward.amh = @(t) 1 - 2*((1-t)^2*log1p(-t) + t) / (3*t^2);
            forward.clayton = @(t) t / (t + 2);
            forward.gumbel = @(t) 1 - 1/t;
            samples = struct(amh=weak, clayton=strong, gumbel=strong);
            for family = ["amh", "clayton", "gumbel"]
                sample = samples.(family);
                tauMatrix = corr(sample, "Type", "Kendall");
                theta = gofcopula.internal.estimation.invertTau( ...
                    gofcopula.CopulaModel(family), sample);
                testCase.verifyEqual(forward.(family)(theta), ...
                    tauMatrix(1,2), "AbsTol", 1e-8, ...
                    sprintf("%s tau inversion failed.", family));
            end
            % Joe and Plackett have no closed forward form: confirm the
            % implied population tau via a large draw from the fitted theta.
            for family = ["joe", "plackett"]
                sample = strong;
                tauMatrix = corr(sample, "Type", "Kendall");
                theta = gofcopula.internal.estimation.invertTau( ...
                    gofcopula.CopulaModel(family), sample);
                big = gofcopula.copulaRandom(family, 60000, theta, ...
                    Stream=RandStream("Threefry","Seed",55));
                bigTau = corr(big, "Type", "Kendall");
                testCase.verifyEqual(bigTau(1,2), tauMatrix(1,2), AbsTol=0.02);
            end
        end

        function replicateFailureRetriesInsteadOfAborting(testCase)
            % A statistic erroring on some replicates must trigger the
            % deterministic redraw, not kill the bootstrap.
            clear("flakyStatistic"); % reset the persistent counter
            U = gofcopula.copulaRandom("clayton", 40, 2.0, ...
                Stream=RandStream("Threefry","Seed",103));
            result = gofcopula.runTest("gofCustomTest", "clayton", U, ...
                Margins="none", M=6, Seed=9, CustomTest=@flakyStatistic);
            testCase.verifyGreaterThanOrEqual(result.Tests.PValue, 0);
            testCase.verifyLessThanOrEqual(result.Tests.PValue, 1);
        end

        function gofDegradesFailedTestsToNaN(testCase)
            U = gofcopula.copulaRandom("clayton", 40, 2.0, ...
                Stream=RandStream("Threefry","Seed",104));
            boom = @(varargin) error("test:Boom", "always fails");
            warning("off", "gofcopula:gof:TestFailed");
            cleaner = onCleanup(@() warning("on", "gofcopula:gof:TestFailed"));
            results = gofcopula.gof(U, Copulas="clayton", ...
                Tests="gofCvM", CustomTests={boom}, M=4, Seed=5);
            rows = results(1).Tests;
            failedRow = rows(~startsWith(rows.Name, "hybrid") & ...
                isnan(rows.PValue), :);
            testCase.verifyEqual(height(failedRow), 1);
            hybrid = rows(startsWith(rows.Name, "hybrid"), :);
            testCase.verifyTrue(all(isnan(hybrid.PValue)));
        end
    end
end

function value = flakyStatistic(U, ~)
% Synthetic failing statistic: errors on every third call so the
% replicate retry logic is exercised.
persistent attempts
if isempty(attempts), attempts = 0; end
attempts = attempts + 1;
if mod(attempts, 3) == 2
    error("test:Flaky", "synthetic replicate failure");
end
value = mean(U(:));
end
