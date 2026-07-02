classdef (TestTags = {'Slow'}) calibrationTest < matlab.unittest.TestCase
    %CALIBRATIONTEST Statistical calibration of the bootstrap p-values.
    %   Under a correctly specified null, corrected-mode p-values must be
    %   (approximately) uniform on (0,1), and the tests must have power
    %   against a misspecified alternative. rCompatible mode intentionally
    %   reproduces the R package's raw-replicate bootstrap, whose null
    %   p-values concentrate near one with rank margins; that behavior is
    %   locked by a regression test so it cannot change silently.

    properties (Constant)
        SampleSize = 200
        BootstrapM = 99
    end

    methods (TestClassSetup)
        function addToolbox(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(root, "src")));
        end
    end

    methods (Test)
        function correctedNullPValuesAreUniform(testCase)
            % Normal-copula data tested against the normal family: the
            % p-value sample must be consistent with U(0,1).
            replications = 150;
            p = zeros(replications, 2);
            for k = 1:replications
                data = gofcopula.copulaRandom("normal", ...
                    calibrationTest.SampleSize, 0.6, ...
                    Stream=RandStream("Threefry", "Seed", 40000+k));
                cvm = gofcopula.runTest("gofCvM", "normal", data, ...
                    Margins="ranks", M=calibrationTest.BootstrapM, ...
                    Seed=41000+k, NumericMode="corrected");
                snb = gofcopula.runTest("gofRosenblattSnB", "normal", data, ...
                    Margins="ranks", M=calibrationTest.BootstrapM, ...
                    Seed=42000+k, NumericMode="corrected");
                p(k,:) = [cvm.Tests.PValue, snb.Tests.PValue];
            end
            for c = 1:2
                [~, ksP] = kstest(p(:,c), "CDF", makedist("Uniform"));
                testCase.verifyGreaterThan(ksP, 0.005, sprintf( ...
                    "Null p-values of column %d are not uniform " + ...
                    "(KS p=%.4f, mean=%.3f).", c, ksP, mean(p(:,c))));
                rejection = mean(p(:,c) < 0.05);
                testCase.verifyGreaterThanOrEqual(rejection, 0.01);
                testCase.verifyLessThanOrEqual(rejection, 0.12);
            end
        end

        function correctedTestHasPowerAgainstClayton(testCase)
            % Clayton(theta=2, tau=0.5) data tested against the normal
            % family must be rejected most of the time.
            replications = 60;
            p = zeros(replications, 1);
            for k = 1:replications
                data = gofcopula.copulaRandom("clayton", ...
                    calibrationTest.SampleSize, 2.0, ...
                    Stream=RandStream("Threefry", "Seed", 43000+k));
                r = gofcopula.runTest("gofCvM", "normal", data, ...
                    Margins="ranks", M=calibrationTest.BootstrapM, ...
                    Seed=44000+k, NumericMode="corrected");
                p(k) = r.Tests.PValue;
            end
            testCase.verifyGreaterThanOrEqual(mean(p < 0.05), 0.85, ...
                sprintf("Power %.2f below 0.85 (median p=%.3f).", ...
                mean(p < 0.05), median(p)));
        end

        function rCompatibleNullBehaviorIsLocked(testCase)
            % Documents the inherited R defect: with rank margins the raw
            % (never re-ranked) replicates make null p-values pile up near
            % one. This regression pin prevents silent behavior changes in
            % rCompatible mode.
            replications = 10;
            p = zeros(replications, 1);
            for k = 1:replications
                data = gofcopula.copulaRandom("normal", ...
                    calibrationTest.SampleSize, 0.6, ...
                    Stream=RandStream("Threefry", "Seed", 45000+k));
                r = gofcopula.runTest("gofCvM", "normal", data, ...
                    Margins="ranks", M=calibrationTest.BootstrapM, ...
                    Seed=46000+k, NumericMode="rCompatible");
                p(k) = r.Tests.PValue;
            end
            testCase.verifyGreaterThan(median(p), 0.9, sprintf( ...
                "rCompatible null p-values moved (median %.3f).", median(p)));
        end
    end
end
