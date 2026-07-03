classdef serialDependenceTest < matlab.unittest.TestCase
    %SERIALDEPENDENCETEST Unit tests for the serial-dependence-robust GoF path.
    %   Covers the decorrelation-length estimator, the row thinner, and the
    %   runTestSerial wrapper. All fixtures are deterministic (Threefry
    %   streams or literal matrices); no bootstrap larger than M=4 is used.

    methods (TestClassSetup)
        function addToolbox(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(root, "src")));
        end
    end

    methods (Test)
        % ---- decorrelationLength ----------------------------------------
        function whiteNoiseGivesUnitInterval(testCase)
            x = randn(RandStream("Threefry", "Seed", 11), 2000, 3);
            [interval, info] = gofcopula.internal.resampling.decorrelationLength(x);
            testCase.verifyEqual(interval, 1);
            testCase.verifyTrue(all(info.iatPerColumn == 1));
        end

        function ar1IntervalExceedsOneAndIsMonotoneInPhi(testCase)
            n = 3000;
            [~, f5] = gofcopula.internal.resampling.decorrelationLength(arSeries(0.5, n, 101));
            [~, f8] = gofcopula.internal.resampling.decorrelationLength(arSeries(0.8, n, 102));
            [i9, f9] = gofcopula.internal.resampling.decorrelationLength(arSeries(0.9, n, 103));
            testCase.verifyGreaterThan(f8.iatPerColumn, f5.iatPerColumn);
            testCase.verifyGreaterThan(f9.iatPerColumn, f8.iatPerColumn);
            testCase.verifyGreaterThan(i9, 1);
        end

        function multiColumnUsesMaxAggregation(testCase)
            x = [arSeries(0.5, 300, 201), arSeries(0.9, 300, 202)];
            [interval, info] = gofcopula.internal.resampling.decorrelationLength(x);
            testCase.verifySize(info.iatPerColumn, [1 2]);
            testCase.verifyGreaterThan(info.iatPerColumn(2), info.iatPerColumn(1));
            testCase.verifyGreaterThan(interval, 1);
        end

        function capBranchLimitsInterval(testCase)
            [interval, info] = gofcopula.internal.resampling.decorrelationLength( ...
                arSeries(0.95, 120, 55));
            testCase.verifyEqual(interval, 2);
            testCase.verifyGreaterThanOrEqual(info.nEffective, 50);
        end

        function minRetainedOneRemovesCap(testCase)
            x = arSeries(0.9, 3000, 105);
            [interval, info] = gofcopula.internal.resampling.decorrelationLength( ...
                x, MinRetained=1);
            testCase.verifyEqual(interval, max(1, round(info.iatPerColumn)));
        end

        function maxLagTruncatesTheSum(testCase)
            x = arSeries(0.9, 3000, 104);
            [~, fDefault] = gofcopula.internal.resampling.decorrelationLength(x);
            [~, fShort] = gofcopula.internal.resampling.decorrelationLength(x, MaxLag=3);
            testCase.verifyLessThan(fShort.iatPerColumn, fDefault.iatPerColumn);
        end

        function constantColumnGivesUnitIat(testCase)
            x = [ones(300, 1), arSeries(0.9, 300, 203)];
            [~, info] = gofcopula.internal.resampling.decorrelationLength(x);
            testCase.verifyEqual(info.iatPerColumn(1), 1);
        end

        function nonFiniteInputRejected(testCase)
            testCase.verifyError( ...
                @() gofcopula.internal.resampling.decorrelationLength([1 2; NaN 3]), ...
                "MATLAB:validators:mustBeFinite");
        end

        function tooFewRowsRejected(testCase)
            testCase.verifyError( ...
                @() gofcopula.internal.resampling.decorrelationLength([1 2]), ...
                "gofcopula:Serial:TooFewRows");
        end

        % ---- thinToIndependence -----------------------------------------
        function intervalOneIsIdentity(testCase)
            x = reshape(1:20, 10, 2);
            [xThin, keepIdx] = gofcopula.internal.resampling.thinToIndependence(x, 1);
            testCase.verifyEqual(xThin, x);
            testCase.verifyEqual(keepIdx, (1:10).');
        end

        function intervalTwoKeepsEveryOtherRow(testCase)
            x = reshape(1:16, 8, 2);
            [xThin, keepIdx] = gofcopula.internal.resampling.thinToIndependence(x, 2);
            testCase.verifyEqual(keepIdx, [1; 3; 5; 7]);
            testCase.verifySize(xThin, [4 2]);
            testCase.verifyEqual(xThin, x([1 3 5 7], :));
        end

        function offsetShiftsStart(testCase)
            x = reshape(1:16, 8, 2);
            [~, keepIdx] = gofcopula.internal.resampling.thinToIndependence(x, 2, 2);
            testCase.verifyEqual(keepIdx, [2; 4; 6; 8]);
        end

        function offsetBeyondRowsErrors(testCase)
            testCase.verifyError( ...
                @() gofcopula.internal.resampling.thinToIndependence([1 2; 3 4], 1, 9), ...
                "gofcopula:Serial:Offset");
        end

        % ---- runTestSerial ----------------------------------------------
        function returnsGofResult(testCase)
            x = fixedUnitData();
            result = gofcopula.runTestSerial("gofCvM", "normal", x, ...
                M=4, Param=0.3, ParamEst=false, DFEst=false, Margins="none", Seed=17);
            testCase.verifyClass(result, "gofcopula.GofResult");
            testCase.verifyEqual(result.Tests.Name, "gofCvM");
        end

        function autoIntervalOnTinyDataMatchesRunTest(testCase)
            x = fixedUnitData();
            [serialResult, serial] = gofcopula.runTestSerial("gofCvM", "normal", x, ...
                M=4, Param=0.3, ParamEst=false, DFEst=false, Margins="none", Seed=17);
            direct = gofcopula.runTest("gofCvM", "normal", x, ...
                M=4, Param=0.3, ParamEst=false, DFEst=false, Margins="none", Seed=17);
            testCase.verifyEqual(serial.thinInterval, 1);
            testCase.verifyEqual(serialResult.Tests.PValue, direct.Tests.PValue, AbsTol=eps);
        end

        function serialStructHasExpectedFields(testCase)
            x = fixedUnitData();
            [~, serial] = gofcopula.runTestSerial("gofCvM", "normal", x, ...
                M=0, Param=0.3, ParamEst=false, DFEst=false, Margins="none");
            expected = ["method", "nObserved", "nThinned", "thinInterval", ...
                "offset", "iatPerColumn", "maxLag", "keepIndices"];
            testCase.verifyTrue(all(ismember(expected, string(fieldnames(serial)))));
            testCase.verifyEqual(serial.nObserved, 5);
        end

        function explicitThinIntervalThins(testCase)
            x = rand(RandStream("Threefry", "Seed", 7), 8, 2);
            [~, serial] = gofcopula.runTestSerial("gofCvM", "normal", x, ...
                M=0, Param=0.3, ParamEst=false, DFEst=false, Margins="none", ThinInterval=2);
            testCase.verifyEqual(serial.thinInterval, 2);
            testCase.verifyEqual(serial.nThinned, 4);
            testCase.verifyTrue(all(isnan(serial.iatPerColumn)));
        end

        function blockSizeTrimsThinnedRows(testCase)
            x = rand(RandStream("Threefry", "Seed", 8), 8, 2);
            [~, serial] = gofcopula.runTestSerial("gofCvM", "normal", x, ...
                M=0, Param=0.3, ParamEst=false, DFEst=false, Margins="none", ...
                ThinInterval=2, BlockSize=3);
            testCase.verifyEqual(serial.nThinned, 3);
        end

        function multiplierMethodErrors(testCase)
            x = fixedUnitData();
            testCase.verifyError(@() gofcopula.runTestSerial("gofCvM", "normal", x, ...
                Method="multiplier"), "gofcopula:Serial:NotImplemented");
        end

        function bootstrapIsSeedReproducible(testCase)
            x = fixedUnitData();
            a = gofcopula.runTestSerial("gofCvM", "normal", x, ...
                M=4, Param=0.3, ParamEst=false, DFEst=false, Margins="none", Seed=23);
            b = gofcopula.runTestSerial("gofCvM", "normal", x, ...
                M=4, Param=0.3, ParamEst=false, DFEst=false, Margins="none", Seed=23);
            testCase.verifyEqual(a.Tests.PValue, b.Tests.PValue, AbsTol=eps);
        end

        % ---- runTestSerial: Method="phase" ------------------------------
        function phaseSurrogatePreservesCorrelation(testCase)
            X = arGaussianMatrix(0.8, 500, 3, 401);
            stream = RandStream("Threefry", "Seed", 9);
            Xs = gofcopula.internal.resampling.phaseSurrogate(X, stream);
            testCase.verifyTrue(isreal(Xs));
            testCase.verifyLessThan(max(abs(corr(X) - corr(Xs)), [], "all"), 1e-10);
        end

        function phaseUsesAllRowsForGaussian(testCase)
            X = arGaussianMatrix(0.8, 200, 2, 402);
            [result, serial] = gofcopula.runTestSerial("gofCvM", "normal", X, ...
                Method="phase", M=9, Seed=5);
            testCase.verifyClass(result, "gofcopula.GofResult");
            testCase.verifyEqual(serial.nThinned, 200);
            testCase.verifyEqual(serial.thinInterval, 1);
            testCase.verifyTrue(all(isnan(serial.iatPerColumn)));
        end

        function phaseIsSeedReproducible(testCase)
            X = arGaussianMatrix(0.7, 200, 2, 403);
            a = gofcopula.runTestSerial("gofCvM", "normal", X, Method="phase", M=9, Seed=8);
            b = gofcopula.runTestSerial("gofCvM", "normal", X, Method="phase", M=9, Seed=8);
            testCase.verifyEqual(a.Tests.PValue, b.Tests.PValue, AbsTol=eps);
        end

        function phaseRejectsNonGaussianFamily(testCase)
            X = arGaussianMatrix(0.5, 100, 2, 404);
            testCase.verifyError(@() gofcopula.runTestSerial("gofCvM", "clayton", X, ...
                Method="phase", M=9), "gofcopula:Serial:PhaseRequiresGaussian");
        end

        function phaseRequiresRanksMargins(testCase)
            X = arGaussianMatrix(0.5, 100, 2, 405);
            testCase.verifyError(@() gofcopula.runTestSerial("gofCvM", "normal", X, ...
                Method="phase", Margins="none", M=9), "gofcopula:Serial:PhaseRequiresRanks");
        end

        function phaseRejectsSuppliedBootstrapSamples(testCase)
            X = arGaussianMatrix(0.5, 100, 2, 406);
            bs = zeros(100, 2, 3);
            testCase.verifyError(@() gofcopula.runTestSerial("gofCvM", "normal", X, ...
                Method="phase", M=3, BootstrapSamples=bs), "gofcopula:Serial:PhaseConflict");
        end
    end
end

function z = arSeries(phi, n, seed)
%ARSERIES Stationary AR(1) column, unit variance, ACF phi^|h|.
stream = RandStream("Threefry", "Seed", seed);
e = randn(stream, n + 200, 1);
z = filter(sqrt(1 - phi^2), [1 -phi], e);
z = z(201:end);
end

function x = fixedUnitData()
%FIXEDUNITDATA Small pseudo-observation matrix inside the unit square.
x = [0.1 0.2; 0.3 0.25; 0.5 0.6; 0.7 0.8; 0.9 0.75];
end

function X = arGaussianMatrix(phi, n, d, seed)
%ARGAUSSIANMATRIX Correlated AR(1) Gaussian columns (exchangeable rho=0.5).
stream = RandStream("Threefry", "Seed", seed);
e = randn(stream, n + 200, d);
z = filter(sqrt(1 - phi^2), [1 -phi], e);
z = z(201:end, :);
R = 0.5 * ones(d) + 0.5 * eye(d);
X = z * chol(R);
end
