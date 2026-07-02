classdef evDerivativesTest < matlab.unittest.TestCase
    %EVDERIVATIVESTEST Analytic Pickands derivatives against finite differences.

    methods (TestClassSetup)
        function addToolbox(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(root, "src")));
        end
    end

    properties (TestParameter)
        derivativeCase = struct( ...
            "galambosSmall", struct(Family="galambos", Theta=0.8, DF=4), ...
            "galambosLarge", struct(Family="galambos", Theta=40, DF=4), ...
            "huslerReiss", struct(Family="huslerreiss", Theta=1.5, DF=4), ...
            "huslerReissWeak", struct(Family="huslerreiss", Theta=0.3, DF=4), ...
            "tawn", struct(Family="tawn", Theta=0.7, DF=4), ...
            "tev", struct(Family="tev", Theta=0.5, DF=5), ...
            "tevNegative", struct(Family="tev", Theta=-0.3, DF=2.5));
    end

    methods (Test)
        function analyticDerivativesMatchFiniteDifferences(testCase, derivativeCase)
            family = derivativeCase.Family;
            theta = derivativeCase.Theta;
            df = derivativeCase.DF;
            w = [0.01, 0.05:0.05:0.95, 0.99].';
            evaluate = @(x) gofcopula.internal.copulas.evDependence(family, x, theta, df);
            h = min(1e-4, 0.2 .* min(w, 1-w));
            fd1 = (evaluate(w-2*h) - 8*evaluate(w-h) + 8*evaluate(w+h) ...
                - evaluate(w+2*h)) ./ (12*h);
            fd2 = (-evaluate(w+2*h) + 16*evaluate(w+h) - 30*evaluate(w) ...
                + 16*evaluate(w-h) - evaluate(w-2*h)) ./ (12*h.^2);
            [A, A1, A2] = gofcopula.internal.copulas.evDependence(family, w, theta, df);
            testCase.verifyEqual(A1, fd1, AbsTol=1e-9);
            testCase.verifyEqual(A2, fd2, AbsTol=1e-6);
            % A convex, A(0)=A(1)=1, max(w,1-w) <= A <= 1.
            testCase.verifyGreaterThanOrEqual(A2, -1e-8); % tail roundoff at extreme theta
            testCase.verifyLessThanOrEqual(A, 1 + 1e-12);
            testCase.verifyGreaterThanOrEqual(A, max(w, 1-w) - 1e-12);
        end
    end
end
