classdef samplersTest < matlab.unittest.TestCase
    %SAMPLERSTEST Fast samplers against the bisection reference and exact laws.

    methods (TestClassSetup)
        function addToolbox(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(root, "src")));
        end
    end

    properties (TestParameter)
        differentialCase = struct( ...
            "clayton2", struct(Family="clayton", Theta=2.0, D=2), ...
            "claytonNegative", struct(Family="clayton", Theta=-0.5, D=2), ...
            "clayton3", struct(Family="clayton", Theta=2.0, D=3), ...
            "gumbel3", struct(Family="gumbel", Theta=1.8, D=3), ...
            "frank2", struct(Family="frank", Theta=4.0, D=2), ...
            "frankNegative", struct(Family="frank", Theta=-3.0, D=2), ...
            "frank3", struct(Family="frank", Theta=4.0, D=3), ...
            "joe2", struct(Family="joe", Theta=2.5, D=2), ...
            "amh2", struct(Family="amh", Theta=0.6, D=2), ...
            "fgm2", struct(Family="fgm", Theta=0.7, D=2), ...
            "galambos2", struct(Family="galambos", Theta=1.3, D=2), ...
            "plackett2", struct(Family="plackett", Theta=5.0, D=2));
        inverseCase = struct( ...
            "clayton", struct(Family="clayton", Thetas=[2.5, 0.7, -0.5]), ...
            "frank", struct(Family="frank", Thetas=[5, 0.5, -3]), ...
            "fgm", struct(Family="fgm", Thetas=[0.7, -0.9]));
    end

    methods (Test)
        function fastSamplerMatchesBisectionReference(testCase, differentialCase)
            family = differentialCase.Family;
            theta = differentialCase.Theta;
            d = differentialCase.D;
            new = gofcopula.internal.copulas.randomCore(family, 20000, d, ...
                theta, 4, "unstructured", 0, RandStream("Threefry","Seed",601));
            reference = gofcopula.internal.copulas.randomCoreBisection(family, ...
                4000, d, theta, 4, "unstructured", 0, RandStream("Threefry","Seed",701));
            % Margins must be uniform.
            for j = 1:d
                z = norminv(min(max(new(:,j), 1e-12), 1-1e-12));
                [~, p] = kstest(z);
                testCase.verifyGreaterThan(p, 1e-3, sprintf( ...
                    "%s margin %d fails uniformity (KS p=%.4f).", family, j, p));
            end
            % Dependence must match the exact reference sampler within
            % Monte Carlo tolerance (reference tau SE ~ 0.011 at n=4000).
            tauNew = corr(new, "Type", "Kendall");
            tauRef = corr(reference, "Type", "Kendall");
            testCase.verifyEqual(tauNew(1,2), tauRef(1,2), AbsTol=0.04);
            if d == 2
                grid = (0.1:0.1:0.9).';
                [G1, G2] = ndgrid(grid, grid);
                points = [G1(:), G2(:)];
                Cn = gofcopula.internal.statistics.empiricalCopula(new, points);
                Cr = gofcopula.internal.statistics.empiricalCopula(reference, points);
                testCase.verifyLessThan(mean(abs(Cn - Cr)), 0.012);
            end
        end

        function closedFormInversesAreExact(testCase, inverseCase)
            [u, w] = ndgrid(0.03:0.061:0.97, 0.03:0.061:0.97);
            for theta = inverseCase.Thetas
                switch inverseCase.Family
                    case "clayton"
                        v = max(1 + u.^(-theta) .* (w.^(-theta/(1+theta)) - 1), ...
                            realmin).^(-1/theta);
                    case "frank"
                        v = -log1p(w .* expm1(-theta) ./ ...
                            ((1-w) .* exp(-theta.*u) + w)) ./ theta;
                    case "fgm"
                        A = theta .* (2.*u - 1);
                        v = w;
                        ok = abs(A) >= 1e-9;
                        v(ok) = ((A(ok)-1) + sqrt((1-A(ok)).^2 + ...
                            4.*A(ok).*w(ok))) ./ (2.*A(ok));
                end
                h = gofcopula.internal.copulas.conditionalCDF( ...
                    inverseCase.Family, [u(:), v(:)], theta, 4, "unstructured");
                testCase.verifyEqual(h, w(:), "AbsTol", 1e-10, sprintf( ...
                    "%s inverse fails at theta=%.2f.", inverseCase.Family, theta));
            end
        end

        function frailtyDistributionsAreExact(testCase)
            n = 50000;
            % Positive stable: Laplace transform exp(-t^alpha).
            alpha = 1/1.8;
            V = gofcopula.internal.copulas.frailty("gumbel", n, 1.8, ...
                RandStream("Threefry","Seed",12345));
            for t = [0.5, 1, 2]
                transform = exp(-t .* V);
                tolerance = 4 * std(transform) / sqrt(n);
                testCase.verifyEqual(mean(transform), exp(-t^alpha), ...
                    AbsTol=tolerance);
            end
            % Logarithmic pmf.
            V = gofcopula.internal.copulas.frailty("frank", n, 3.0, ...
                RandStream("Threefry","Seed",22));
            p = -expm1(-3);
            kMax = 12;
            observed = histcounts(min(V, kMax+1), 0.5:1:kMax+1.5);
            pmf = p.^(1:kMax) ./ (-(1:kMax) .* log1p(-p));
            pmf(kMax+1) = 1 - sum(pmf);
            chi2 = sum((observed - n.*pmf).^2 ./ (n.*pmf));
            testCase.verifyLessThan(chi2, chi2inv(0.999, kMax));
            % Sibuya pmf.
            alphaJ = 1/2.5;
            V = gofcopula.internal.copulas.frailty("joe", n, 2.5, ...
                RandStream("Threefry","Seed",23));
            kMax = 15;
            observed = histcounts(min(V, kMax+1), 0.5:1:kMax+1.5);
            pmfS = arrayfun(@(k) (-1)^(k+1) * gamma(alphaJ+1) / ...
                (gamma(k+1) * gamma(alphaJ-k+1)), 1:kMax);
            pmfS(kMax+1) = 1 - sum(pmfS);
            chi2 = sum((observed - n.*pmfS).^2 ./ (n.*pmfS));
            testCase.verifyLessThan(chi2, chi2inv(0.999, kMax));
            % Gamma via its CDF.
            V = gofcopula.internal.copulas.frailty("clayton", n, 2.0, ...
                RandStream("Threefry","Seed",24));
            [~, p] = kstest(V, "CDF", [V, gammainc(V, 0.5, "lower")]);
            testCase.verifyGreaterThan(p, 1e-3);
        end

        function samplingIsSeedDeterministic(testCase)
            for family = ["clayton", "gumbel", "frank", "joe", "amh", "galambos"]
                a = gofcopula.copulaRandom(family, 100, defaultTheta(family), ...
                    Stream=RandStream("Threefry","Seed",5));
                b = gofcopula.copulaRandom(family, 100, defaultTheta(family), ...
                    Stream=RandStream("Threefry","Seed",5));
                testCase.verifyEqual(a, b);
            end
        end

        function samplingSpeedBudget(testCase)
            % Regression guard for the 2700x sampler speedup.
            timer = tic;
            gofcopula.copulaRandom("clayton", 10000, 2.0, ...
                Stream=RandStream("Threefry","Seed",1));
            gofcopula.copulaRandom("joe", 10000, 2.0, ...
                Stream=RandStream("Threefry","Seed",2));
            elapsed = toc(timer);
            testCase.verifyLessThan(elapsed, 0.5);
        end
    end
end

function theta = defaultTheta(family)
switch family
    case {"gumbel", "joe"}, theta = 1.8;
    case "amh", theta = 0.5;
    case "galambos", theta = 1.2;
    otherwise, theta = 2.0;
end
end
