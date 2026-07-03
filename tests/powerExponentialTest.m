classdef powerExponentialTest < matlab.unittest.TestCase
    %POWEREXPONENTIALTEST Tests for the power-exponential (PE) copula family.
    %   The PE copula generalizes the Gaussian copula with a shape beta (stored in
    %   DegreesOfFreedom); beta = 1 recovers the Gaussian copula. These tests pin
    %   the reduction to Gaussian, the marginal transform, shape recovery, the
    %   Rosenblatt transform, and end-to-end bootstrap execution.

    methods (TestClassSetup)
        function addToolbox(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(root,"src")));
        end
    end

    methods (Test)
        function densityReducesToGaussianAtBetaOne(testCase)
            R = [1 0.5; 0.5 1]; rng(1); U = rand(2000,2);
            lc = gofcopula.copulaPDF("powerexp",U,0.5,DF=1,Dispersion="unstructured",Log=true);
            lg = log(copulapdf('Gaussian',U,R));
            testCase.verifyLessThan(max(abs(lc-lg)), 1e-3);
        end

        function cdfReducesToGaussianAtBetaOne(testCase)
            R = [1 0.5; 0.5 1]; rng(2); U = rand(300,2);
            C = gofcopula.copulaCDF("powerexp",U,0.5,DF=1,Dispersion="unstructured");
            testCase.verifyLessThan(max(abs(C-copulacdf('Gaussian',U,R))), 8e-3);
        end

        function samplerMatchesGaussianKendallTau(testCase)
            rng(3); U = gofcopula.copulaRandom("powerexp",100000,0.6, ...
                Dimension=2,DF=1,Dispersion="unstructured");
            tau = corr(U,'type','Kendall');
            testCase.verifyLessThan(abs(tau(1,2)-(2/pi)*asin(0.6)), 0.01);
        end

        function marginalIntegratesToOne(testCase)
            tk = gofcopula.internal.elliptical.peMarginals(2,0.6,SampleSize=5000);
            testCase.verifyLessThan(abs((tk.CDF(1e6)-tk.CDF(-1e6))-1), 1e-4);
        end

        function marginalGeneratorMatchesNormalAtBetaOne(testCase)
            tk = gofcopula.internal.elliptical.peMarginals(2,1.0,SampleSize=20000);
            t = [0 0.5 1 2 4 8]';
            % Generator is analytic; the sinh grid interpolates between nodes, so
            % off-node accuracy is ~1e-8 rather than machine precision.
            testCase.verifyLessThan(max(abs(tk.LogMarginalGenerator(t) ...
                - (-0.5*log(2*pi) - t/2))), 1e-6);
        end

        function copulaEntropyShapeTerm(testCase)
            testCase.verifyLessThan(abs(peKappa(1.0)), 1e-3);   % kappa(2,1) = 0
            testCase.verifyLessThan(peKappa(0.7), 0);           % kappa < 0 for beta != 1
            testCase.verifyLessThan(peKappa(1.6), 0);
        end

        function estimatorRecoversBeta(testCase)
            model = gofcopula.CopulaModel("powerexp",Theta=0.5,Dispersion="unstructured");
            rng(11); U = gofcopula.copulaRandom("powerexp",3000,0.5, ...
                Dimension=2,DF=0.7,Dispersion="unstructured");
            fitted = gofcopula.internal.estimation.estimateModel(model,U);
            testCase.verifyLessThan(abs(fitted.Beta-0.7), 0.15);
            testCase.verifyLessThan(abs(fitted.Theta-0.5), 0.05);
        end

        function defaultBetaIsOne(testCase)
            model = gofcopula.CopulaModel("powerexp",Theta=0.5);
            testCase.verifyEqual(model.Beta, 1);
            testCase.verifyEqual(model.DegreesOfFreedom, 1);
        end

        function rosenblattProducesUniform(testCase)
            rng(4); U = gofcopula.copulaRandom("powerexp",3000,0.5, ...
                Dimension=2,DF=0.7,Dispersion="unstructured");
            Z = gofcopula.rosenblatt("powerexp",U,0.5,DF=0.7,Dispersion="unstructured");
            [~,p1] = kstest(Z(:,1),CDF=makedist("Uniform"));
            [~,p2] = kstest(Z(:,2),CDF=makedist("Uniform"));
            testCase.verifyGreaterThan(p1, 0.01);
            testCase.verifyGreaterThan(p2, 0.01);
            testCase.verifyLessThan(abs(corr(Z(:,1),Z(:,2))), 0.06);
        end

        function capabilityTableListsPowerexp(testCase)
            T = gofcopula.CopulaTestTable();
            testCase.verifyTrue(ismember("powerexp", string(T.Properties.VariableNames)));
            supported = gofcopula.gofTest4Copula("powerexp",2);
            testCase.verifyTrue(all(ismember( ...
                ["gofPIOSTn","gofKendallCvM","gofRosenblattSnB","gofCvM","gofWhite"], supported)));
        end

        function bootstrapTestsRun(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            X = load(fullfile(root,"data","IndexReturns2D.mat")).IndexReturns2D;
            for t = ["gofKendallCvM","gofKernel","gofPIOSRn","gofCvM","gofRosenblattSnB"]
                fn = str2func("gofcopula."+t);
                r = fn("powerexp",X,M=9,Seed=1);
                testCase.verifyClass(r,"gofcopula.GofResult");
                testCase.verifyGreaterThanOrEqual(r.Tests.PValue(1),0);
                testCase.verifyLessThanOrEqual(r.Tests.PValue(1),1);
            end
        end

        function rosenblattUniformInThreeDimensions(testCase)
            % Exercises the d > 2 conditional path (p >= 2 marginalization).
            rho = 0.4; rng(5);
            U = gofcopula.copulaRandom("powerexp",1500,rho, ...
                Dimension=3,DF=0.8,Dispersion="exchangeable");
            Z = gofcopula.rosenblatt("powerexp",U,rho,DF=0.8,Dispersion="exchangeable");
            for k = 1:3
                [~,p] = kstest(Z(:,k),CDF=makedist("Uniform"));
                testCase.verifyGreaterThan(p, 0.005);
            end
        end

        function marginalRespectsMaxRadiusOption(testCase)
            tk = gofcopula.internal.elliptical.peMarginals(2,0.7, ...
                MaxRadius=15,QuadraturePoints=3000);
            testCase.verifyGreaterThan(tk.MaxSquaredRadius, 190);
            testCase.verifyLessThan(abs((tk.CDF(1e6)-tk.CDF(-1e6))-1), 1e-3);
        end

        function ellipticalToolkitMatchesGaussian(testCase)
            % Cross-validation of the ported ElliptCopulas generator toolkit
            % against standard-normal closed forms (beta = 1 elliptical case).
            grid = (0:0.01:60)';
            gd = (2*pi)^(-3/2)*exp(-grid/2);
            g1 = gofcopula.internal.elliptical.convertGdToG1(grid,gd,3);
            ok = ~isnan(g1);
            testCase.verifyLessThan( ...
                max(abs(g1(ok)-(2*pi)^(-1/2)*exp(-grid(ok)/2))), 1e-5);
            gm = gofcopula.internal.elliptical.convertGdToGm(grid,gd,3,1);
            testCase.verifyEqual(gm(ok), g1(ok), AbsTol=1e-12);
            Fg1 = gofcopula.internal.elliptical.convertG1ToFg1(grid,g1);
            testCase.verifyLessThan(abs(Fg1(1.0)-normcdf(1.0)), 5e-3);
            Qg1 = gofcopula.internal.elliptical.convertG1ToQg1(grid,g1);
            testCase.verifyLessThan(abs(Qg1(0.975)-norminv(0.975)), 1.5e-2);
        end
    end
end

function k = peKappa(beta)
% Copula entropy shape term kappa(2,beta) = C(2,beta) - 2 H1(2,beta).
tk = gofcopula.internal.elliptical.peMarginals(2,beta,SampleSize=20000);
xg = (0:0.003:sqrt(tk.MaxSquaredRadius)*0.99)';
lg1 = tk.LogMarginalGenerator(xg.^2); f1 = exp(lg1);
H1 = -2*trapz(xg, f1.*lg1);
C2 = 1/beta + log(2)/beta + log(pi) + gammaln(1+1/beta);
k = C2 - 2*H1;
end
