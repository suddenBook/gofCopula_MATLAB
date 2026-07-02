classdef copulaCoreTest < matlab.unittest.TestCase
    %COPULACORETEST Public-interface tests for copula numerical primitives.

    properties (TestParameter)
        familyCase = struct( ...
            "normal", struct(family="normal", theta=0.35, df=4), ...
            "student", struct(family="t", theta=-0.25, df=6), ...
            "clayton", struct(family="clayton", theta=1.4, df=4), ...
            "gumbel", struct(family="gumbel", theta=1.7, df=4), ...
            "frank", struct(family="frank", theta=2.2, df=4), ...
            "joe", struct(family="joe", theta=1.8, df=4), ...
            "amh", struct(family="amh", theta=0.45, df=4), ...
            "galambos", struct(family="galambos", theta=1.1, df=4), ...
            "huslerReiss", struct(family="huslerReiss", theta=1.2, df=4), ...
            "tawn", struct(family="tawn", theta=0.7, df=4), ...
            "tev", struct(family="tev", theta=0.3, df=5), ...
            "fgm", struct(family="fgm", theta=-0.6, df=4), ...
            "plackett", struct(family="plackett", theta=2.5, df=4));
    end

    methods (TestClassSetup)
        function addToolboxToPath(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(root,"src")));
        end
    end

    methods (Test)
        function cdfHasCopulaBoundsAndMargins(testCase, familyCase)
            U = [0.23 0.61; 0.37 1; 0 0.8; 1 1];
            C = gofcopula.copulaCDF(familyCase.family,U,familyCase.theta, ...
                DF=familyCase.df);
            testCase.verifyGreaterThanOrEqual(C, zeros(4,1));
            testCase.verifyLessThanOrEqual(C, min(U,[],2) + 2e-12);
            testCase.verifyEqual(C(2), U(2,1), AbsTol=2e-9);
            testCase.verifyEqual(C(3), 0);
            testCase.verifyEqual(C(4), 1, AbsTol=2e-9);
        end

        function pdfIsPositiveAndLogConsistent(testCase, familyCase)
            U = [0.19 0.31; 0.46 0.73; 0.81 0.58];
            density = gofcopula.copulaPDF(familyCase.family,U,familyCase.theta, ...
                DF=familyCase.df);
            logDensity = gofcopula.copulaPDF(familyCase.family,U,familyCase.theta, ...
                DF=familyCase.df,Log=true);
            testCase.verifyGreaterThan(density, zeros(3,1));
            testCase.verifyEqual(log(density), logDensity, AbsTol=2e-12);
        end

        function randomOutputAndRosenblattAreValid(testCase, familyCase)
            firstStream = RandStream("mt19937ar",Seed=9321);
            secondStream = RandStream("mt19937ar",Seed=9321);
            first = gofcopula.copulaRandom(familyCase.family,16,familyCase.theta, ...
                DF=familyCase.df,Stream=firstStream);
            second = gofcopula.copulaRandom(familyCase.family,16,familyCase.theta, ...
                DF=familyCase.df,Stream=secondStream);
            transformed = gofcopula.rosenblatt( ...
                familyCase.family,first,familyCase.theta,DF=familyCase.df);
            testCase.verifyEqual(first,second);
            testCase.verifySize(first,[16 2]);
            testCase.verifyGreaterThan(first,zeros(16,2));
            testCase.verifyLessThan(first,ones(16,2));
            testCase.verifyGreaterThanOrEqual(transformed,zeros(16,2));
            testCase.verifyLessThanOrEqual(transformed,ones(16,2));
            testCase.verifyEqual(transformed(:,1),first(:,1));
        end

        function gaussianQuadrantProbability(testCase)
            rho = 0.6;
            actual = gofcopula.copulaCDF("normal",[0.5 0.5],rho);
            expected = 0.25 + asin(rho)/(2*pi);
            testCase.verifyEqual(actual,expected,AbsTol=2e-10);
        end

        function studentQuadrantProbability(testCase)
            rho = -0.4;
            actual = gofcopula.copulaCDF("t",[0.5 0.5],rho,DF=7);
            expected = 0.25 + asin(rho)/(2*pi);
            testCase.verifyEqual(actual,expected,AbsTol=2e-8);
        end

        function claytonClosedForm(testCase)
            actual = gofcopula.copulaCDF("clayton",[0.5 0.5],2);
            testCase.verifyEqual(actual,1/sqrt(7),AbsTol=2e-14);
        end

        function independenceLimits(testCase)
            U = [0.2 0.7; 0.6 0.4; 0.8 0.9];
            expected = prod(U,2);
            testCase.verifyEqual(gofcopula.copulaCDF("normal",U,0),expected,AbsTol=2e-9);
            testCase.verifyEqual(gofcopula.copulaCDF("clayton",U,0),expected,AbsTol=2e-14);
            testCase.verifyEqual(gofcopula.copulaCDF("frank",U,0),expected,AbsTol=2e-14);
            testCase.verifyEqual(gofcopula.copulaCDF("gumbel",U,1),expected,AbsTol=2e-14);
            testCase.verifyEqual(gofcopula.copulaCDF("joe",U,1),expected,AbsTol=2e-14);
            testCase.verifyEqual(gofcopula.copulaCDF("amh",U,0),expected,AbsTol=2e-14);
            testCase.verifyEqual(gofcopula.copulaCDF("galambos",U,0),expected,AbsTol=2e-14);
            testCase.verifyEqual(gofcopula.copulaCDF("tawn",U,0),expected,AbsTol=2e-14);
            testCase.verifyEqual(gofcopula.copulaCDF("fgm",U,0),expected,AbsTol=2e-14);
            testCase.verifyEqual(gofcopula.copulaCDF("plackett",U,1),expected,AbsTol=2e-14);
        end

        function rotationsFollowDefinitions(testCase)
            U = [0.23 0.67; 0.74 0.42];
            base90 = gofcopula.copulaCDF("clayton",[1-U(:,1),U(:,2)],1.5);
            base180 = gofcopula.copulaCDF("clayton",1-U,1.5);
            base270 = gofcopula.copulaCDF("clayton",[U(:,1),1-U(:,2)],1.5);
            testCase.verifyEqual(gofcopula.copulaCDF("clayton",U,1.5,Rotation=90), ...
                U(:,2)-base90,AbsTol=2e-14);
            testCase.verifyEqual(gofcopula.copulaCDF("clayton",U,1.5,Rotation=180), ...
                sum(U,2)-1+base180,AbsTol=2e-14);
            testCase.verifyEqual(gofcopula.copulaCDF("clayton",U,1.5,Rotation=270), ...
                U(:,1)-base270,AbsTol=2e-14);
        end

        function negativeBivariateClaytonHasCorrectSupport(testCase)
            U = [0.04 0.09; 0.81 0.64];
            expected = max(sqrt(U(:,1))+sqrt(U(:,2))-1,0).^2;
            actual = gofcopula.copulaCDF("clayton",U,-0.5);
            density = gofcopula.copulaPDF("clayton",U,-0.5);
            testCase.verifyEqual(actual,expected,AbsTol=2e-14);
            testCase.verifyEqual(density(1),0);
            testCase.verifyGreaterThan(density(2),0);
        end

        function supportsMultivariateArchimedeanAndElliptical(testCase)
            U = [0.2 0.5 0.8; 0.7 0.4 0.6];
            normal = gofcopula.copulaCDF("normal",U,0.25,Dispersion="exchangeable");
            clayton = gofcopula.copulaCDF("clayton",U,1.2);
            joeDensity = gofcopula.copulaPDF("joe",U,1.5);
            testCase.verifySize(normal,[2 1]);
            testCase.verifySize(clayton,[2 1]);
            testCase.verifyGreaterThan(joeDensity,zeros(2,1));
        end

        function rejectsOutsideUnitCube(testCase)
            testCase.verifyError(@() gofcopula.copulaCDF("normal",[0.2 1.1],0.3), ...
                "gofcopula:copula:OutsideUnitCube");
        end

        function rejectsInvalidParameter(testCase)
            testCase.verifyError(@() gofcopula.copulaPDF("gumbel",[0.2 0.8],0.8), ...
                "gofcopula:copula:InvalidParameter");
        end

        function rejectsExoticHigherDimension(testCase)
            testCase.verifyError(@() gofcopula.copulaCDF("fgm",ones(2,3)/2,0.2), ...
                "gofcopula:copula:UnsupportedDimension");
        end

        function rejectsRotatedHigherDimension(testCase)
            testCase.verifyError(@() gofcopula.copulaCDF( ...
                "clayton",ones(2,3)/2,1.2,Rotation=90), ...
                "gofcopula:copula:RotationDimension");
        end
    end
end
