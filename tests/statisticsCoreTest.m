classdef statisticsCoreTest < matlab.unittest.TestCase
    %STATISTICSCORETEST Hand-computable tests for statistic dispatch.

    methods (TestClassSetup)
        function addToolbox(testCase)
            root=fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(root,"src")));
        end
    end

    methods (Test)
        function empiricalCopulaCountsAndOffset(testCase)
            sample=[0.25 0.50; 0.75 0.25];
            points=[0.25 0.50; 1.00 1.00; 0.50 0.25];

            actual=gofcopula.internal.statistics.empiricalCopula(sample,points);
            offset=gofcopula.internal.statistics.empiricalCopula(sample,points,1);

            testCase.verifyEqual(actual,[0.5;1;0],AbsTol=eps);
            testCase.verifyEqual(offset,[1/3;2/3;0],AbsTol=eps);
        end

        function snCMatchesClosedForm(testCase)
            transformed=[0.25 0.50; 0.75 0.25];

            actual=gofcopula.internal.statistics.computeStatistic( ...
                "RosenblattSnC",transformed,[],struct(InputTransformed=true));

            testCase.verifyEqual(actual,0.23828125,AbsTol=1e-14);
        end

        function snBMatchesClosedForm(testCase)
            transformed=[0.25 0.50; 0.75 0.25];

            actual=gofcopula.internal.statistics.computeStatistic( ...
                "SnB",transformed,[],struct(InputTransformed=true));

            testCase.verifyEqual(actual,0.0718315972222222,AbsTol=1e-14);
        end

        function gammaAndersonDarlingMatchesAnchor(testCase)
            transformed=[0.25 0.50; 0.75 0.25];

            actual=gofcopula.internal.statistics.computeStatistic( ...
                "RosenblattGamma",transformed,[],struct(InputTransformed=true));

            testCase.verifyEqual(actual,0.589936609793515,AbsTol=2e-14);
        end

        function kendallKSStepFunctionUsesBothSides(testCase)
            U=[0.25 0.50; 0.75 0.25];
            options=struct(ReferenceCopulaValues=[0.25;0.75]);

            actual=gofcopula.internal.statistics.computeStatistic( ...
                "KendallKS",U,[],options);

            testCase.verifyEqual(actual,sqrt(0.5),AbsTol=1e-14);
        end

        function kendallCvMMatchesRCompatibleFormula(testCase)
            U=[0.25 0.50; 0.75 0.25];
            options=struct(ReferenceCopulaValues=[0.25;0.75]);

            actual=gofcopula.internal.statistics.computeStatistic( ...
                "KendallCvM",U,[],options);

            testCase.verifyEqual(actual,1/6,AbsTol=1e-14);
        end

        function piosRnUsesInformationMatrixRatio(testCase)
            U=[0.2 0.3;0.6 0.8];
            options=struct("ScoreFunction",@(~,~)[1;2], ...
                "HessianFunction",@(~,~)reshape([-2,-3],1,1,2));

            actual=gofcopula.internal.statistics.computeStatistic( ...
                "PIOSRn",U,[],options);

            testCase.verifyEqual(actual,0,AbsTol=eps);
        end

        function inverseGeneratorsHaveKnownValues(testCase)
            clayton=gofcopula.internal.transforms.inverseGenerator( ...
                "clayton",0.5,1);
            gumbel=gofcopula.internal.transforms.inverseGenerator( ...
                "gumbel",0.5,2);

            testCase.verifyEqual(clayton,1,AbsTol=eps);
            testCase.verifyEqual(gumbel,log(2)^2,AbsTol=eps);
        end

        function independenceRosenblattIsIdentity(testCase)
            U=[0.15 0.25;0.35 0.80;0.65 0.45;0.90 0.70];
            model=gofcopula.CopulaModel("normal",Theta=0, ...
                EstimateTheta=false,EstimateDegreesOfFreedom=false);

            transformed=gofcopula.internal.transforms.rosenblattTransform(U,model);

            testCase.verifyEqual(transformed,U,AbsTol=2e-15);
        end

        function archimedeanTransformMatchesClaytonAnchor(testCase)
            U=[0.25 0.50;0.75 0.25];
            model=gofcopula.CopulaModel("clayton",Theta=1, ...
                EstimateTheta=false,EstimateDegreesOfFreedom=false);

            transformed=gofcopula.internal.transforms.archimedeanTransform(U,model);

            testCase.verifyEqual(transformed,[0.75 1/3;0.10 2/3],AbsTol=2e-15);
        end

        function kernelAliasesMatchDeterministicAnchor(testCase)
            U=[0.15 0.25;0.35 0.80;0.65 0.45;0.90 0.70];
            model=gofcopula.CopulaModel("clayton",Theta=1, ...
                EstimateTheta=false,EstimateDegreesOfFreedom=false);
            options=struct("MJ",4,"IntegrationNodes",3,"KernelScale",0.5, ...
                "ModelSample",[0.20 0.30;0.40 0.80;0.60 0.40;0.85 0.75]);

            actual=gofcopula.internal.statistics.computeStatistic( ...
                "Kernel",U,model,options);

            testCase.verifyEqual(actual,0.198983787150003,AbsTol=1e-14);
        end

        function piosTnRefitsEachBlockThroughEstimatorContract(testCase)
            U=[0.15 0.25;0.35 0.80;0.65 0.45;0.90 0.70];
            model=gofcopula.CopulaModel("clayton",Theta=1, ...
                EstimateTheta=false,EstimateDegreesOfFreedom=false);
            options=struct("BlockSize",1,"Estimator",@(~,m)m, ...
                "ParameterCount",1);

            actual=gofcopula.internal.statistics.computeStatistic( ...
                "PIOSTn",U,model,options);

            testCase.verifyEqual(actual,-1,AbsTol=eps);
        end

        function whiteMatchesScalarAnchor(testCase)
            U=[0.2 0.3;0.6 0.8];
            options=struct("ScoreFunction",@(~,~)[1;2], ...
                "HessianFunction",@(~,~)reshape([-0.5,-3],1,1,2), ...
                "DerivativeMatrix",0,"NumericMode","corrected");

            actual=gofcopula.internal.statistics.computeStatistic( ...
                "White",U,[],options);

            testCase.verifyEqual(actual,18,AbsTol=1e-13);
        end

        function whiteIncludesNuisanceParameterCorrection(testCase)
            U=[0.2 0.3;0.6 0.8];
            model=struct("Theta",1);
            options=struct( ...
                "ScoreFunction",@(~,m)[m.Theta;2*m.Theta], ...
                "HessianFunction",@(~,m)reshape( ...
                    [-0.5*m.Theta,-3*m.Theta],1,1,2), ...
                "NumericMode","corrected");

            actual=gofcopula.internal.statistics.computeStatistic( ...
                "White",U,model,options);

            testCase.verifyEqual(actual,98/121,AbsTol=1e-10);
        end

        function customStatisticContract(testCase)
            U=[0.2 0.3;0.6 0.8];

            actual=gofcopula.internal.statistics.computeStatistic( ...
                @(x,~)sum(x,"all"),U,[],struct());

            testCase.verifyEqual(actual,1.9,AbsTol=eps);
        end

        function outsideUnitCubeIsRejected(testCase)
            call=@()gofcopula.internal.statistics.computeStatistic( ...
                "SnC",[0.2 1.1;0.3 0.4],[],struct(InputTransformed=true));

            testCase.verifyError(call,"gofcopula:statistics:OutsideUnitCube");
        end
    end
end
