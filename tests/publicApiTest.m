classdef publicApiTest < matlab.unittest.TestCase
    %PUBLICAPITEST Contract tests for the namespaced source API.

    methods (TestClassSetup)
        function addToolbox(testCase)
            root=fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(root,"src")));
        end
    end

    methods (Test)
        function capabilityTableHasExpectedShape(testCase)
            actual=gofcopula.CopulaTestTable();
            testCase.verifySize(actual,[16,13]);
        end

        function claytonSupportsKernelInTwoDimensions(testCase)
            actual=gofcopula.gofTest4Copula("clayton",2);
            testCase.verifyTrue(ismember("gofKernel",actual));
        end

        function kernelRejectsThreeDimensions(testCase)
            x=reshape(linspace(0.1,0.9,30),10,3);
            testCase.verifyError(@()gofcopula.gofKernel("clayton",x,M=0, ...
                Param=1,ParamEst=false,Margins="none"), ...
                "gofcopula:UnsupportedCombination");
        end

        function fixedGaussianReturnsResult(testCase)
            x=[0.1 0.2;0.3 0.25;0.5 0.6;0.7 0.8;0.9 0.75];
            actual=gofcopula.gofCvM("normal",x,M=0,Param=0.3, ...
                ParamEst=false,DFEst=false,Margins="none");
            testCase.verifyClass(actual,"gofcopula.GofResult");
            testCase.verifyEqual(actual.Tests.Name,"gofCvM");
            testCase.verifyTrue(isnan(actual.Tests.PValue));
            testCase.verifyGreaterThanOrEqual(actual.Tests.Statistic,0);
        end

        function bootstrapIsRepeatable(testCase)
            x=[0.1 0.2;0.3 0.25;0.5 0.6;0.7 0.8;0.9 0.75];
            first=gofcopula.gofCvM("normal",x,M=4,Param=0.3, ...
                ParamEst=false,DFEst=false,Margins="none",Seed=17);
            second=gofcopula.gofCvM("normal",x,M=4,Param=0.3, ...
                ParamEst=false,DFEst=false,Margins="none",Seed=17);
            testCase.verifyEqual(first.Tests.PValue,second.Tests.PValue,AbsTol=eps);
            testCase.verifyEqual(first.Tests.Statistic,second.Tests.Statistic,AbsTol=1e-14);
        end

        function correctedBootstrapCannotReturnZero(testCase)
            x=[0.1 0.2;0.3 0.25;0.5 0.6;0.7 0.8;0.9 0.75];
            actual=gofcopula.gofCvM("normal",x,M=4,Param=0.3, ...
                ParamEst=false,DFEst=false,Margins="none",Seed=17, ...
                NumericMode="corrected");
            testCase.verifyGreaterThanOrEqual(actual.Tests.PValue,0.2);
        end

        function modelOrchestratorReturnsResult(testCase)
            x=[0.1 0.2;0.3 0.25;0.5 0.6;0.7 0.8;0.9 0.75];
            model=gofcopula.CopulaModel("normal",Theta=0.3, ...
                EstimateTheta=false,EstimateDegreesOfFreedom=false);
            actual=gofcopula.gofco(model,x,Tests="gofCvM",M=0,Margins="none");
            testCase.verifyClass(actual,"gofcopula.GofResult");
        end

        function hybridUtilitiesReturnRequestedRows(testCase)
            tests=table(["first";"second"],[0.02;0.4],[1.2;0.8], ...
                VariableNames=["Name","PValue","Statistic"]);
            source=gofcopula.GofResult(Method="fixture",Copula="normal",Tests=tests);
            generated=gofcopula.gofGetHybrid(source);
            testCase.verifyClass(generated,"gofcopula.GofResult");
            testCase.verifyEqual(generated.Copula,"normal");
            hybridRow=generated.Tests(startsWith(generated.Tests.Name,"hybrid"),:);
            testCase.verifyEqual(hybridRow.PValue,0.04,AbsTol=1e-14);
            selected=gofcopula.gofOutputHybrid(generated,NumberOfTests=2);
            testCase.verifyClass(selected,"gofcopula.GofResult");
            kept=selected.Tests(startsWith(selected.Tests.Name,"hybrid"),:);
            testCase.verifyEqual(kept.PValue,0.04,AbsTol=1e-14);
            named=gofcopula.gofGetHybrid(source,PValues=0.5, ...
                PValueNames="external_study");
            testCase.verifyTrue(any(named.Tests.Name=="external_study"));
        end

        function parallelBootstrapMatchesSerial(testCase)
            testCase.assumeFalse(isempty(ver("parallel")));
            x=[0.1 0.2;0.3 0.25;0.5 0.6;0.7 0.8;0.9 0.75];
            serial=gofcopula.gofCvM("normal",x,M=4,Param=0.3, ...
                ParamEst=false,DFEst=false,Margins="none",Seed=31,Processes=1);
            parallel=gofcopula.gofCvM("normal",x,M=4,Param=0.3, ...
                ParamEst=false,DFEst=false,Margins="none",Seed=31,Processes=2);
            testCase.verifyEqual(parallel.Tests.PValue,serial.Tests.PValue,AbsTol=eps);
            testCase.verifyEqual(parallel.Tests.Statistic,serial.Tests.Statistic,AbsTol=1e-14);
        end
    end
end
