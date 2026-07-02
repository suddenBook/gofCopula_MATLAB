classdef dataTest < matlab.unittest.TestCase
    %DATATEST Verify the exact structure and reference values of MAT datasets.

    properties (SetAccess = private)
        DataFolder (1,1) string
    end

    methods (TestClassSetup)
        function locateData(testCase)
            projectRoot = fileparts(fileparts(mfilename("fullpath")));
            testCase.DataFolder = fullfile(projectRoot,"data");
        end
    end

    methods (Test)
        function banksPreserveYearsNamesAndValues(testCase)
            actual = load(fullfile(testCase.DataFolder,"Banks.mat"));

            testCase.verifyEqual(actual.Banks.Years,(2004:2012).',AbsTol=eps);
            testCase.verifyEqual(actual.Banks.VariableNames,["C" "BoA"]);
            testCase.verifySize(actual.Banks.Data,[9 1]);
            testCase.verifyEqual(cellfun(@(x) size(x,1),actual.Banks.Data), ...
                [251;251;250;250;252;251;251;251;249],AbsTol=eps);
            testCase.verifyEqual(actual.Banks.Data{1}(1,:), ...
                [1.5913595871077859 0.33736033833593981],AbsTol=1e-15);
            testCase.verifyEqual(actual.Banks.Data{9}(end,:), ...
                [0.680728212704383 1.0159275538131531],AbsTol=1e-15);
            testCase.verifyEqual(actual.Provenance.OriginalPackage,"gofCopula 0.4-3");
        end

        function cryptocurrenciesPreserveYearsNamesAndValues(testCase)
            actual = load(fullfile(testCase.DataFolder,"CryptoCurrencies.mat"));

            testCase.verifyEqual(actual.CryptoCurrencies.Years,(2015:2018).',AbsTol=eps);
            testCase.verifyEqual(actual.CryptoCurrencies.VariableNames, ...
                ["Bitcoin" "Litecoin"]);
            testCase.verifySize(actual.CryptoCurrencies.Data,[4 1]);
            testCase.verifyEqual(cellfun(@(x) size(x,1),actual.CryptoCurrencies.Data), ...
                [364;365;364;364],AbsTol=eps);
            testCase.verifyEqual(actual.CryptoCurrencies.Data{1}(1,:), ...
                [0.10034642701666457 -0.19216390431709432],AbsTol=1e-15);
            testCase.verifyEqual(actual.CryptoCurrencies.Data{4}(end,:), ...
                [-0.80055657994515195 -0.88638304404631141],AbsTol=1e-15);
            testCase.verifyEqual(actual.Provenance.OriginalFile,"CryptoCurrencies.RData");
        end

        function indexReturns2DPreserveMatrixAndNames(testCase)
            actual = load(fullfile(testCase.DataFolder,"IndexReturns2D.mat"));

            testCase.verifySize(actual.IndexReturns2D,[100 2]);
            testCase.verifyEqual(actual.VariableNames,["DAX" "SMI"]);
            testCase.verifyEqual(actual.IndexReturns2D(1,:), ...
                [0.00058993035994880927 0.00079646374327069225],AbsTol=1e-18);
            testCase.verifyEqual(actual.IndexReturns2D(end,:), ...
                [0.021922152290178687 0.016245785397567047],AbsTol=1e-18);
            testCase.verifyEqual(actual.Provenance.OriginalFile,"IndexReturns2D.RData");
        end

        function indexReturns3DPreserveMatrixAndNames(testCase)
            actual = load(fullfile(testCase.DataFolder,"IndexReturns3D.mat"));

            testCase.verifySize(actual.IndexReturns3D,[200 3]);
            testCase.verifyEqual(actual.VariableNames,["DAX" "SMI" "CAC"]);
            testCase.verifyEqual(actual.IndexReturns3D(1,:), ...
                [0.0034901798653415028 0.0038721683866089762 0],AbsTol=1e-18);
            testCase.verifyEqual(actual.IndexReturns3D(end,:), ...
                [0.021922152290178687 0.016245785397567047 0.010897713145171295], ...
                AbsTol=1e-18);
            testCase.verifyEqual(actual.Provenance.OriginalFile,"IndexReturns3D.RData");
        end
    end
end
