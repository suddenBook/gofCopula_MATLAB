function families = gofCopula4Test(test)
%GOFCOPULA4TEST Return copula families supported by a test.
arguments
    test {mustBeTextScalar}
end
allFamilies = gofcopula.internal.utilities.families();
keep = arrayfun(@(x) gofcopula.internal.utilities.maxDimension(test,x) > 0,allFamilies);
families = allFamilies(keep);
end
