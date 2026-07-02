function result = CopulaTestTable()
%COPULATESTTABLE Supported copula/test combinations and dimensionalities.
families = gofcopula.internal.utilities.families();
tests = ["gofCvM","gofKS","gofKendallCvM","gofKendallKS", ...
    "gofRosenblattSnB","gofRosenblattSnC","gofRosenblattGamma", ...
    "gofRosenblattChisq","gofKernel","gofWhite","gofPIOSTn", ...
    "gofPIOSRn","gofArchmSnB","gofArchmSnC","gofArchmGamma", ...
    "gofArchmChisq"];
values = strings(numel(tests),numel(families));
for i = 1:numel(tests)
    for j = 1:numel(families)
        maxD = gofcopula.internal.utilities.maxDimension(tests(i),families(j));
        if maxD == 0
            values(i,j) = "-";
        elseif isinf(maxD)
            values(i,j) = ">=2";
        else
            values(i,j) = string(maxD);
        end
    end
end
result = array2table(values,RowNames=cellstr(tests),VariableNames=cellstr(families));
end
