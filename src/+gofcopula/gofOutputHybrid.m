function results = gofOutputHybrid(result, options)
%GOFOUTPUTHYBRID Select stored hybrid results, keeping single tests.
%   Returns one GofResult per input element containing all single tests
%   plus the selected hybrid rows, with metadata preserved (mirroring R's
%   gofOutputHybrid). Tests selects hybrids that mention any of the given
%   test names; NumberOfTests selects by combination size.
arguments
    result gofcopula.GofResult
    options.Tests string = strings(0,1)
    options.NumberOfTests {mustBeNumeric,mustBeInteger,mustBePositive} = []
end
results = gofcopula.GofResult.empty(0, numel(result));
for k = 1:numel(result)
    singles = result(k).Tests(~startsWith(result(k).Tests.Name,"hybrid"),:);
    hybrid = result(k).Tests(startsWith(result(k).Tests.Name,"hybrid"),:);
    if ~isempty(options.Tests)
        keep = false(height(hybrid),1);
        for j = 1:numel(options.Tests)
            keep = keep | contains(hybrid.Name, options.Tests(j));
        end
        hybrid = hybrid(keep,:);
    end
    if ~isempty(options.NumberOfTests)
        hybrid = hybrid(ismember(count(hybrid.Name,",")+1, options.NumberOfTests),:);
    end
    results(k) = gofcopula.GofResult(Method=result(k).Method, ...
        Copula=result(k).Copula, Margins=result(k).Margins, ...
        MarginParameters=result(k).MarginParameters, ...
        Theta=result(k).Theta, DegreesOfFreedom=result(k).DegreesOfFreedom, ...
        Rotation=result(k).Rotation, Tests=[singles; hybrid], ...
        NumericMode=result(k).NumericMode);
end
end
