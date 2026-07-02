function results = gofGetHybrid(result, options)
%GOFGETHYBRID Compute hybrid p-values, optionally with external p-values.
%   Returns one GofResult per input element carrying the single tests, any
%   external p-values (named via PValueNames), and the requested hybrid
%   combinations, with the input's metadata preserved (mirroring R's
%   gofGetHybrid, which returns a full gofCOP object).
arguments
    result gofcopula.GofResult
    options.PValues {mustBeNumeric,mustBeReal,mustBeFinite} = []
    options.PValueNames string = strings(0,1)
    options.NumberOfTests {mustBeNumeric,mustBeInteger,mustBePositive} = []
end
results = gofcopula.GofResult.empty(0, numel(result));
for k = 1:numel(result)
    base = result(k).Tests(~startsWith(result(k).Tests.Name,"hybrid"),:);
    for j = 1:numel(options.PValues)
        if numel(options.PValueNames) >= j
            name = options.PValueNames(j);
        else
            name = "external" + j;
        end
        base(end+1,:) = {name, options.PValues(j), NaN}; %#ok<AGROW>
    end
    hybrid = gofcopula.internal.utilities.hybridRows(base);
    if ~isempty(options.NumberOfTests)
        counts = count(hybrid.Name, ",") + 1;
        hybrid = hybrid(ismember(counts, options.NumberOfTests), :);
    end
    results(k) = gofcopula.GofResult(Method="Hybrid combination of tests", ...
        Copula=result(k).Copula, Margins=result(k).Margins, ...
        MarginParameters=result(k).MarginParameters, ...
        Theta=result(k).Theta, DegreesOfFreedom=result(k).DegreesOfFreedom, ...
        Rotation=result(k).Rotation, Tests=[base; hybrid], ...
        NumericMode=result(k).NumericMode);
end
end
