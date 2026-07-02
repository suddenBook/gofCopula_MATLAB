function rows = hybridRows(tests)
%HYBRIDROWS Bonferroni hybrid p-values for every subset of two or more tests.
q = height(tests);
rows = table(Size=[0,3],VariableTypes=["string","double","double"], ...
    VariableNames=["Name","PValue","Statistic"]);
for subsetSize = 2:q
    combinations = nchoosek(1:q,subsetSize);
    for i = 1:size(combinations,1)
        idx = combinations(i,:);
        name = "hybrid(" + join(tests.Name(idx),",") + ")";
        subsetP = tests.PValue(idx);
        if any(isnan(subsetP))
            p = NaN; % a failed test poisons its combinations, as in R
        else
            p = min(1,subsetSize*min(subsetP));
        end
        rows(end+1,:) = {name,p,NaN}; %#ok<AGROW>
    end
end
end
