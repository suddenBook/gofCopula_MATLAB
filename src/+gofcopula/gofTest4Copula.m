function tests = gofTest4Copula(copula, d)
%GOFTEST4COPULA Return tests supporting a family and dimension.
arguments
    copula {mustBeTextScalar} = "clayton"
    d (1,1) {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger,mustBeGreaterThanOrEqual(d,2)} = 2
end
family = lower(string(copula));
if family == "gaussian", family = "normal"; end
mustBeMember(family,gofcopula.internal.utilities.families());
tbl = gofcopula.CopulaTestTable();
names = string(tbl.Properties.RowNames);
keep = arrayfun(@(x) gofcopula.internal.utilities.maxDimension(x,family) >= d,names);
tests = names(keep);
end
