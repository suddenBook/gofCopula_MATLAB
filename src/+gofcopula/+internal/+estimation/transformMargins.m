function [u, parameters] = transformMargins(x, margins)
%TRANSFORMMARGINS Convert columns to open-unit-interval pseudo-observations.
arguments
    x {mustBeFloat,mustBeReal,mustBeFinite,mustBeMatrix}
    margins string = "ranks"
end
[n,d] = size(x);
if isscalar(margins), margins = repmat(margins,1,d); end
if numel(margins) ~= d
    error("gofcopula:Margins:WrongCount", ...
        "Margins must be scalar or have one entry per data column.");
end
valid = ["ranks","beta","cauchy","chisq","f","gamma", ...
    "lnorm","norm","t","weibull","exp","none"];
if ~all(ismember(lower(margins),valid))
    error("gofcopula:Margins:Unsupported", "Unsupported margin distribution.");
end
u = zeros(n,d); parameters = cell(1,d);
for j = 1:d
    family = lower(margins(j)); values = x(:,j);
    if family == "ranks"
        % R's ecdf uses the upper rank for ties, not the average tied rank.
        [sorted, order] = sort(values);
        if all(diff(sorted) > 0)
            % No ties: ordinal ranks equal upper ranks. This fast path
            % matters because the corrected-mode bootstrap ranks every
            % replicate sample.
            ranks = zeros(n,1);
            ranks(order) = (1:n).';
            u(:,j) = ranks ./ (n+1);
        else
            u(:,j) = arrayfun(@(z) sum(values <= z),values) ./ (n+1);
        end
        parameters{j} = [];
    elseif family == "none"
        u(:,j) = values; parameters{j} = [];
    else
        [u(:,j),parameters{j}] = fitOneMargin(values,family);
    end
end
if any(u(:) < 0 | u(:) > 1)
    error("gofcopula:Margins:OutsideUnitInterval", ...
        "Data with Margins='none' must lie in [0,1].");
end
% Protect inverse CDFs when callers provide exact boundaries.
u = min(max(u,eps),1-eps);
end

function [u,pars] = fitOneMargin(x,family)
positive = ["beta","chisq","f","gamma","lnorm","weibull"];
if ismember(family,positive) && any(x <= 0)
    error("gofcopula:Margins:PositiveSupport", ...
        "%s margins require positive observations.",family);
end
if family == "exp" && any(x < 0)
    error("gofcopula:Margins:NonnegativeSupport", ...
        "Exponential margins require nonnegative observations.");
end
switch family
    case "norm"
        pd = fitdist(x,"Normal");
    case "beta"
        pd = fitdist(x,"Beta");
    case "gamma"
        pd = fitdist(x,"Gamma");
    case "lnorm"
        pd = fitdist(x,"Lognormal");
    case "weibull"
        pd = fitdist(x,"Weibull");
    case "exp"
        pd = fitdist(x,"Exponential");
    case "t"
        pd = fitdist(x,"tLocationScale");
    case "cauchy"
        objective = @(p) -sum(log(cauchypdf(x,p(1),exp(p(2)))));
        p = fminsearch(objective,[median(x),log(max(iqr(x)/2,eps))]);
        u = cauchycdf(x,p(1),exp(p(2)));
        pars = struct(Family=family,Location=p(1),Scale=exp(p(2)));
        return
    case "chisq"
        df = fminbnd(@(v)-sum(log(chi2pdf(x,v))),eps,1e4);
        u = chi2cdf(x,df); pars = struct(Family=family,DegreesOfFreedom=df);
        return
    case "f"
        p = fminsearch(@(q)-sum(log(max(fpdf(x,exp(q(1)),exp(q(2))),realmin))),[1,1]);
        a=exp(p(1)); b=exp(p(2)); u=fcdf(x,a,b);
        pars=struct(Family=family,DF1=a,DF2=b); return
end
u = cdf(pd,x);
pars = struct(Family=family,Distribution=pd);
end

function y = cauchypdf(x,location,scale)
y = 1 ./ (pi*scale*(1+((x-location)/scale).^2));
end

function y = cauchycdf(x,location,scale)
y = 0.5 + atan((x-location)/scale)/pi;
end
