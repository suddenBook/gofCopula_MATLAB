function model = adjustDegreesOfFreedom(testName, model, mode)
%ADJUSTDEGREESOFFREEDOM R-compatible Student df rules of gofCopula 0.4-3.
%   In NumericMode "rCompatible" the R package modifies the fitted t/tev
%   degrees of freedom in place (original R internal_bootstrap.R lines 51-62 and
%   173-183): the CvM and KS statistics are evaluated with df ceiled to an
%   integer, and the PIOS Tn statistic caps df at 60. The modified model is
%   also the one replicates are simulated from and refits start at, exactly
%   as in R. "corrected" mode leaves the fractional df untouched.

if mode ~= "rCompatible" || ~isa(model, "gofcopula.CopulaModel")
    return
end
if ~ismember(model.Family, ["t", "tev"])
    return
end
if ismember(testName, ["gofCvM", "gofKS"])
    model = model.withFit(model.Theta, ceil(model.DegreesOfFreedom));
elseif testName == "gofPIOSTn"
    model = model.withFit(model.Theta, min(model.DegreesOfFreedom, 60));
end
end
