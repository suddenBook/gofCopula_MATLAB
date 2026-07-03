function maxD = maxDimension(testName, family)
%MAXDIMENSION Largest supported dimension; Inf means any d >= 2.
testName = string(testName); family = lower(string(family));
allFamilies = gofcopula.internal.utilities.families();
if ~ismember(family, allFamilies)
    maxD = 0; return
end
wide = ["normal","t","clayton","gumbel","frank","joe","powerexp"];
bivOnly = ["amh","galambos","huslerreiss","tawn","tev","fgm","plackett"];
if ismember(testName,["gofCvM","gofKS"])
    maxD = inf;
    if ismember(family,bivOnly), maxD = 2; end
elseif ismember(testName,["gofKendallCvM","gofKendallKS","gofCustomTest"])
    maxD = inf;
    if ismember(family,bivOnly), maxD = 2; end
elseif startsWith(testName,"gofRosenblatt")
    maxD = inf;
    if ismember(family,["amh","galambos","fgm","plackett"]), maxD = 2; end
    if ismember(family,["huslerreiss","tawn","tev"]), maxD = 0; end
elseif testName == "gofKernel"
    maxD = 2;
elseif testName == "gofWhite"
    maxD = 2 * ismember(family,wide);
elseif ismember(testName,["gofPIOSTn","gofPIOSRn"])
    if family == "t", maxD = 2;
    elseif ismember(family,["normal","clayton","gumbel","frank","joe"]), maxD = 3;
    elseif family == "powerexp", maxD = inf;
    elseif ismember(family,["amh","galambos","fgm","plackett"]), maxD = 2;
    else, maxD = 0;
    end
elseif startsWith(testName,"gofArchm")
    if ismember(family,["clayton","gumbel","frank","joe"]), maxD = inf;
    elseif family == "amh", maxD = 2;
    else, maxD = 0;
    end
else
    maxD = 0;
end
end
