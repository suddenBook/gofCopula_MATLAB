function u = replicatePipeline(simulated, margins, marginParameters)
%REPLICATEPIPELINE Apply the observed-data margins pipeline to a replicate.
%   Implements the Genest--Remillard (2008) parametric bootstrap for
%   estimated margins: each simulated copula sample is passed through the
%   SAME pipeline the observed data went through, so replicate statistics
%   share the observed statistic's null distribution.
%     "ranks" -> pseudo-observations of the simulated sample;
%     "none"  -> identity (the observed pipeline was the identity);
%     parametric margins -> map to the data scale with the margins fitted
%       on the observed sample, then refit the margins on the replicate.
%   Used only in NumericMode "corrected"; "rCompatible" reproduces R, which
%   computes replicate statistics on the raw copula draws.

d = size(simulated, 2);
margins = string(margins);
if isscalar(margins)
    margins = repmat(margins, 1, d);
end
if all(margins == "none")
    u = simulated;
    return
end

x = simulated;
for j = 1:d
    family = lower(margins(j));
    if family == "ranks" || family == "none"
        continue % transformMargins ranks/passes these columns below
    end
    if numel(marginParameters) < j || isempty(marginParameters{j})
        error("gofcopula:Bootstrap:MissingMarginParameters", ...
            "Fitted margin parameters for column %d are required.", j);
    end
    x(:,j) = gofcopula.internal.estimation.inverseMargins( ...
        simulated(:,j), margins(j), marginParameters{j});
end
u = gofcopula.internal.estimation.transformMargins(x, margins);
end
