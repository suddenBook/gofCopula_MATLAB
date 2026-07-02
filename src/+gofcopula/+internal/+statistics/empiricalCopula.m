function values = empiricalCopula(sample, points, offset)
%EMPIRICALCOPULA Evaluate the empirical distribution at arbitrary points.
%   VALUES = EMPIRICALCOPULA(SAMPLE, POINTS) returns
%      n^(-1) sum_i I(SAMPLE(i,:) <= POINTS(j,:)).
%   OFFSET changes the denominator to size(SAMPLE,1)+OFFSET. This matches
%   copula::F.n(..., offset=OFFSET); in particular OFFSET=1 is used by the
%   Hering--Hofert Archimedean transform.

arguments
    sample (:,:) double {mustBeReal,mustBeFinite}
    points (:,:) double {mustBeReal,mustBeFinite}
    offset (1,1) double {mustBeReal,mustBeFinite,mustBeNonnegative} = 0
end

if size(sample,2) ~= size(points,2)
    error("gofcopula:statistics:DimensionMismatch", ...
        "Sample and evaluation points must have the same number of columns.");
end
if isempty(sample)
    error("gofcopula:statistics:EmptySample", "Sample must not be empty.");
end

n = size(sample,1);
values = zeros(size(points,1),1);
% Evaluate blocks by implicit expansion. The fixed cap keeps the temporary
% n-by-d-by-block logical array modest even for the 10,000-point Kendall
% reference sample, while avoiding one interpreted loop per evaluation point.
blockSize = 256;
for first = 1:blockSize:size(points,1)
    last = min(first+blockSize-1,size(points,1));
    block = reshape(points(first:last,:).',[1,size(points,2),last-first+1]);
    counts = sum(all(sample <= block,2),1);
    values(first:last) = reshape(counts,[],1) / (n + offset);
end
end
