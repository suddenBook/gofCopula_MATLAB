% Build the marginal cumulative distribution function F_g1 from g_1.
%
% The 1-D marginal density of an elliptical distribution with generator g_1
% is f_1(x) = g_1(x^2). This helper builds f_1 on a symmetric x-grid by
% reflecting g_1 around zero (with x = sqrt(xi) on the positive side), then
% cumulatively integrates to obtain the CDF F_g1. A function handle is
% returned that linearly interpolates F_g1 at arbitrary x values.
%
% Ported from ElliptCopulas::Convert_g1_To_Fg1 (numericalAnalysis.R, lines 76-91).
%
% Inputs:
%   grid - Equally spaced nonnegative grid (row or column), grid(1) = 0.
%   g_1  - Values of the 1-D marginal generator on the grid.
%
% Output:
%   Fg1  - Function handle: F_g1(x) for scalar or vector x.
%          Values outside the x-support are extrapolated linearly via
%          the edge pads at +/- 2*max(grid) used inside (matches R).

function Fg1 = convertG1ToFg1(grid, g_1)

    grid = grid(:);
    g_1  = g_1(:);
    n1   = numel(grid);
    if n1 ~= numel(g_1)
        error('convertG1ToFg1:SizeMismatch', 'grid and g_1 must be the same length.');
    end

    % Drop endpoints (endpoints often carry NaN/edge-effects from g_1 construction).
    interior = 2:(n1 - 1);
    g_1_interior  = g_1(interior);
    grid_interior = grid(interior);

    % Symmetric x-grid: f_1(x) = g_1(x^2), so at x = +/- sqrt(xi), f = g_1(xi).
    % Outer pads at +/- 2*max(grid) with value 0 let cumsum taper off cleanly.
    xmax = 2 * max(grid);
    x_nodes = [-xmax; -flipud(sqrt(grid_interior)); sqrt(grid_interior); xmax];
    y_nodes = [0;     flipud(g_1_interior);         g_1_interior;        0];

    % De-duplicate at x = 0 if present (grid_interior(1) could equal 0 in edge cases).
    [x_nodes_u, ia] = uniqueSorted(x_nodes);
    y_nodes_u = accumByKey(ia, y_nodes);

    f_handle = griddedInterpolant(x_nodes_u, y_nodes_u, 'linear', 'nearest');

    % Evaluate on a new x-grid (uniform in x-space) and cumulatively sum.
    % The R reference uses nvlle_grid = [-grid_interior; 0; grid_interior],
    % interpreting grid values as x-coordinates. We match that choice.
    new_grid = [-flipud(grid_interior); 0; grid_interior];
    f_vals   = f_handle(new_grid);
    Y        = cumsum(f_vals);
    if Y(end) > 0
        Y = Y / Y(end);
    end

    % Ensure strict monotonicity for the inverse (used by Qg1). Here for CDF
    % linear interpolation we also want a valid gridded interpolant.
    [new_grid_u, Y_u] = enforceStrictMonotone(new_grid, Y);

    Fg1_grid = griddedInterpolant(new_grid_u, Y_u, 'linear', 'nearest');
    Fg1 = @(x) Fg1_grid(x);
end


function [x_u, idx_first] = uniqueSorted(x)
% Unique entries of a sorted vector, keeping the first occurrence index for
% each unique value.
    [x_sorted, order] = sort(x);
    [x_u, ia] = unique(x_sorted, 'stable');
    idx_first = order(ia);
end


function y_u = accumByKey(first_idx, y)
    y_u = y(first_idx);
end


function [x_out, y_out] = enforceStrictMonotone(x, y)
% Drop duplicate x (keeping the mean of coincident y), then add a tiny
% increment to any y that is not strictly increasing so the result is a
% valid griddedInterpolant target.
    [x_sorted, order] = sort(x);
    y_sorted = y(order);

    % Collapse exact duplicates by averaging.
    [x_uniq, ~, ic] = unique(x_sorted);
    y_uniq = accumarray(ic, y_sorted, [], @mean);

    % Enforce strict monotonicity with a small additive ramp when flat.
    dy = diff(y_uniq);
    if any(dy <= 0)
        eps_step = 1e-12 * max(1, max(abs(y_uniq)));
        for k = 2:numel(y_uniq)
            if y_uniq(k) <= y_uniq(k - 1)
                y_uniq(k) = y_uniq(k - 1) + eps_step;
            end
        end
    end

    x_out = x_uniq;
    y_out = y_uniq;
end
