% Build the marginal quantile function Q_g1 from g_1.
%
% Q_g1 is the inverse of F_g1: given u in (0, 1), it returns x such that
% F_g1(x) = u. Used inside MECIP to map pseudo-observations U_hat to
% elliptical-latent Z_hat = Q_g1(U_hat) at each iteration.
%
% Ported from ElliptCopulas::Convert_g1_To_Qg1 (numericalAnalysis.R, lines 99-112).
%
% Inputs:
%   grid - Equally spaced nonnegative grid (grid(1) = 0).
%   g_1  - Values of the 1-D marginal generator on the grid.
%
% Output:
%   Qg1  - Function handle. Qg1(u) for u in (0, 1) returns x on the real line.
%          Linear interpolation; clipped to the interior support to avoid
%          returning NaN at u exactly 0 or 1.

function Qg1 = convertG1ToQg1(grid, g_1)

    grid = grid(:);
    g_1  = g_1(:);
    n1   = numel(grid);
    if n1 ~= numel(g_1)
        error('convertG1ToQg1:SizeMismatch', 'grid and g_1 must be the same length.');
    end

    % Drop endpoints; reflect to build f_1(x) = g_1(x^2) on a symmetric x-grid.
    interior = 2:(n1 - 1);
    g_1_interior  = g_1(interior);
    grid_interior = grid(interior);

    xmax = 2 * max(grid);
    x_nodes = [-xmax; -flipud(sqrt(grid_interior)); sqrt(grid_interior); xmax];
    y_nodes = [0;     flipud(g_1_interior);         g_1_interior;        0];

    [x_nodes_u, ia] = unique(x_nodes, 'first');
    y_nodes_u = y_nodes(ia);
    f_handle = griddedInterpolant(x_nodes_u, y_nodes_u, 'linear', 'nearest');

    % Evaluate f on a uniform-in-x grid (matches R), then cumsum to a CDF.
    new_grid = [-flipud(grid_interior); 0; grid_interior];
    f_vals   = f_handle(new_grid);
    Y        = cumsum(f_vals);
    if Y(end) > 0
        Y = Y / Y(end);
    else
        error('convertG1ToQg1:ZeroCDF', ...
              'Cumulative density is zero; generator has no positive mass on the grid.');
    end

    % Build the inverse map Y -> x. Enforce strict monotonicity on Y.
    [Y_u, x_u] = enforceStrictMonotone(Y, new_grid);

    Qg1_grid = griddedInterpolant(Y_u, x_u, 'linear', 'nearest');
    eps_clip = min(1e-10, Y_u(2) - Y_u(1));
    Qg1 = @(u) Qg1_grid(min(max(u, Y_u(1) + eps_clip), Y_u(end) - eps_clip));
end


function [x_out, y_out] = enforceStrictMonotone(x, y)
% Sort by x, average y over duplicate x, and bump any non-increasing x.
    [x_sorted, order] = sort(x);
    y_sorted = y(order);

    [x_uniq, ~, ic] = unique(x_sorted);
    y_uniq = accumarray(ic, y_sorted, [], @mean);

    dx = diff(x_uniq);
    if any(dx <= 0)
        eps_step = 1e-12 * max(1, max(abs(x_uniq)));
        for k = 2:numel(x_uniq)
            if x_uniq(k) <= x_uniq(k - 1)
                x_uniq(k) = x_uniq(k - 1) + eps_step;
            end
        end
    end

    x_out = x_uniq;
    y_out = y_uniq;
end
