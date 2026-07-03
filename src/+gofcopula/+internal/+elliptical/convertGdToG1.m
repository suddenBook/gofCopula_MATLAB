% Convert a d-dimensional elliptical generator g_d to its 1-D marginal g_1.
%
% For an elliptical distribution with generator g_d, the density of a single
% standardized marginal is f_1(x_1) = g_1(x_1^2), with
%   g_1(t) = (pi^((d-1)/2) / Gamma((d-1)/2)) * int_0^inf u^((d-1)/2 - 1) * g_d(t + u) du.
%
% Ported from ElliptCopulas::Convert_gd_To_g1 (numericalAnalysis.R, lines 56-68).
% Grid-based API: values in, values out, on the same regular grid.
%
% Numerical note: for d = 2 the weight u^(-1/2) has an integrable
% singularity at u = 0. R zeros the first weight, which undercounts the
% first box by ~2*sqrt(step)*g_d(t) (about 15% at step = 0.01 for Gaussian).
% We instead substitute v = sqrt(u) so that
%   int_0^inf u^(-1/2) g_d(t+u) du = 2 * int_0^inf g_d(t + v^2) dv
% and use trapezoidal quadrature on a regular v-grid; this removes the
% singularity and matches the analytical answer to roundoff inside the grid
% support. For d >= 3 the weight is regular at u = 0 and we use the simple
% Riemann form (identical to R up to the (k-1)/k scaling) for code clarity.
%
% Inputs:
%   grid    - Equally spaced nonnegative grid (row or column), grid(1) = 0 assumed.
%   g_d     - Values of g_d on the grid (same size).
%   d       - Dimension (integer >= 2).
%   force_R - Optional struct (default empty). When .singular_I2_rectangle
%             is true and d == 2, reproduces R's flaw: zero the singular
%             boundary weight and use a rectangle sum instead of the
%             v = sqrt(u) substitution. Used only by the R-vs-MATLAB
%             ablation diagnostic (B3 / Audit P6); production keeps it off.
%
% Output:
%   g_1   - Values of g_1 on the same grid. The last point is returned as NaN
%           (no integration tail available), matching the R reference.
%
% Reference: Derumigny & Fermanian (2022), Appendix.

function g_1 = convertGdToG1(grid, g_d, d, force_R)

    if nargin < 4 || isempty(force_R), force_R = struct(); end

    grid = grid(:);
    g_d  = g_d(:);
    n1   = numel(grid);
    if n1 ~= numel(g_d)
        error('convertGdToG1:SizeMismatch', 'grid and g_d must be the same length.');
    end
    if d < 2
        error('convertGdToG1:BadDim', 'd must be >= 2.');
    end

    step = median(diff(grid));
    constant = pi^((d - 1) / 2) / gamma((d - 1) / 2);

    g_1 = nan(n1, 1);

    if d == 2
        if isfield(force_R, 'singular_I2_rectangle') && force_R.singular_I2_rectangle
            % R's flaw: rectangle on u^{-1/2} g_d(t+u) with the singular
            % first weight zeroed (Convert_gd_To_g1, w_puiss[1] = 0).
            w = grid.^(-1/2);
            w(1) = 0;
            for i1 = 1:(n1 - 1)
                k = n1 - i1 + 1;
                g_1(i1) = step * sum(w(1:k) .* g_d(i1:n1));
            end
        else
            % v = sqrt(u) substitution; trapezoidal quadrature over v.
            % g_d evaluated at t + v^2 via linear interpolation.
            u_max   = grid(end);
            n_v     = 2 * n1;               % finer v-grid so t + v^2 hits many grid points
            g_d_fn  = griddedInterpolant(grid, g_d, 'linear', 'nearest');
            for i1 = 1:(n1 - 1)
                t = grid(i1);
                v_max = sqrt(u_max - t);
                v_vec = linspace(0, v_max, n_v)';
                integrand = g_d_fn(t + v_vec.^2);
                g_1(i1) = 2 * trapz(v_vec, integrand);
            end
        end
    else
        % Regular weight for d >= 3 (u^((d-3)/2), finite at u = 0).
        % Use trapezoidal quadrature rather than R's mean*length rectangle
        % rule so range/resolution sensitivity is not polluted by avoidable
        % first-order quadrature bias.
        w_puiss = grid.^((d - 1) / 2 - 1);
        for i1 = 1:(n1 - 1)
            k = n1 - i1 + 1;
            u_grid = grid(1:k);
            g_1(i1) = trapz(u_grid, w_puiss(1:k) .* g_d(i1:n1));
        end
    end

    g_1 = g_1 * constant;
end
