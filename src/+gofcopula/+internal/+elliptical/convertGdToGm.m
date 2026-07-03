% Convert an n-dimensional elliptical generator g_n to an m-dimensional
% marginal generator g_m (1 <= m < n).
%
% For an elliptical distribution with generator g_n, every coordinate subset
% of size m is again elliptical with a (generally different-family) marginal
% generator g_m given by
%   g_m(t) = (pi^(p/2) / Gamma(p/2)) * int_0^inf u^(p/2 - 1) * g_n(t + u) du,
% with p = n - m the dimension drop (Hindriks et al. 2025, Section 2.3;
% Cambanis-Huang-Simons 1981).
%
% This is the general-(n -> m) form of convertGdToG1.m, which is the special
% case m = 1 (p = n - 1). The numerical schemes are inherited verbatim from
% convertGdToG1 so that, for m = 1, convertGdToGm(grid, g_n, n, 1) reproduces
% convertGdToG1(grid, g_n, n) bit-for-bit (used as a Phase A self-consistency
% check).
%
% Numerical note: the weight u^(p/2 - 1) has an integrable singularity at
% u = 0 exactly when p = 1 (exponent -1/2). For p = 1 we substitute v = sqrt(u)
% so that
%   int_0^inf u^(-1/2) g_n(t+u) du = 2 * int_0^inf g_n(t + v^2) dv
% and use trapezoidal quadrature on a regular v-grid (this removes the
% singularity and matches the analytical answer to roundoff inside the grid
% support). For p >= 2 the weight is regular at u = 0 and the simple
% trapezoidal Riemann form is used. For the d = 3 branch the only subsets are
% m = 2 (p = 1, singular) and m = 1 (p = 2, regular flat weight u^0).
%
% Inputs:
%   grid    - Equally spaced nonnegative grid (row or column), grid(1) = 0
%             assumed.
%   g_n     - Values of g_n on the grid (same size).
%   n       - Source dimension (integer >= 2).
%   m       - Target marginal dimension (integer, 1 <= m < n).
%   force_R - Optional struct (default empty). When .singular_I2_rectangle
%             is true and p == 1, reproduces R's flaw (zero the singular
%             boundary weight and use a rectangle sum instead of the
%             v = sqrt(u) substitution). Used only by the R-vs-MATLAB
%             ablation diagnostic; production keeps it off.
%
% Output:
%   g_m   - Values of g_m on the same grid. The last point is returned as NaN
%           (no integration tail available), matching convertGdToG1.
%
% Reference: Derumigny & Fermanian (2022), Appendix; Hindriks et al. (2025),
% Section 2.3.

function g_m = convertGdToGm(grid, g_n, n, m, force_R)

    if nargin < 5 || isempty(force_R), force_R = struct(); end

    grid = grid(:);
    g_n  = g_n(:);
    n1   = numel(grid);
    if n1 ~= numel(g_n)
        error('convertGdToGm:SizeMismatch', 'grid and g_n must be the same length.');
    end
    if ~isscalar(n) || n < 2 || mod(n, 1) ~= 0
        error('convertGdToGm:BadSourceDim', 'n must be an integer >= 2.');
    end
    if ~isscalar(m) || m < 1 || m >= n || mod(m, 1) ~= 0
        error('convertGdToGm:BadTargetDim', ...
              'm must be an integer with 1 <= m < n (got m = %g, n = %g).', m, n);
    end

    p = n - m;                              % dimension drop (>= 1)
    step = median(diff(grid));
    constant = pi^(p / 2) / gamma(p / 2);

    g_m = nan(n1, 1);

    if p == 1
        % Exponent p/2 - 1 = -1/2: integrable singularity at u = 0.
        if isfield(force_R, 'singular_I2_rectangle') && force_R.singular_I2_rectangle
            % R's flaw: rectangle on u^{-1/2} g_n(t+u) with the singular
            % first weight zeroed.
            w = grid.^(-1/2);
            w(1) = 0;
            for i1 = 1:(n1 - 1)
                k = n1 - i1 + 1;
                g_m(i1) = step * sum(w(1:k) .* g_n(i1:n1));
            end
        else
            % v = sqrt(u) substitution; trapezoidal quadrature over v.
            % g_n evaluated at t + v^2 via linear interpolation.
            u_max   = grid(end);
            n_v     = 2 * n1;               % finer v-grid so t + v^2 hits many grid points
            g_n_fn  = griddedInterpolant(grid, g_n, 'linear', 'nearest');
            for i1 = 1:(n1 - 1)
                t = grid(i1);
                v_max = sqrt(u_max - t);
                v_vec = linspace(0, v_max, n_v)';
                integrand = g_n_fn(t + v_vec.^2);
                g_m(i1) = 2 * trapz(v_vec, integrand);
            end
        end
    else
        % Regular weight for p >= 2 (u^(p/2 - 1), finite at u = 0).
        % Trapezoidal quadrature rather than a mean*length rectangle rule so
        % range/resolution sensitivity is not polluted by avoidable
        % first-order quadrature bias.
        w_puiss = grid.^(p / 2 - 1);
        for i1 = 1:(n1 - 1)
            k = n1 - i1 + 1;
            u_grid = grid(1:k);
            g_m(i1) = trapz(u_grid, w_puiss(1:k) .* g_n(i1:n1));
        end
    end

    g_m = g_m * constant;
end
