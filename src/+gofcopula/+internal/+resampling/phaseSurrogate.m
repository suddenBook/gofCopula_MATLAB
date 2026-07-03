function xs = phaseSurrogate(x, stream)
%PHASESURROGATE Coherent phase-randomized surrogate of a time-series matrix.
%   XS = PHASESURROGATE(X, STREAM) returns a surrogate of the [n x d] matrix X
%   that keeps each column's power spectrum (hence its autocorrelation) and the
%   cross-spectrum between columns (hence the sample correlation matrix R), but
%   replaces the phase spectrum with ONE shared random, conjugate-symmetric
%   phase applied to every column.
%
%   The surrogate is a draw from a stationary Gaussian process with X's full
%   second-order structure: its copula is Gaussian and it carries X's serial
%   dependence, while any non-Gaussian copula structure in X is destroyed. It
%   is the Gaussian-copula null used by
%   gofcopula.runTestSerial(..., Method="phase").
%
%   Because the phase is shared across columns, the cross-spectrum
%   X_a(f).*conj(X_b(f)) is preserved exactly, so by Parseval's theorem the
%   sample covariance -- and therefore R -- is preserved to rounding.
%
%   STREAM is a RandStream, so the surrogate is reproducible for a given seed.

n = size(x, 1);
xf = fft(x, [], 1);
phase = zeros(n, 1);
if mod(n, 2) == 0
    m = (n - 2) / 2;
    rp = 2 * pi * rand(stream, m, 1) - pi;
    phase(2:n) = [rp; 0; -flipud(rp)];   % Nyquist bin kept real
else
    m = (n - 1) / 2;
    rp = 2 * pi * rand(stream, m, 1) - pi;
    phase(2:n) = [rp; -flipud(rp)];
end
% Shared phase across columns preserves the cross-spectrum; the DC bin keeps
% its (zero) phase, so the column means are preserved.
xs = real(ifft(abs(xf) .* exp(1i * (angle(xf) + phase)), [], 1));
end
