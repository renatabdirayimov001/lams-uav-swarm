function c = fit_calibration(x, y, cDefault)
%FIT_CALIBRATION  Fit the lifetime-error model of Sec. 4.2.
%   x : predicted lifetimes at link birth, y : realised lifetimes.
%   Model:  T = a + b*That + eps,   eps ~ N(0, (s0 + k*That)^2).
%   Returns c = [a b s0 k].
x = x(:); y = y(:);
X = [ones(size(x)) x];
if numel(x) < 5 || rank(X) < 2
    c = cDefault; return;
end
ab = X \ y;
ab(2) = min(max(ab(2), 0), 1.5);          % keep a physically sensible slope
ab(1) = mean(y - ab(2)*x);
res = y - X*ab;
a = abs(res) * sqrt(pi/2);                % E|eps| = sigma*sqrt(2/pi)
w = X \ a;
s0 = min(max(w(1), 0.5), 30);
k  = min(max(w(2), 0), 1);
c = [ab(1) ab(2) s0 k];
if any(~isfinite(c)), c = cDefault; end
end
