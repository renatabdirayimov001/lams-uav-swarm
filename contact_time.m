function C = contact_time(dp, dv, R, H)
%CONTACT_TIME  Predicted time until two nodes come within range R.
%   dp, dv : n-by-3 relative positions / velocities (constant-velocity model).
%   Returns 0 if already in range and Inf if no contact within horizon H.
a = sum(dv.^2, 2); b = sum(dp.*dv, 2); c = sum(dp.^2, 2) - R^2;
C = inf(size(a));
C(c <= 0) = 0;
idx = find(c > 0 & a > 1e-6 & b < 0);       % outside and approaching
disc = b(idx).^2 - a(idx).*c(idx);
ok = disc >= 0;
C(idx(ok)) = (-b(idx(ok)) - sqrt(disc(ok))) ./ a(idx(ok));
C(C > H) = inf;
end
