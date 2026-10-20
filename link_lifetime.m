function T = link_lifetime(dp, dv, R, H)
%LINK_LIFETIME  Predicted residual lifetime of existing links (eq. 5).
%   dp, dv : n-by-3 relative positions / velocities, R : effective range,
%   H : horizon. Returns the time until ||dp + dv*tau|| = R, capped at H.
a = sum(dv.^2, 2); b = sum(dp.*dv, 2); c = sum(dp.^2, 2) - R^2;
T = H * ones(size(a));
m = a > 1e-6;
disc = b(m).^2 - a(m).*c(m);
tau = (-b(m) + sqrt(max(disc, 0))) ./ a(m);
T(m) = tau;
T(c > 0 & b >= 0) = 0;          % already outside and separating
T = min(max(T, 0), H);
end
