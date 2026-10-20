function RT = route_trees(A, TL, RK, d, beta)
%ROUTE_TREES  Destination-rooted routing trees used by every scheme.
%   A  : known adjacency (logical N-by-N), TL : predicted link lifetimes,
%   RK : link rates (bit/s), d : destination, beta : weight of 1/T.
%   Returns for each tree m = 1 (min-hop), 2 (max-lifetime), 3 (lifetime-
%   weighted) and each node i: next hop nh(m,i), hops h(m,i), bottleneck
%   lifetime life(m,i) and bottleneck rate rate(m,i). nh = 0: no route.
N = size(A,1);
RT.nh = zeros(3,N); RT.h = inf(3,N); RT.life = zeros(3,N); RT.rate = zeros(3,N);

% ---- tree 1: minimum hop (BFS) ----
nh = zeros(1,N); hop = inf(1,N); hop(d) = 0; order = d; head = 1;
while head <= numel(order)
    v = order(head); head = head + 1;
    nb = find(A(:,v)' & isinf(hop));
    hop(nb) = hop(v) + 1; nh(nb) = v; order = [order nb]; %#ok<AGROW>
end
RT = fill_tree(RT, 1, nh, order, TL, RK, d);

% ---- tree 2: maximum bottleneck lifetime (widest path) ----
b = -inf(1,N); b(d) = inf; hb = inf(1,N); hb(d) = 0;
nh = zeros(1,N); done = false(1,N); order = [];
for it = 1:N
    cand = b; cand(done) = -inf;
    [bv, v] = max(cand - 1e-6*hb);  % tie-break on fewer hops
    if isinf(b(v)) && b(v) < 0, break; end
    if bv == -inf, break; end
    done(v) = true; order = [order v]; %#ok<AGROW>
    nb = find(A(:,v)' & ~done);
    c = min(b(v), TL(nb,v)');
    upd = c > b(nb) | (c == b(nb) & hb(v)+1 < hb(nb));
    b(nb(upd)) = c(upd); hb(nb(upd)) = hb(v) + 1; nh(nb(upd)) = v;
end
RT = fill_tree(RT, 2, nh, order, TL, RK, d);

% ---- tree 3: lifetime-weighted shortest path (Dijkstra) ----
W = 1 + beta ./ max(TL, 0.1);
dist = inf(1,N); dist(d) = 0; nh = zeros(1,N); done = false(1,N); order = [];
for it = 1:N
    cand = dist; cand(done) = inf;
    [dv, v] = min(cand);
    if isinf(dv), break; end
    done(v) = true; order = [order v]; %#ok<AGROW>
    nb = find(A(:,v)' & ~done);
    c = dist(v) + W(nb,v)';
    upd = c < dist(nb);
    dist(nb(upd)) = c(upd); nh(nb(upd)) = v;
end
RT = fill_tree(RT, 3, nh, order, TL, RK, d);
end

function RT = fill_tree(RT, m, nh, order, TL, RK, d)
% Propagate hop count, bottleneck lifetime and bottleneck rate from d outwards.
N = numel(nh);
h = inf(1,N); life = zeros(1,N); rate = zeros(1,N);
h(d) = 0; life(d) = inf; rate(d) = inf;
for v = order
    if v == d, continue; end
    u = nh(v);
    if u == 0 || isinf(h(u)), continue; end
    h(v) = h(u) + 1;
    life(v) = min(TL(v,u), life(u));
    rate(v) = min(RK(v,u), rate(u));
end
nh(d) = 0;
RT.nh(m,:) = nh .* isfinite(h); RT.h(m,:) = h;
RT.life(m,:) = life; RT.rate(m,:) = rate;
end
