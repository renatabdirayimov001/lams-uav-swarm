function TR = gen_trace(P, seed)
%GEN_TRACE  Mobility, channel and traffic trace shared by all schemes.
%   The same trace is replayed for every scheme (common random numbers),
%   so differences between schemes come from the forwarding logic only.
set_seed(seed);
N = P.N; T = round(P.T_sim / P.dt); dt = P.dt;

% ---------------- RPGM mobility ----------------
grp = mod(0:N-1, P.G) + 1;                      % group of each UAV
ref = [rand(P.G,1)*P.area(1), rand(P.G,1)*P.area(2), ...
       P.zLim(1) + rand(P.G,1)*diff(P.zLim)];
tgt = new_waypoint(P, P.G);
spd = P.vMin + rand(P.G,1)*(P.vMax - P.vMin);
pauseLeft = zeros(P.G,1);
off = (rand(N,3)*2 - 1) .* [P.rGroup P.rGroup 30] * 0.5;
u = zeros(N,3);                                  % smooth offset velocity

pos = zeros(N,3,T,'single');
for t = 1:T
    for g = 1:P.G
        if pauseLeft(g) > 0
            pauseLeft(g) = pauseLeft(g) - dt; continue;
        end
        d = tgt(g,:) - ref(g,:); L = norm(d); step = spd(g)*dt;
        if L <= step
            ref(g,:) = tgt(g,:);
            tgt(g,:) = new_waypoint(P,1);
            spd(g) = P.vMin + rand*(P.vMax - P.vMin);
            pauseLeft(g) = rand*P.pauseMax;
        else
            ref(g,:) = ref(g,:) + d/L*step;
        end
    end
    % Gauss-Markov offset velocity (smooth intra-group motion) with a weak
    % restoring term that keeps members near their reference point
    u = P.gmAlpha*u + sqrt(1-P.gmAlpha^2)*P.gmSigma*randn(N,3).*[1 1 0.3] ...
        - 0.02*off/dt;
    off = off + u*dt;
    r = sqrt(sum(off(:,1:2).^2,2));
    sc = min(1, P.rGroup ./ max(r,1e-9));
    off(:,1:2) = off(:,1:2) .* sc;
    off(:,3) = max(min(off(:,3), 30), -30);
    pos(:,:,t) = ref(grp,:) + off;
end
vel = zeros(N,3,T,'single');
vel(:,:,2:T) = (pos(:,:,2:T) - pos(:,:,1:T-1)) / dt;
vel(:,:,1) = vel(:,:,2);

% ---------------- air-to-air channel ----------------
N0 = -174 + 10*log10(P.W) + P.NF;
K = 10^(P.Kdb/10);
meanSNR = zeros(N,N,T,'single');
instSNR = zeros(N,N,T,'single');
X = triu(randn(N)*P.sigmaS, 1);
prevRel = [];
for t = 1:T
    p = double(pos(:,:,t));
    D = sqrt(max(sq_dist(p),1));
    % Gudmundson shadowing driven by relative displacement of each pair
    rel = cat(3, p(:,1)-p(:,1)', p(:,2)-p(:,2)', p(:,3)-p(:,3)');
    if ~isempty(prevRel)
        dd = sqrt(sum((rel - prevRel).^2, 3));
        rho = exp(-dd / P.dCorr);
        X = triu(rho.*X + sqrt(1-rho.^2).*P.sigmaS.*randn(N), 1);
    end
    prevRel = rel;
    Xs = X + X';
    m = P.Pt - P.PL0 - 10*P.alpha*log10(D) - Xs - N0;
    m(1:N+1:end) = -inf;
    h = sqrt(K/(K+1)) + sqrt(1/(2*(K+1)))*(randn(N)+1i*randn(N));
    g = triu(abs(h).^2, 1); g = g + g'; g(1:N+1:end) = 1;
    meanSNR(:,:,t) = m;
    instSNR(:,:,t) = m + 10*log10(g);
end

% ---------------- traffic ----------------
nodes = randperm(N);
src = zeros(P.F,1); dst = zeros(P.F,1);
for f = 1:P.F
    s = randi(N); d = randi(N);
    while d == s || grp(d) == grp(s) && rand < 0.75   % mostly inter-group flows
        d = randi(N);
    end
    src(f) = s; dst(f) = d;
end
pk_t = []; pk_f = [];
for f = 1:P.F
    t0 = P.T_warm + rand/P.lambda;
    n = ceil((P.T_genEnd - P.T_warm)*P.lambda*1.5) + 20;
    tt = t0 + cumsum(-log(rand(n,1))/P.lambda);
    tt = tt(tt < P.T_genEnd);
    pk_t = [pk_t; tt]; pk_f = [pk_f; f*ones(numel(tt),1)]; %#ok<AGROW>
end
[pk_t, o] = sort(pk_t); pk_f = pk_f(o);

TR.pos = pos; TR.vel = vel; TR.meanSNR = meanSNR; TR.instSNR = instSNR;
TR.grp = grp; TR.src = src; TR.dst = dst;
TR.pk_t = pk_t; TR.pk_src = src(pk_f); TR.pk_dst = dst(pk_f); TR.pk_flow = pk_f;
TR.N0 = N0; TR.T = T;
TR.Reff = 10^((P.Pt - P.PL0 - N0 - P.gammaTh - P.M)/(10*P.alpha));
end

function w = new_waypoint(P, n)
w = [rand(n,1)*P.area(1), rand(n,1)*P.area(2), P.zLim(1) + rand(n,1)*diff(P.zLim)];
end
