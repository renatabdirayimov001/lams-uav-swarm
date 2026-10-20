function R = run_sim(P, TR, scheme, seed)
%RUN_SIM  Replay trace TR with one forwarding scheme and return metrics.
%   scheme: 'LAMS', 'LAMS-noAlt', 'LAMS-noCal', 'TH', 'RC', 'MH', 'EPI', 'SNW'
%   seed  : seed for scheme-internal randomness (processing order).
set_seed(seed + 7919);
N = P.N; T = TR.T; dt = P.dt; S = P.S; H = P.H; Reff = TR.Reff;
isLams  = strncmp(scheme, 'LAMS', 4);
useAlt  = ~strcmp(scheme, 'LAMS-noAlt');
useCal  = ~strcmp(scheme, 'LAMS-noCal');
multi   = any(strcmp(scheme, {'EPI','SNW'}));
updEvery = max(1, round(P.T_upd/dt));

% ---------------- packets ----------------
Np = numel(TR.pk_t);
pt = TR.pk_t; pdst = TR.pk_dst; psrc = TR.pk_src;
state  = zeros(Np,1);           % 0 not yet, 1 active, 2 delivered, 3 dropped
holder = zeros(Np,1);           % single-copy schemes
arr    = zeros(Np,1);           % arrival time at current holder (FIFO key)
tdel   = nan(Np,1);
hops   = zeros(Np,1);
lat    = zeros(Np,1); latSlot = zeros(Np,1);
noRoute = nan(Np,1);
usedSCF = false(Np,1); usedAlt = false(Np,1);
dropReason = zeros(Np,1);       % 1 buffer at source, 2 TTL, 3 no route (MH)
if multi
    Hc = false(N, Np);          % copy matrix
    Lt = zeros(N, Np);          % Spray-and-Wait tokens
end
occ = zeros(N,1);
nextPk = 1;

% ---------------- knowledge ----------------
KP = zeros(N,N,3); KV = zeros(N,N,3); KT = -inf(N,N);   % i's view of j
A = false(N); Aprev = false(N); LQ = []; TL = zeros(N); RK = zeros(N); tUpd = 0;
comp = zeros(1,N);
trees = cell(1,N);              % trees{d}
born = zeros(N); pred0 = zeros(N);
calX = zeros(N, P.Wc); calY = zeros(N, P.Wc); calN = zeros(N,1);
cal = repmat(P.calDefault, N, 1);
allPred = []; allObs = [];
modeState = zeros(N,N);         % current mode of (node, destination)
Dall = inf(N,N); scfStamp = zeros(1,N);   % cached SCF delays Dall(d,x)
nSwitch = 0;

% ---------------- counters ----------------
txData = 0; ctrl = 0; occSum = 0; occCnt = 0;
txMode = zeros(1,4);            % transmissions in M1, M2, SCF custody, other

for t = 1:T
    tnow = (t-1)*dt;
    pos = double(TR.pos(:,:,t)); vel = double(TR.vel(:,:,t));
    mS = TR.meanSNR(:,:,t); iS = TR.instSNR(:,:,t);

    % ---- inject new packets at their sources ----
    while nextPk <= Np && pt(nextPk) < tnow
        p = nextPk; s = psrc(p);
        if occ(s) < P.B
            state(p) = 1; occ(s) = occ(s) + 1; arr(p) = pt(p);
            if multi, Hc(s,p) = true; Lt(s,p) = P.snwL; else, holder(p) = s; end
        else
            state(p) = 3; dropReason(p) = 1;
        end
        nextPk = nextPk + 1;
    end

    % ---- link-state update (HELLO + TC) ----
    if mod(t-1, updEvery) == 0
        Aprev = A;
        % link quality: EWMA of mean SNR over HELLO periods (OLSR-like hysteresis)
        if isempty(LQ), LQ = double(mS); else, LQ = P.lqAlpha*LQ + (1-P.lqAlpha)*double(mS); end
        A = LQ >= P.gammaTh; A(1:N+1:end) = false;
        comp = components(A);
        for c = unique(comp)
            mem = find(comp == c);
            for i = mem
                KP(i,mem,:) = reshape(pos(mem,:), [1 numel(mem) 3]);
                KV(i,mem,:) = reshape(vel(mem,:), [1 numel(mem) 3]);
                KT(i,mem) = tnow;
            end
            if multi
                ctrl = ctrl + numel(mem);              % HELLO + summary vectors
            else
                ctrl = ctrl + numel(mem) + numel(mem)^2;  % HELLO + TC flooding
            end
        end
        [ii, jj] = find(triu(A,1));
        TL = zeros(N);
        if ~isempty(ii)
            lt = link_lifetime(pos(jj,:)-pos(ii,:), vel(jj,:)-vel(ii,:), Reff, H);
            TL(sub2ind([N N],ii,jj)) = lt; TL = TL + TL';
        end
        RK = rate_of(mS, P) .* A;
        % calibration samples: predicted (at birth) vs realised lifetime
        nb = triu(A & ~Aprev, 1); [bi,bj] = find(nb);
        born(sub2ind([N N],bi,bj)) = tnow; pred0(sub2ind([N N],bi,bj)) = TL(sub2ind([N N],bi,bj));
        en = triu(Aprev & ~A, 1); [ei,ej] = find(en);
        for k = 1:numel(ei)
            a_ = ei(k); b_ = ej(k);
            if born(a_,b_) <= 0, continue; end        % born before first update
            x = pred0(a_,b_); y = tnow - born(a_,b_);
            allPred(end+1) = x; allObs(end+1) = y; %#ok<AGROW>
            for v = [a_ b_]
                calN(v) = calN(v) + 1; s_ = mod(calN(v)-1, P.Wc) + 1;
                calX(v,s_) = x; calY(v,s_) = y;
            end
        end
        if useCal
            for v = find(calN' >= 10)
                n = min(calN(v), P.Wc);
                cal(v,:) = fit_calibration(calX(v,1:n), calY(v,1:n), P.calDefault);
            end
        end
        % routing trees for destinations with active packets
        trees = cell(1,N);
        if ~multi
            for d = unique(pdst(state == 1))'
                trees{d} = route_trees(A, TL, RK, d, P.beta);
            end
        end
        tUpd = tnow;
    end
    age = tnow - tUpd;

    % ---- TTL expiry ----
    act = find(state == 1);
    exp_ = act(tnow - pt(act) > P.TTL);
    if ~isempty(exp_)
        state(exp_) = 3; dropReason(exp_) = 2;
        if multi
            occ = occ - sum(Hc(:,exp_), 2); Hc(:,exp_) = false; Lt(:,exp_) = 0;
        else
            for p = exp_', occ(holder(p)) = occ(holder(p)) - 1; holder(p) = 0; end
        end
    end

    air = dt*ones(N,1);
    failed = false(N);

    if ~multi
        % ================= single-copy schemes =================
        rank_ = randperm(N);
        for rnd = 1:P.maxRounds
            act = find(state == 1);
            if isempty(act), break; end
            M_ = sortrows([rank_(holder(act))' pdst(act) arr(act) act]);
            moved = 0;
            starts = [1; find(any(diff(M_(:,1:2),1,1),2))+1; size(M_,1)+1];
            for sg = 1:numel(starts)-1
                q = M_(starts(sg):starts(sg+1)-1, 4);
                i = holder(q(1)); d = pdst(q(1)); RT = trees{d};
                if isempty(RT)
                    RT.nh = zeros(3,N); RT.h = inf(3,N); RT.life = zeros(3,N); RT.rate = zeros(3,N);
                end
                inComp = RT.nh(1,i) > 0;
                DscfI = []; uBest = -1;
                kk = 1; guard = 0;
                while kk <= numel(q) && guard < 3*numel(q) + 10
                    guard = guard + 1;
                    p = q(kk);
                    if state(p) ~= 1 || holder(p) ~= i, kk = kk + 1; continue; end
                    action = 3; j = 0; mode = 4;   % 1 tx, 2 stop, 3 next, 4 drop
                    switch scheme
                        case {'LAMS','LAMS-noAlt','LAMS-noCal'}
                            if useAlt, ms = 1:3; else, ms = 1; end
                            Jbest = inf; mBest = 0;
                            c_ = cal(i,:);
                            for m = ms
                                nhm = RT.nh(m,i);
                                if nhm == 0 || failed(i,nhm), continue; end
                                Tm = RT.life(m,i) - age; rm = min(RT.rate(m,i), P.Rmax);
                                tau = kk*S/rm + RT.h(m,i)*(S/rm + P.dProc);
                                sig = c_(3) + c_(4)*max(Tm,0);
                                Pk = 0.5*erfc(-((c_(1) + c_(2)*Tm - tau)/sig)/sqrt(2));
                                J = tau + P.lam*(1 - Pk);
                                if J < Jbest, Jbest = J; mBest = m; end
                            end
                            if isempty(DscfI)
                                [Dall, scfStamp(d)] = get_scf(Dall, scfStamp(d), t, d, comp, pos, vel, KP, KV, KT, tnow, Reff, P, RT);
                                DscfI = Dall(d,i);
                            end
                            ttlRem = P.TTL - (tnow - pt(p));
                            Pscf = double(DscfI <= min(H, ttlRem)) + P.p0*double(DscfI > min(H, ttlRem));
                            Jscf = DscfI + P.lam*(1 - Pscf);
                            if isinf(DscfI), Jscf = H + P.lam*(1 - P.p0); end
                            if Jscf < Jbest, mSel = 4; Jsel = Jscf; else, mSel = mBest; Jsel = Jbest; end
                            % hysteresis against oscillation
                            cur = modeState(i,d);
                            if cur ~= 0 && cur ~= mSel
                                if cur == 4, Jcur = Jscf;
                                else
                                    Jcur = inf; nhm = RT.nh(cur,i);
                                    if nhm > 0 && ~failed(i,nhm), Jcur = mode_cost(RT, cur, i, kk, age, cal(i,:), P); end
                                end
                                if ~(Jsel < (1 - P.eta)*Jcur), mSel = cur; end
                            end
                            if mSel ~= cur
                                if cur ~= 0, nSwitch = nSwitch + 1; end
                                modeState(i,d) = mSel;
                            end
                            if mSel >= 1 && mSel <= 3 && RT.nh(mSel,i) > 0 && ~failed(i,RT.nh(mSel,i))
                                action = 1; j = RT.nh(mSel,i); mode = 1 + (mSel > 1);
                            else
                                usedSCF(q(kk:end)) = true;
                                if uBest < 0
                                    uBest = best_carrier(i, DscfI, Dall(d,:), A, failed, occ, P);
                                end
                                if uBest > 0, action = 1; j = uBest; mode = 3; else, action = 2; end
                            end
                        case 'TH'
                            nhm = RT.nh(1,i);
                            if nhm > 0 && ~failed(i,nhm) && RT.life(1,i) - age >= P.thTmin
                                action = 1; j = nhm; mode = 1;
                            else
                                usedSCF(q(kk:end)) = true;
                                if isempty(DscfI)
                                    [Dall, scfStamp(d)] = get_scf(Dall, scfStamp(d), t, d, comp, pos, vel, KP, KV, KT, tnow, Reff, P, RT);
                                    DscfI = Dall(d,i);
                                end
                                if uBest < 0
                                    uBest = best_carrier(i, DscfI, Dall(d,:), A, failed, occ, P);
                                end
                                if uBest > 0, action = 1; j = uBest; mode = 3; else, action = 2; end
                            end
                        case 'RC'
                            nhm = RT.nh(3,i);
                            if nhm > 0 && ~failed(i,nhm)
                                action = 1; j = nhm; mode = 1;
                            else
                                g = greedy_next(i, d, A, failed, occ, pos, KP, KV, KT, tnow, P);
                                if g > 0, action = 1; j = g; mode = 4;
                                else, action = 2; usedSCF(q(kk:end)) = true; end
                            end
                        case 'MH'
                            nhm = RT.nh(1,i);
                            if nhm > 0
                                noRoute(p) = nan;
                                if ~failed(i,nhm), action = 1; j = nhm; mode = 1; else, action = 2; end
                            else
                                if isnan(noRoute(p)), noRoute(p) = tnow; end
                                if tnow - noRoute(p) > P.T_hold, action = 4; else, action = 3; end
                            end
                    end
                    if action == 2, break; end
                    if action == 3, kk = kk + 1; continue; end
                    if action == 4
                        state(p) = 3; dropReason(p) = 3; occ(i) = occ(i) - 1; holder(p) = 0;
                        kk = kk + 1; continue;
                    end
                    % ---- transmission attempt i -> j ----
                    r_ = rate_of(max(iS(i,j), mS(i,j)), P);
                    s_ = S / r_;
                    if air(i) < s_ || air(j) < s_ || (j ~= d && occ(j) >= P.B), break; end
                    air(i) = air(i) - s_; air(j) = air(j) - s_; txData = txData + 1;
                    if iS(i,j) >= P.gammaTh
                        txMode(mode) = txMode(mode) + 1;
                        if mode == 2, usedAlt(p) = true; end
                        if latSlot(p) ~= t, lat(p) = 0; latSlot(p) = t; end
                        lat(p) = lat(p) + s_ + P.dProc;
                        hops(p) = hops(p) + 1; occ(i) = occ(i) - 1; moved = moved + 1;
                        if j == d
                            state(p) = 2; tdel(p) = tnow + lat(p); holder(p) = 0;
                        else
                            holder(p) = j; occ(j) = occ(j) + 1; arr(p) = tnow + lat(p);
                        end
                        kk = kk + 1;
                    else
                        failed(i,j) = true; failed(j,i) = true;
                        uBest = -1;
                    end
                end
            end
            if moved == 0, break; end
        end
    else
        % ================= multi-copy schemes (EPI, SNW) =================
        [ii, jj] = find(triu(A,1));
        o = randperm(numel(ii)); ii = ii(o); jj = jj(o);
        act = find(state == 1)';
        for e = 1:numel(ii)
            for dirn = 1:2
                if dirn == 1, i = ii(e); j = jj(e); else, i = jj(e); j = ii(e); end
                if failed(i,j), continue; end
                cand = act(Hc(i,act) & ~Hc(j,act) & state(act)' == 1);
                if isempty(cand), continue; end
                if strcmp(scheme,'SNW')
                    cand = cand(Lt(i,cand) > 1 | pdst(cand)' == j);
                end
                if isempty(cand), continue; end
                toDst = pdst(cand)' == j;
                cand = [cand(toDst) cand(~toDst)];
                for p = cand
                    if state(p) ~= 1, continue; end
                    r_ = rate_of(max(iS(i,j), mS(i,j)), P); s_ = S / r_;
                    if air(i) < s_ || air(j) < s_, break; end
                    if pdst(p) ~= j && occ(j) >= P.B, break; end
                    air(i) = air(i) - s_; air(j) = air(j) - s_; txData = txData + 1;
                    if iS(i,j) < P.gammaTh, failed(i,j) = true; failed(j,i) = true; break; end
                    txMode(4) = txMode(4) + 1;
                    if latSlot(p) ~= t, lat(p) = 0; latSlot(p) = t; end
                    if pdst(p) == j
                        state(p) = 2; tdel(p) = tnow + lat(p) + s_ + P.dProc;
                        occ = occ - Hc(:,p); Hc(:,p) = false; Lt(:,p) = 0;   % carriers learn of delivery lazily; freed here for simplicity of accounting
                    else
                        Hc(j,p) = true; occ(j) = occ(j) + 1;
                        if strcmp(scheme,'SNW')
                            give = floor(Lt(i,p)/2); Lt(j,p) = give; Lt(i,p) = Lt(i,p) - give;
                        end
                    end
                end
            end
        end
    end

    % ---- statistics ----
    if tnow >= P.T_warm
        occSum = occSum + mean(occ); occCnt = occCnt + 1;
    end
end

% ---------------- metrics ----------------
gen = state > 0 | pt < P.T_sim;           % all generated packets
nGen = sum(gen);
dl = state == 2; nDel = sum(dl);
dly = tdel(dl) - pt(dl);
R.scheme = scheme;
R.nGen = nGen; R.nDel = nDel;
R.pdr = nDel / max(nGen,1);
R.delayMean = mean(dly);
R.delayP95 = prctile_(dly, 95);
R.txPerDel = txData / max(nDel,1);
R.ctrlPerDel = ctrl / max(nDel,1);
R.occMean = occSum / max(occCnt,1);
R.dropBuf = sum(dropReason == 1)/nGen;
R.dropTTL = sum(dropReason == 2)/nGen;
R.dropNoRoute = sum(dropReason == 3)/nGen;
R.undelivered = sum(state == 1)/nGen;
R.shareSCF = sum(dl & usedSCF)/max(nDel,1);
R.shareAlt = sum(dl & usedAlt)/max(nDel,1);
R.switchPerMin = nSwitch / ((P.T_sim - P.T_warm)/60) / N;
R.txMode = txMode;
R.hopsMean = mean(hops(dl));
R.calPred = allPred; R.calObs = allObs;
R.calFinal = cal;
end

% =====================================================================
function c = components(A)
N = size(A,1); c = zeros(1,N); k = 0;
for s = 1:N
    if c(s), continue; end
    k = k + 1; c(s) = k; st = s;
    while ~isempty(st)
        v = st(end); st(end) = [];
        nb = find(A(v,:) & c == 0); c(nb) = k; st = [st nb]; %#ok<AGROW>
    end
end
end

function r = rate_of(snr, P)
r = min(P.Rmax, P.W*log2(1 + 10.^((snr - P.gapdB)/10)));
r = max(r, 1e3);
end

function J = mode_cost(RT, m, i, kk, age, c_, P)
Tm = RT.life(m,i) - age; rm = min(RT.rate(m,i), P.Rmax);
tau = kk*P.S/rm + RT.h(m,i)*(P.S/rm + P.dProc);
sig = c_(3) + c_(4)*max(Tm,0);
Pk = 0.5*erfc(-((c_(1) + c_(2)*Tm - tau)/sig)/sqrt(2));
J = tau + P.lam*(1 - Pk);
end

function [p, v] = view_of(i, x, pos, vel, KP, KV, KT, tnow)
% Node i's estimate of node x's current position and velocity.
if x == i
    p = pos(i,:); v = vel(i,:);
elseif isinf(KT(i,x))
    p = nan(1,3); v = zeros(1,3);
else
    v = reshape(KV(i,x,:),1,3);
    p = reshape(KP(i,x,:),1,3) + v*(tnow - KT(i,x));
end
end

function [Dall, stamp] = get_scf(Dall, stamp, t, d, comp, pos, vel, KP, KV, KT, tnow, Reff, P, RT)
% Expected SCF delay Dall(d,x) for every node x carrying a packet for d
% (Sec. 4.3). Nodes of one connected component share the same link-state
% view, so the view of the first member is used for the whole component.
if stamp == t, return; end
stamp = t; N = numel(comp); D = inf(1,N);
for c = unique(comp)
    mem = find(comp == c); rep = mem(1);
    if isinf(KT(rep,d)), continue; end
    age = (tnow - KT(rep,:))'; age(isinf(age)) = 0;
    Pv = reshape(KP(rep,:,:),N,3) + reshape(KV(rep,:,:),N,3).*age;
    Vv = reshape(KV(rep,:,:),N,3);
    Pv(mem,:) = pos(mem,:); Vv(mem,:) = vel(mem,:);
    pd = Pv(d,:); vd = Vv(d,:);
    if any(mem == d) || RT.nh(1,rep) > 0
        % connected to d: carrying helps only by reaching d directly
        x = mem(mem ~= d);
        dist = sqrt(sum((Pv(x,:) - pd).^2, 2));
        cx = contact_time(repmat(pd,numel(x),1) - Pv(x,:), repmat(vd,numel(x),1) - Vv(x,:), Reff, P.H);
        cx(dist <= Reff) = inf;
        D(x) = cx';
    else
        known = find(~isinf(KT(rep,:)));
        near = known(sqrt(sum((Pv(known,:) - pd).^2,2)) <= Reff);
        tg = unique([d near(~ismember(near, mem))]);
        nm = numel(mem); nt = numel(tg);
        [xm, yt] = meshgrid(mem, tg); xm = xm(:); yt = yt(:);
        cx = contact_time(Pv(yt,:) - Pv(xm,:), Vv(yt,:) - Vv(xm,:), Reff, P.H);
        cx = cx + P.tauBar*(yt ~= d);
        cx = reshape(cx, nt, nm);
        Dm = min(cx, [], 1);
        % no predicted contact within H: rank by distance (greedy fallback)
        fb = isinf(Dm);
        if any(fb)
            dist = sqrt(sum((Pv(mem(fb),:) - pd).^2, 2))';
            Dm(fb) = P.H + dist / P.vMax;
        end
        D(mem) = Dm;
    end
end
Dall(d,:) = D;
end

function u = best_carrier(i, DscfI, Drow, A, failed, occ, P)
% Neighbour predicted to deliver sooner by SCF (custody transfer), or 0.
u = 0;
nb = find(A(i,:) & ~failed(i,:) & occ' < P.B);
if isempty(nb), return; end
[best, k] = min(Drow(nb));
if best < DscfI - P.etaC, u = nb(k); end
end

function g = greedy_next(i, d, A, failed, occ, pos, KP, KV, KT, tnow, P)
% Greedy geographic forwarding: neighbour closest to d's estimated position.
g = 0;
[pd, ~] = view_of(i, d, pos, zeros(size(pos)), KP, KV, KT, tnow);
if any(isnan(pd)), return; end
best = norm(pos(i,:) - pd);
nb = find(A(i,:) & ~failed(i,:) & occ' < P.B);
for v = nb
    dv = norm(pos(v,:) - pd);
    if dv < best, best = dv; g = v; end
end
end

function q = prctile_(x, p)
if isempty(x), q = nan; return; end
x = sort(x(:)); k = max(1, ceil(p/100*numel(x))); q = x(k);
end
