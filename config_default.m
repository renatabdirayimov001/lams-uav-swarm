function P = config_default()
%CONFIG_DEFAULT  Default parameters of the LAMS simulator.
%   All times in seconds, distances in metres, powers in dB/dBm.

% ---------------- simulation ----------------
P.dt        = 0.5;      % slot length
P.T_sim     = 600;      % total simulated time
P.T_warm    = 60;       % no traffic before this time
P.T_genEnd  = 480;      % no new packets after this time (lets packets drain)
P.T_upd     = 1.0;      % link-state / route update period (HELLO + TC)
P.maxRounds = 10;       % forwarding rounds per slot (multi-hop pipelining)

% ---------------- mobility (RPGM) ----------------
P.N      = 40;          % number of UAVs
P.G      = 4;           % number of groups
P.area   = [2000 2000]; % x-y extent
P.zLim   = [100 300];   % altitude range of group reference points
P.vMin   = 10;          % reference-point speed range
P.vMax   = 25;
P.rGroup = 150;         % max member offset from reference point
P.gmAlpha = 0.95;       % Gauss-Markov memory of member offset velocity (per slot)
P.gmSigma = 3;          % stationary std of member offset velocity (m/s)
P.pauseMax = 10;        % max pause at waypoints

% ---------------- channel (air-to-air) ----------------
P.Pt     = 5;           % transmit power (dBm)
P.PL0    = 40;          % path loss at 1 m, 2.4 GHz (dB)
P.alpha  = 2.5;         % path-loss exponent
P.sigmaS = 4;           % shadowing std (dB)
P.dCorr  = 50;          % shadowing decorrelation distance (m)
P.Kdb    = 10;          % Rician K-factor (dB)
P.W      = 2e6;         % bandwidth (Hz)
P.NF     = 7;           % noise figure (dB)
P.gammaTh = 5;          % SNR threshold (dB)
P.gapdB  = 3;           % SNR gap of the MCS (dB)
P.Rmax   = 6e6;         % max PHY rate (bit/s)
P.M      = 3;           % range margin used for R_eff (dB)
P.lqAlpha = 0.6;        % EWMA weight of link-quality estimate (HELLO hysteresis)
P.rateSNR = 'inst';    % airtime SNR: 'inst' (instantaneous) or 'max' (v1.0, optimistic)

% ---------------- traffic ----------------
P.F      = 10;          % number of flows
P.lambda = 1;           % bundles per second per flow
P.S      = 400e3;       % bundle size (bits) = 50 kB
P.B      = 100;         % buffer per node (bundles)
P.TTL    = 120;         % time-to-live (s)
P.dProc  = 0.005;       % per-hop processing delay (s)

% ---------------- schemes ----------------
P.T_hold = 2;           % MH: max time without route before drop
P.snwL   = 8;           % Spray-and-Wait copies (binary)
P.thTmin = 0.5;         % TH: fixed lifetime threshold, replaced by the tune_th result below
f_ = fullfile(fileparts(mfilename('fullpath')), 'results', 'th_best.mat');
if exist(f_, 'file'), L_ = load(f_); P.thTmin = L_.best; end

% ---------------- LAMS ----------------
P.H      = 60;          % prediction horizon
P.lam    = 30;          % failure penalty lambda (s)
P.eta    = 0.1;         % hysteresis
P.etaC   = 2;           % custody-transfer margin (s)
P.beta   = 5;           % weight of 1/T in lifetime-weighted route cost
P.tauBar = 1;           % intra-component delivery time after a contact
P.p0     = 0.2;         % SCF success prob. beyond the horizon
P.Wc     = 50;          % calibration window (link breaks)
P.calDefault = [0 1 2 0.3]; % [a b sigma0 kappa] before calibration / noCal
end
