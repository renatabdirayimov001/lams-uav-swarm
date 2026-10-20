%TUNE_TH  Grid search of the fixed threshold T_min of the TH baseline.
%   The best value (highest mean PDR at default parameters) should be
%   written into config_default.m (P.thTmin) before run_experiments.
clear; clc;
here = fileparts(mfilename('fullpath')); if isempty(here), here = pwd; end
addpath(here);
grid_ = [0.5 1 2 4 8]; seeds = 101:103;     % tuning seeds differ from test seeds
pdr = zeros(numel(grid_), numel(seeds));
for k = 1:numel(seeds)
    P = config_default(); TR = gen_trace(P, seeds(k));
    for g = 1:numel(grid_)
        P.thTmin = grid_(g);
        R = run_sim(P, TR, 'TH', seeds(k)); pdr(g,k) = R.pdr;
    end
end
m = mean(pdr, 2);
for g = 1:numel(grid_), fprintf('T_min = %4.1f s : PDR = %.3f\n', grid_(g), m(g)); end
[~, b] = max(m); fprintf('Best T_min = %g s -> set P.thTmin in config_default.m\n', grid_(b));
