%QUICK_TEST  One seed, all schemes, default parameters.
%   Run this first to check that everything works on your machine.
clear; clc;
here = fileparts(mfilename('fullpath')); if isempty(here), here = pwd; end
addpath(here);
P = config_default();
TR = gen_trace(P, 1);
fprintf('R_eff = %.0f m, %d packets generated, T_min(TH) = %g s\n', TR.Reff, numel(TR.pk_t), P.thTmin);
[~, ~, schemes] = sweep_defs();
fprintf('%-17s %6s %9s %9s %8s %8s\n', 'scheme', 'PDR', 'delay[s]', 'p95[s]', 'tx/del', 'buffer');
for s = 1:numel(schemes)
    t0 = tic;
    R = run_sim(P, TR, schemes{s}, 1);
    fprintf('%-17s %6.3f %9.1f %9.1f %8.2f %8.1f   (%.0f s)\n', schemes{s}, R.pdr, ...
        R.delayMean, R.delayP95, R.txPerDel, R.occMean, toc(t0));
end
