function best = tune_th()
%TUNE_TH  Grid search of the fixed threshold T_min of the TH baseline.
%   Tuning seeds (101-103) differ from the test seeds (1-10). The best value
%   (highest mean PDR at default parameters) is saved to results/th_best.mat,
%   which config_default.m reads, so all later runs use it automatically.
here = fileparts(mfilename('fullpath')); if isempty(here), here = pwd; end
addpath(here); outDir = fullfile(here, 'results'); if ~exist(outDir,'dir'), mkdir(outDir); end
grid_ = [0 0.1 0.25 0.5 1 2 4 8]; seeds = 101:103;
pdr = zeros(numel(grid_), numel(seeds));
for k = 1:numel(seeds)
    P = config_default(); TR = gen_trace(P, seeds(k));
    for g = 1:numel(grid_)
        P.thTmin = grid_(g);
        R = run_sim(P, TR, 'TH', seeds(k)); pdr(g,k) = R.pdr;
    end
    fprintf('tuning seed %d done\n', seeds(k));
end
m = mean(pdr, 2);
for g = 1:numel(grid_), fprintf('T_min = %5.2f s : PDR = %.4f\n', grid_(g), m(g)); end
[~, b] = max(m); best = grid_(b);
save(fullfile(outDir, 'th_best.mat'), 'best', 'grid_', 'pdr', '-v7');
fprintf('Best T_min = %g s (saved to results/th_best.mat)\n', best);
end
