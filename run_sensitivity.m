function run_sensitivity()
%RUN_SENSITIVITY  Default scenario with the v1.0 airtime model (rateSNR = 'max').
%   Compares against the default point of the main sweeps (rateSNR = 'inst').
here = fileparts(mfilename('fullpath')); if isempty(here), here = pwd; end
addpath(here); outDir = fullfile(here, 'results');
[~, seeds, ~, metrics] = sweep_defs();
schemes = {'LAMS-noAlt','TH','RC','MH'};
f = fullfile(outDir, 'sens_rateMax.mat'); if exist(f, 'file'), return; end
M = nan(numel(metrics), numel(schemes), numel(seeds));
for k = 1:numel(seeds)
    P = config_default(); P.rateSNR = 'max';
    M(:,:,k) = run_point(P, seeds(k), schemes, metrics);
    fprintf('sensitivity seed %d done\n', seeds(k));
end
save(f, 'M', 'schemes', 'metrics', 'seeds', '-v7');
end
