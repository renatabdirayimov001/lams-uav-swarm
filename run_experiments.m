function run_experiments(which)
%RUN_EXPERIMENTS  All parameter sweeps for the paper (resumable).
%   run_experiments            - run every sweep
%   run_experiments('N')       - run one sweep: 'N', 'speed', 'load', 'shadow'
%   run_experiments('cal')     - only the calibration data (Fig. 1)
%
%   Every point (one parameter value, all seeds, all schemes) is saved as
%   soon as it finishes in results/pt_<sweep>_<k>.mat. If the session stops,
%   simply call run_experiments again: finished points are skipped.
%   When all points of a sweep exist they are merged into sweep_<name>.mat
%   and results/summary.csv is rewritten.
if nargin < 1, which = 'all'; end
here = fileparts(mfilename('fullpath')); if isempty(here), here = pwd; end
addpath(here);
outDir = fullfile(here, 'results'); if ~exist(outDir, 'dir'), mkdir(outDir); end

seeds   = 1:10;                     % increase to 1:20 for the final version
schemes = {'LAMS','LAMS-noAlt','LAMS-noCal','TH','RC','MH','EPI','SNW'};
metrics = {'pdr','delayMean','delayP95','txPerDel','ctrlPerDel','occMean', ...
           'shareSCF','shareAlt','switchPerMin','dropTTL','dropBuf','dropNoRoute'};

% name, parameter field, values
sweeps = {
    'N',      'N',      [20 30 40 50 60];
    'speed',  'vMax',   [15 25 35];
    'load',   'lambda', [0.5 1 2 4];
    'shadow', 'sigmaS', [2 4 8];
};

for s = 1:size(sweeps,1)
    name = sweeps{s,1}; field = sweeps{s,2}; vals = sweeps{s,3};
    if ~(strcmp(which,'all') || strcmp(which,name)), continue; end
    f = fullfile(outDir, ['sweep_' name '.mat']);
    if exist(f, 'file'), fprintf('skip %s (done)\n', name); continue; end
    for v = 1:numel(vals)
        fp = fullfile(outDir, sprintf('pt_%s_%d.mat', name, v));
        if exist(fp, 'file'), fprintf('skip %s = %g (done)\n', field, vals(v)); continue; end
        P = config_default(); P.(field) = vals(v);
        if strcmp(field, 'N'), P.G = max(2, round(vals(v)/10)); end
        fprintf('[%s] %s = %g ... ', datestr(now,'HH:MM'), field, vals(v)); t0 = tic;
        Mv = nan(numel(metrics), numel(schemes), numel(seeds));
        parfor k = 1:numel(seeds)
            Mv(:,:,k) = run_point(P, seeds(k), schemes, metrics);
        end
        save(fp, 'Mv');
        fprintf('%.0f s (saved)\n', toc(t0));
    end
    % merge points into one sweep file
    M = nan(numel(metrics), numel(schemes), numel(seeds), numel(vals));
    for v = 1:numel(vals)
        S = load(fullfile(outDir, sprintf('pt_%s_%d.mat', name, v)));
        M(:,:,:,v) = S.Mv;
    end
    save(f, 'M', 'vals', 'field', 'schemes', 'metrics', 'seeds');
    fprintf('sweep %s finished\n', name);
    write_summary(outDir);
end

% calibration data (Fig. 1): predicted vs realised link lifetime
if strcmp(which,'all') || strcmp(which,'cal')
    f = fullfile(outDir, 'calibration.mat');
    if ~exist(f, 'file')
        sig = [2 4 8]; cal = struct('sigmaS', {}, 'pred', {}, 'obs', {}, 'fit', {});
        for k = 1:numel(sig)
            P = config_default(); P.sigmaS = sig(k);
            TR = gen_trace(P, 1); R = run_sim(P, TR, 'LAMS', 1);
            cal(k).sigmaS = sig(k); cal(k).pred = R.calPred; cal(k).obs = R.calObs;
            cal(k).fit = fit_calibration(R.calPred, R.calObs, P.calDefault);
        end
        save(f, 'cal');
        fprintf('calibration data saved\n');
    end
end

write_summary(outDir);
fprintf('Done. Run plot_results when all sweeps are finished.\n');
end
