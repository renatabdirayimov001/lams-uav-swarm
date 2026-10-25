function run_worker(wid, nw)
%RUN_WORKER  Run every nw-th point of all sweeps (worker wid of nw, 0-based).
%   Start nw MATLAB/Octave sessions with wid = 0..nw-1 to use several cores
%   without the Parallel Computing Toolbox. Finished points are skipped, so
%   the call can be repeated after an interruption. Merge with run_experiments.
here = fileparts(mfilename('fullpath')); if isempty(here), here = pwd; end
addpath(here); outDir = fullfile(here, 'results'); if ~exist(outDir,'dir'), mkdir(outDir); end
[sweeps, seeds, schemes, metrics] = sweep_defs();
jobs = {};
for s = 1:size(sweeps,1)
    for v = 1:numel(sweeps{s,3}), for k = 1:numel(seeds), jobs{end+1} = [s v k]; end, end %#ok<AGROW>
end
for q = 1:numel(jobs)
    if mod(q-1, nw) ~= wid, continue; end
    s = jobs{q}(1); v = jobs{q}(2); k = jobs{q}(3);
    name = sweeps{s,1}; field = sweeps{s,2}; val = sweeps{s,3}(v);
    fp = fullfile(outDir, sprintf('seed_%s_%d_%d.mat', name, v, seeds(k)));
    if exist(fp, 'file'), continue; end
    P = config_default(); P.(field) = val;
    if strcmp(field, 'N'), P.G = max(2, round(val/10)); end
    t0 = tic; Mk = run_point(P, seeds(k), schemes, metrics); %#ok<NASGU>
    save(fp, 'Mk', '-v7');
    fprintf('[%s] %s=%g seed %d: %.0f s\n', datestr(now,'HH:MM'), field, val, seeds(k), toc(t0));
end
end
