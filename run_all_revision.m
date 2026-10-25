%RUN_ALL_REVISION  Everything needed for the revised paper, in order.
%   1) tune the TH threshold   2) all sweeps (10 seeds, 8 schemes)
%   3) calibration data (Fig. 1)   4) v1.0 airtime sensitivity check
%   Resumable: if MATLAB stops, run this script again; finished parts are skipped.
%   When it prints ALL DONE, zip the whole "results" folder and send it.
here = fileparts(mfilename('fullpath')); if isempty(here), here = pwd; end
addpath(here); cd(here);
if ~exist(fullfile(here, 'results', 'th_best.mat'), 'file'), tune_th(); end
run_experiments('all');
run_sensitivity();
fprintf('\nALL DONE. Zip the folder %s and send it.\n', fullfile(here, 'results'));
