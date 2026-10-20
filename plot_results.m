%PLOT_RESULTS  Figures for the paper from results/*.mat (PNG + EPS in figures/).
clear; close all;
here = fileparts(mfilename('fullpath')); if isempty(here), here = pwd; end
addpath(here);
inDir = fullfile(here, 'results'); figDir = fullfile(here, 'figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end

% fixed colour per scheme (colour-blind-checked categorical order) + marker
col = [42 120 214; 235 104 52; 27 175 122; 237 161 0; ...
       232 123 164; 0 131 0; 98 80 214; 227 73 72] / 255;
names = {'LAMS','LAMS-noAlt','LAMS-noCal','TH','RC','MH','EPI','SNW'};
mk    = {'o','s','d','^','v','>','<','p'};
ls    = {'-','--','--','-','-','-','-','-'};
labels = struct('N','Number of UAVs N', 'vMax','Max. speed v_{max} (m/s)', ...
                'lambda','Load \lambda_f (bundles/s per flow)', 'sigmaS','Shadowing \sigma_s (dB)');
ylab = struct('pdr','Packet delivery ratio', 'delayMean','Mean delay (s)', ...
              'delayP95','95th-percentile delay (s)', 'txPerDel','Transmissions per delivered packet', ...
              'occMean','Mean buffer occupancy (bundles)', 'shareSCF','Share delivered via SCF');

sweeps = {'N', {'pdr','delayMean','txPerDel'};
          'speed', {'pdr','delayMean'};
          'load', {'pdr','delayMean','occMean'};
          'shadow', {'pdr','delayMean'}};
fig = 0;
for s = 1:size(sweeps,1)
    f = fullfile(inDir, ['sweep_' sweeps{s,1} '.mat']);
    if ~exist(f, 'file'), continue; end
    S = load(f);
    for mm = sweeps{s,2}
        m = find(strcmp(S.metrics, mm{1}));
        fig = fig + 1; h = figure('Color','w','Position',[100 100 520 380]); hold on;
        for k = 1:numel(S.schemes)
            c = find(strcmp(names, S.schemes{k}));
            x = squeeze(S.M(m,k,:,:));               % seeds x values
            if isvector(x), x = x(:)'; end
            mu = zeros(1,size(x,2)); ci = mu;
            for v = 1:size(x,2)
                xv = x(isfinite(x(:,v)), v);
                mu(v) = mean(xv); ci(v) = 1.96*std(xv)/sqrt(max(numel(xv),1));
            end
            eb = errorbar(S.vals, mu, ci);
            set(eb, 'LineStyle', ls{c}, 'Color', col(c,:), 'Marker', mk{c}, ...
                'MarkerSize', 6, 'MarkerFaceColor', col(c,:), 'LineWidth', 1.5);
        end
        grid on; box on; set(gca, 'FontSize', 10); try, set(gca, 'GridAlpha', 0.15); catch, end
        xlabel(labels.(S.field)); ylabel(ylab.(mm{1}));
        if strcmp(mm{1}, 'pdr'), ylim([0 1]); end
        lg = legend(S.schemes, 'Location', 'best'); set(lg, 'FontSize', 8);
        try, set(lg, 'NumColumns', 2); catch, end
        base = fullfile(figDir, sprintf('fig_%s_%s', sweeps{s,1}, mm{1}));
        print(h, [base '.png'], '-dpng', '-r300');
        print(h, [base '.eps'], '-depsc');
    end
end

% ---- calibration figure ----
f = fullfile(inDir, 'calibration.mat');
if exist(f, 'file')
    load(f, 'cal');
    h = figure('Color','w','Position',[100 100 900 300]);
    for k = 1:numel(cal)
        subplot(1, numel(cal), k); hold on;
        sc = scatter(cal(k).pred, cal(k).obs, 6, col(1,:), 'filled');
        try, set(sc, 'MarkerFaceAlpha', 0.25); catch, end
        xx = linspace(0, max(cal(k).pred), 50); c = cal(k).fit;
        plot(xx, c(1) + c(2)*xx, '-', 'Color', col(2,:), 'LineWidth', 2);
        plot(xx, c(1) + c(2)*xx + (c(3) + c(4)*xx), '--', 'Color', col(2,:), 'LineWidth', 1);
        plot(xx, max(c(1) + c(2)*xx - (c(3) + c(4)*xx), 0), '--', 'Color', col(2,:), 'LineWidth', 1);
        plot(xx, xx, ':', 'Color', [0.4 0.4 0.4]);
        grid on; box on; try, set(gca, 'GridAlpha', 0.15); catch, end
        title(sprintf('\\sigma_s = %g dB', cal(k).sigmaS));
        xlabel('Predicted lifetime (s)'); if k == 1, ylabel('Realised lifetime (s)'); end
    end
    print(h, fullfile(figDir, 'fig_calibration.png'), '-dpng', '-r300');
    print(h, fullfile(figDir, 'fig_calibration.eps'), '-depsc');
end
fprintf('Figures written to %s\n', figDir);
