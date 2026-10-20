function write_summary(outDir)
%WRITE_SUMMARY  Mean and 95% CI of every metric to results/summary.csv.
fs = dir(fullfile(outDir, 'sweep_*.mat'));
fid = fopen(fullfile(outDir, 'summary.csv'), 'w');
fprintf(fid, 'sweep,param,value,scheme,metric,mean,ci95,nSeeds\n');
for f = 1:numel(fs)
    S = load(fullfile(outDir, fs(f).name));
    name = strrep(strrep(fs(f).name, 'sweep_', ''), '.mat', '');
    for v = 1:numel(S.vals)
        for s = 1:numel(S.schemes)
            for m = 1:numel(S.metrics)
                x = squeeze(S.M(m, s, :, v)); x = x(isfinite(x));
                n = numel(x); mu = mean(x);
                ci = 1.96*std(x)/sqrt(max(n,1));
                fprintf(fid, '%s,%s,%g,%s,%s,%.6g,%.6g,%d\n', name, S.field, ...
                        S.vals(v), S.schemes{s}, S.metrics{m}, mu, ci, n);
            end
        end
    end
end
fclose(fid);
end
