function write_summary(outDir)
%WRITE_SUMMARY  Mean and 95% Student-t CI of every metric to results/summary.csv.
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
                tq = [12.706 4.303 3.182 2.776 2.571 2.447 2.365 2.306 2.262 2.228 2.201 2.179 2.160 2.145 2.131 2.120 2.110 2.101 2.093];
                if n >= 2 && n <= 20, t = tq(n-1); else, t = 1.96; end
                ci = t*std(x)/sqrt(max(n,1));   % 95% Student-t CI
                fprintf(fid, '%s,%s,%g,%s,%s,%.6g,%.6g,%d\n', name, S.field, ...
                        S.vals(v), S.schemes{s}, S.metrics{m}, mu, ci, n);
            end
        end
    end
end
fclose(fid);
end
