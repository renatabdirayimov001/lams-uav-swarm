function Mk = run_point(P, seed, schemes, metrics)
%RUN_POINT  Run all schemes on one trace (same seed) and collect metrics.
TR = gen_trace(P, seed);
Mk = nan(numel(metrics), numel(schemes));
for s = 1:numel(schemes)
    R = run_sim(P, TR, schemes{s}, seed);
    for m = 1:numel(metrics)
        Mk(m,s) = R.(metrics{m});
    end
end
end
