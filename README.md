# LAMS: Lifetime-Aware Mode Selection for UAV Swarms

MATLAB simulator and results for the paper
**"Lifetime-Aware Adaptive Switching Between Multi-Hop Forwarding and Store-Carry-Forward in UAV Swarm Networks"** (under review).

LAMS decides, for every queued packet, whether a UAV should forward it over the current multi-hop route or carry it (store-carry-forward). The decision compares the backlog drain time with the predicted route lifetime, uses an online-calibrated model of the prediction error, and estimates the carrying delay from predicted contacts.

## Requirements
MATLAB R2019b or newer. No toolboxes are needed; the Parallel Computing Toolbox (optional) speeds up `run_experiments`. The code also runs in GNU Octave 8, more slowly.

## Quick start
```matlab
quick_test          % one seed, all schemes (~1 min)
tune_th             % tunes the threshold of the TH baseline (result: 0.5 s)
run_experiments     % all sweeps, 10 seeds; resumable
plot_results        % figures into figures/
```

## Scheme names: simulator vs. paper
| Simulator name | Name in the paper | Description |
|---|---|---|
| `LAMS-noAlt` | **LAMS** | Proposed method (route / carry) |
| `LAMS` | LAMS-Alt | Ablation: with alternative routes |
| `LAMS-noCal` | LAMS-Alt-noCal | Ablation: without calibration |
| `TH` | TH | Fixed lifetime threshold |
| `RC` | RC | Reactive cascade in the style of R-T-OLSR |
| `MH` | MH | Pure multi-hop (AODV/OLSR style) |
| `EPI` | EPI | Epidemic routing |
| `SNW` | SnW | Binary Spray-and-Wait (L = 8) |

## Files
| File | Purpose |
|---|---|
| `config_default.m` | All parameters (Table 2 of the paper) |
| `gen_trace.m` | RPGM group mobility, air-to-air channel (path loss, Gudmundson shadowing, Rician fading), traffic |
| `run_sim.m` | Discrete-time simulator and all forwarding schemes |
| `route_trees.m` | Min-hop, max-lifetime and lifetime-weighted routing trees |
| `link_lifetime.m`, `contact_time.m` | Predicted link lifetime and contact time |
| `fit_calibration.m` | Online calibration of the lifetime-error model |
| `run_experiments.m`, `run_point.m`, `write_summary.m` | Parameter sweeps and CSV summary |
| `plot_results.m`, `quick_test.m`, `tune_th.m` | Figures, smoke test, baseline tuning |
| `results/summary.csv` | Mean and 95% CI of every metric (10 seeds) used in the paper |
| `figures/` | Figures of the paper |

## Citation
If you use this code, please cite the paper (reference will be added after publication).

## License
MIT
