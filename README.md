# LAMS: Lifetime-Aware Mode Selection for UAV Swarms

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23632607.svg)](https://doi.org/10.5281/zenodo.23632607)

MATLAB simulator and results for the paper
**"Lifetime-Aware Adaptive Switching Between Multi-Hop Forwarding and Store-Carry-Forward in UAV Swarm Networks"** (under review).

LAMS decides, for every queued packet, whether a UAV should forward it over the current multi-hop route or carry it (store-carry-forward). The decision compares the backlog drain time with the predicted route lifetime, uses an online-calibrated model of the prediction error, and estimates the carrying delay from predicted contacts.

## Requirements
MATLAB R2019b or newer (or GNU Octave 8). No toolboxes are needed; the Parallel Computing Toolbox (optional) speeds up `run_experiments`. The code also runs in GNU Octave 8, more slowly.

## Quick start
```matlab
quick_test          % one seed, all schemes
run_all_revision    % TH tuning + all sweeps (10 seeds) + calibration data + sensitivity check; resumable
```
`run_all_revision` writes everything to `results/`. To use several cores without the
Parallel Computing Toolbox, start `run_worker(0,2)` and `run_worker(1,2)` in two sessions
first, then call `run_all_revision`.

## Scheme names: simulator vs. paper
| Simulator name | Name in the paper | Description |
|---|---|---|
| `LAMS-noAlt` | **LAMS** | Proposed method (route / carry) |
| `LAMS` | LAMS-Alt | Ablation: with alternative routes |
| `LAMS-noAlt-noCal` | LAMS-noCal | Ablation: LAMS without calibration |
| `TH` | TH | Fixed lifetime threshold |
| `RC` | RC | Reactive cascade in the style of R-T-OLSR |
| `MH` | MH | Pure multi-hop (AODV/OLSR style) |
| `EPI` | EPI | Epidemic routing |
| `SNW` | SnW | Binary Spray-and-Wait (L = 8) |

## Changes in v1.1
- Airtime is computed from the instantaneous SNR (`P.rateSNR = 'inst'`); v1.0 used the larger of instantaneous and mean SNR (`'max'`, kept for the sensitivity check).
- New ablation `LAMS-noAlt-noCal` (LAMS without online calibration).
- TH threshold tuned on the extended grid {0, 0.1, 0.25, 0.5, 1, 2, 4, 8} s.
- Per-seed results are kept, so paired comparisons can be computed.

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
| `results/summary.csv` | Mean and 95% Student-t CI of every metric (10 seeds) |
| `results/sweep_*.mat`, `calibration.mat`, `sens_rateMax.mat`, `th_best.mat` | Per-seed results of v1.1 (produced with MATLAB) |
| `results/paired.csv` | Paired differences of LAMS vs. every other scheme with 95% t-intervals and p-values |
| `analysis/analyze.py` | Python script that produces the figures, `paired.csv` and all numbers of the paper from `results/*.mat` |
| `figures/` | Figures of the paper (`fig_lams_architecture.png`: framework overview; icons and city illustration by Freepik) |

## Reproducing the figures
```bash
pip install numpy scipy matplotlib
python3 analysis/analyze.py results figures
```

## Citation
If you use this code, please cite the archived version: R. Abdirayimov, *LAMS: Lifetime-Aware Mode Selection for UAV Swarms, Simulator and Results*, Zenodo, v1.1.0, 2026. doi:10.5281/zenodo.23632607. The paper reference will be added after publication.

## License
MIT
