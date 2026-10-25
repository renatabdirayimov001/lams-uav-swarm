"""Tables, paired statistics and figures for the LAMS paper from results/*.mat.

Usage: python3 analyze.py <results_dir> <out_dir>
Writes: stats.json (all numbers quoted in the paper), paired.csv, fig_*.pdf/png
"""
import sys, os, json
import numpy as np
import scipy.io as sio
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)

NAME = {"LAMS-noAlt": "LAMS", "LAMS-noAlt-noCal": "LAMS-noCal", "LAMS": "LAMS-Alt",
        "TH": "TH", "RC": "RC", "MH": "MH", "EPI": "EPI", "SNW": "SnW"}
ORDER = ["LAMS-noAlt", "LAMS-noAlt-noCal", "LAMS", "TH", "RC", "MH", "EPI", "SNW"]
COL = {"LAMS-noAlt": "#2a78d6", "LAMS-noAlt-noCal": "#1baf7a", "LAMS": "#eb6834", "TH": "#eda100",
       "RC": "#e87ba4", "MH": "#008300", "EPI": "#6250d6", "SNW": "#e34948"}
MK = {"LAMS-noAlt": "o", "LAMS-noAlt-noCal": "d", "LAMS": "s", "TH": "^", "RC": "v", "MH": ">",
      "EPI": "<", "SNW": "p"}
LS = {"LAMS-noAlt-noCal": "--", "LAMS": "--"}


def cellstr(x):
    x = np.atleast_1d(np.squeeze(x))
    return [str(np.squeeze(v)) for v in x.ravel()]


def load_sweep(name):
    d = sio.loadmat(os.path.join(RES, f"sweep_{name}.mat"))
    M = np.asarray(d["M"], float)                      # metrics x schemes x seeds x vals
    if M.ndim == 3:
        M = M[..., None]
    return dict(M=M, vals=np.ravel(d["vals"]).astype(float), field=str(np.squeeze(d["field"])),
                schemes=cellstr(d["schemes"]), metrics=cellstr(d["metrics"]))


def get(S, metric, scheme):
    """seeds x vals"""
    return S["M"][S["metrics"].index(metric), S["schemes"].index(scheme)]


def tci(x):
    x = np.asarray(x, float); x = x[np.isfinite(x)]; n = len(x)
    if n < 2:
        return float(np.mean(x)), float("nan")
    return float(x.mean()), float(stats.t.ppf(0.975, n - 1) * x.std(ddof=1) / np.sqrt(n))


def paired(a, b, rel=False):
    """mean of a-b (or relative a/b-1 in %) with paired 95% t-interval and p-value"""
    a = np.asarray(a, float); b = np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b); a, b = a[ok], b[ok]
    d = (a / b - 1) * 100 if rel else (a - b)
    m, h = tci(d)
    p = float(stats.ttest_rel(a, b).pvalue) if len(d) > 1 and np.std(d) > 0 else float("nan")
    return dict(mean=m, lo=m - h, hi=m + h, p=p, n=int(len(d)))


sweeps = {k: load_sweep(k) for k in ["N", "speed", "load", "shadow"]
          if os.path.exists(os.path.join(RES, f"sweep_{k}.mat"))}
st = {"sweeps": {}, "paired": {}}

# ---------------- per-point means and paired comparisons ----------------
rows = []
for sw, S in sweeps.items():
    st["sweeps"][sw] = {"field": S["field"], "vals": S["vals"].tolist(), "mean": {}, "ci": {}}
    for sc in S["schemes"]:
        for m in S["metrics"]:
            X = get(S, m, sc)
            mc = [tci(X[:, v]) for v in range(X.shape[1])]
            st["sweeps"][sw]["mean"].setdefault(sc, {})[m] = [a for a, _ in mc]
            st["sweeps"][sw]["ci"].setdefault(sc, {})[m] = [b for _, b in mc]
    for ref in ["RC", "TH", "MH", "EPI", "SNW", "LAMS-noAlt-noCal", "LAMS"]:
        for m, rel in [("pdr", False), ("txPerDel", True), ("delayMean", True), ("delayP95", True),
                       ("occMean", True)]:
            A, B = get(S, m, "LAMS-noAlt"), get(S, m, ref)
            for v in range(A.shape[1]):
                r = paired(A[:, v], B[:, v], rel)
                if m == "pdr":
                    r = {k: (x * 100 if k in ("mean", "lo", "hi") else x) for k, x in r.items()}
                r.update(sweep=sw, value=float(S["vals"][v]), ref=ref, metric=m,
                         unit="%" if rel else "pp")
                rows.append(r)
st["paired"] = rows
import csv
with open(os.path.join(OUT, "paired.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["sweep", "value", "ref", "metric", "unit", "mean", "lo", "hi", "p", "n"])
    w.writeheader(); [w.writerow({k: r[k] for k in w.fieldnames}) for r in rows]

# sensitivity (v1.0 airtime model) at the default point
f = os.path.join(RES, "sens_rateMax.mat")
if os.path.exists(f) and "N" in sweeps:
    d = sio.loadmat(f); M = np.asarray(d["M"], float); scs = cellstr(d["schemes"]); mets = cellstr(d["metrics"])
    S = sweeps["N"]; v = list(S["vals"]).index(40.0)
    sens = {}
    for sc in scs:
        sens[sc] = {}
        for m in ["pdr", "delayMean", "txPerDel"]:
            old = M[mets.index(m), scs.index(sc)]; new = get(S, m, sc)[:, v]
            sens[sc][m] = dict(max=tci(old), inst=tci(new))
    # paired LAMS-RC and LAMS-TH under the v1.0 model
    for ref in ["RC", "TH"]:
        for m, rel in [("pdr", False), ("txPerDel", True)]:
            r = paired(M[mets.index(m), scs.index("LAMS-noAlt")], M[mets.index(m), scs.index(ref)], rel)
            sens.setdefault("paired_max", {})[f"{ref}_{m}"] = r
    st["sensitivity"] = sens

# TH tuning
f = os.path.join(RES, "th_best.mat")
if os.path.exists(f):
    d = sio.loadmat(f)
    st["th_tuning"] = dict(grid=np.ravel(d["grid_"]).tolist(), pdr=np.mean(d["pdr"], 1).tolist(),
                           best=float(np.squeeze(d["best"])))

# ---------------- figures ----------------
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.labelsize": 8.5,
                     "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 7.5,
                     "pdf.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.alpha": 0.25, "lines.linewidth": 1.4, "lines.markersize": 4.5})
XL = {"N": "Number of UAVs $N$", "vMax": "$v_{\\max}$ (m/s)",
      "lambda": "Load $\\lambda_f$ (bundles/s per flow)", "sigmaS": "Shadowing $\\sigma_s$ (dB)"}
YL = {"pdr": "Packet delivery ratio", "delayMean": "Mean end-to-end delay (s)",
      "txPerDel": "Tx per delivered bundle", "occMean": "Mean buffer occupancy (bundles)"}


def panel(ax, S, m, tag, logy=False):
    for sc in ORDER:
        if sc not in S["schemes"]:
            continue
        mu = st["sweeps"][[k for k, v in sweeps.items() if v is S][0]]["mean"][sc][m]
        ci = st["sweeps"][[k for k, v in sweeps.items() if v is S][0]]["ci"][sc][m]
        ax.errorbar(S["vals"], mu, yerr=ci, color=COL[sc], marker=MK[sc], ls=LS.get(sc, "-"),
                    capsize=2, elinewidth=0.8, label=NAME[sc])
    ax.set_xlabel(XL[S["field"]]); ax.set_ylabel(YL[m]); ax.set_xticks(S["vals"])
    if m == "pdr":
        ax.set_ylim(0, 1.02)
    if logy:
        ax.set_yscale("log")
        from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator
        t = [2, 3, 5, 10, 20, 30]; ax.yaxis.set_major_locator(FixedLocator(t))
        ax.yaxis.set_major_formatter(FixedFormatter([str(x) for x in t])); ax.yaxis.set_minor_locator(NullLocator())
    ax.text(-0.02, 1.03, tag, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", ha="right")


def multi(sw, ms, fname, logs=(), size=(7.16, 2.25), ncol=8):
    S = sweeps[sw]
    fig, axs = plt.subplots(1, len(ms), figsize=size)
    for a, m, t in zip(np.atleast_1d(axs), ms, "abc"):
        panel(a, S, m, f"({t})", logy=m in logs)
    h, l = np.atleast_1d(axs)[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=ncol, frameon=False, bbox_to_anchor=(0.5, 1.07 if ncol == 8 else 1.16),
               handlelength=2.2, columnspacing=1.0)
    fig.tight_layout(w_pad=1.2)
    for e in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"{fname}.{e}"), dpi=300, bbox_inches="tight")
    plt.close(fig)


if "N" in sweeps:
    multi("N", ["pdr", "delayMean", "txPerDel"], "fig_density", logs=("txPerDel",))
if "load" in sweeps:
    multi("load", ["pdr", "delayMean", "occMean"], "fig_load")
if "speed" in sweeps:
    multi("speed", ["pdr", "txPerDel"], "fig_speed", logs=("txPerDel",), size=(3.5, 2.0), ncol=4)

# losses at the default point (N = 40 of the N sweep)
if "N" in sweeps:
    S = sweeps["N"]; v = list(S["vals"]).index(40.0)
    fig, ax = plt.subplots(figsize=(3.5, 2.1))
    x = np.arange(len(ORDER)); bottom = np.zeros(len(ORDER))
    for m, lab, c in [("dropTTL", "TTL expiry", "#2a78d6"), ("dropBuf", "Source-buffer admission loss", "#eda100"),
                      ("dropNoRoute", "No route (MH)", "#6250d6")]:
        y = np.array([np.nanmean(get(S, m, sc)[:, v]) for sc in ORDER])
        ax.bar(x, y, 0.62, bottom=bottom, color=c, label=lab, edgecolor="white", linewidth=0.4)
        bottom += y
    ax.set_xticks(x); ax.set_xticklabels([NAME[s] for s in ORDER], rotation=30, ha="right")
    ax.set_ylabel("Share of generated bundles"); ax.set_ylim(0, max(0.40, np.ceil(bottom.max() * 20) / 20))
    ax.legend(loc="upper left", frameon=False, fontsize=7)
    ax.grid(axis="x", visible=False)
    st["losses_default"] = {NAME[s]: float(b) for s, b in zip(ORDER, bottom)}
    fig.tight_layout()
    for e in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"fig_losses.{e}"), dpi=300, bbox_inches="tight")
    plt.close(fig)

# calibration scatter
f = os.path.join(RES, "calibration.mat")
if os.path.exists(f):
    d = sio.loadmat(f, squeeze_me=True, struct_as_record=False)
    cal = np.atleast_1d(d["cal"])
    fig, axs = plt.subplots(1, len(cal), figsize=(7.16, 1.75), sharey=True)
    st["calibration"] = []
    for a, c, t in zip(np.atleast_1d(axs), cal, "abc"):
        pr, ob = np.ravel(c.pred).astype(float), np.ravel(c.obs).astype(float); fit = np.ravel(c.fit)
        a.scatter(pr, ob, s=3, color="#2a78d6", alpha=0.25, linewidths=0)
        xx = np.linspace(0, max(pr.max(), 1), 50)
        a.plot(xx, fit[0] + fit[1] * xx, color="#eb6834", lw=1.6)
        a.plot(xx, fit[0] + fit[1] * xx + fit[2] + fit[3] * xx, "--", color="#eb6834", lw=0.9)
        a.plot(xx, np.maximum(fit[0] + fit[1] * xx - fit[2] - fit[3] * xx, 0), "--", color="#eb6834", lw=0.9)
        a.plot(xx, xx, ":", color="0.4", lw=0.9)
        a.set_title(f"({t}) $\\sigma_s$ = {float(c.sigmaS):g} dB", fontsize=8.5)
        a.set_xlabel("Predicted lifetime $\\hat T^{(0)}$ (s)")
        st["calibration"].append(dict(sigmaS=float(c.sigmaS), n=int(len(pr)), fit=fit.tolist()))
    np.atleast_1d(axs)[0].set_ylabel("Observed lifetime (s)")
    fig.tight_layout()
    for e in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"fig_calibration.{e}"), dpi=300, bbox_inches="tight")
    plt.close(fig)

json.dump(st, open(os.path.join(OUT, "stats.json"), "w"), indent=1)
print("written to", OUT)
