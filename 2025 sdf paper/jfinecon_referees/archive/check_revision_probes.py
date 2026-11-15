"""
Exploratory probes behind revision_plan.tex (15 November 2026).

Data: stored global Preqin fund-level file (equal fund weights, commitment-scaled),
vintages <= 2010, Buyout and Venture Capital; one-factor model with alpha = 0;
monthly q5 market factor. These are design probes, not results for the paper.

Probe 1: is the "flat Buyout profile" under the cross-fund product criterion robust
         to vintage weights, a minimum fund count, and a predetermined origin/window?
Probe 2: what does a mean-moment (square-of-the-average) criterion give, and how
         large and persistent are vintage-level pooled pricing errors?
Probe 3: joint (alpha, beta) grid minimisers for squared loss and the product criterion.

Needs numpy and pandas only.
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
FACTOR_CSV = ROOT / "2025 sdf paper/presentation/nav_returns/q5_factors_monthly_2022_fallback.csv"
FUND_CSV = ROOT / "r_project/empirical/data_prepared_2026_b/preqin_cashflows_EW_NAV.csv"
HORIZONS = [0, 2.5, 5, 10, 15, 20, 30]
GRID = np.linspace(-0.5, 3.5, 401)


def load():
    q = pd.read_csv(FACTOR_CSV)
    q["Date"] = pd.date_range("1967-01-01", periods=len(q), freq="ME")
    q["RF"] = q.R_F / 100.0
    q["MKT"] = q.R_MKT / 100.0
    q = q.set_index("Date")[["RF", "MKT"]]
    fl = pd.read_csv(FUND_CSV, parse_dates=["Date"])
    return q, fl


def build(fl, q, asset, origin="sample"):
    """Per vintage: (rf, mkt, fund-by-month cash-flow matrix).

    origin="sample": grid from the earliest to the latest cash-flow date in the vintage
                     (what the paper's pipeline does).
    origin="pre":    grid from the December before the vintage year to the common
                     calendar end of the factor file (predetermined origin and window).
    """
    sub = fl[(fl.type == asset) & (fl.Vintage <= 2010)]
    out = {}
    for v, g in sub.groupby("Vintage"):
        fids = g["Fund.ID"].unique()
        if origin == "sample":
            d0, d1 = g.Date.min(), g.Date.max()
        else:
            d0, d1 = pd.Timestamp(year=int(v) - 1, month=12, day=31), q.index.max()
            g = g[g.Date >= d0]
        grid = pd.date_range(d0, d1, freq="ME")
        idx = {d: i for i, d in enumerate(grid)}
        fmap = {f: i for i, f in enumerate(fids)}
        cf = np.zeros((len(fids), len(grid)))
        for f, d, c in zip(g["Fund.ID"], g.Date, g.CF):
            if d in idx:
                cf[fmap[f], idx[d]] += c
        fq = q.reindex(grid)
        out[int(v)] = (fq.RF.values, fq.MKT.values, cf)
    return out


def fund_errors(vf, beta, mm, alpha=0.0, fix_a=None):
    """Horizon-scaled fund pricing errors A_{v,H} * PV_{v,i} for every vintage."""
    res = []
    for v, (rf, mk, cf) in vf.items():
        g = 1.0 + alpha + rf + beta * mk
        if np.any(g <= 0):
            return None
        ret = np.cumprod(g)
        pv = (cf / ret[None, :]).sum(1)
        if mm == 0:
            a = 1.0
        else:
            r = ret if fix_a is None else np.cumprod(1.0 + rf + fix_a * mk)
            a = r[: min(mm, len(r))].mean()
        res.append((v, a * pv))
    return res


def crit(vf, beta, mm, kind, wt="eq", minn=1, alpha=0.0):
    """kind="LS": squared constant-composition pooled error; kind="U": cross-fund products."""
    p = fund_errors(vf, beta, mm, alpha)
    if p is None:
        return np.inf
    tot = wsum = 0.0
    for _, e in p:
        n = len(e)
        if kind == "U":
            if n < max(2, minn):
                continue
            val = (e.sum() ** 2 - (e ** 2).sum()) / (n * (n - 1))
        else:
            if n < minn:
                continue
            val = e.mean() ** 2
        w = 1.0 if wt == "eq" else float(n)
        tot += w * val
        wsum += w
    return tot / wsum


def argmin(f):
    return GRID[int(np.argmin([f(b) for b in GRID]))]


def bisect(f, a, b):
    fa = f(a)
    for _ in range(50):
        m = 0.5 * (a + b)
        fm = f(m)
        if (fa < 0) == (fm < 0):
            a, fa = m, fm
        else:
            b = m
    return 0.5 * (a + b)


def roots(f):
    vals = np.array([f(b) for b in GRID])
    s = np.where(np.sign(vals[:-1]) != np.sign(vals[1:]))[0]
    return [round(bisect(f, GRID[i], GRID[i + 1]), 2) for i in s]


def pooled(vf, beta, mm, fix_a=None):
    p = fund_errors(vf, beta, mm, fix_a=fix_a)
    v = np.array([x[0] for x in p])
    e = np.array([x[1].mean() for x in p])
    n = np.array([len(x[1]) for x in p])
    return v, e, n


def main():
    q, fl = load()
    for asset in ["BO", "VC"]:
        vs = build(fl, q, asset, "sample")
        vp = build(fl, q, asset, "pre")
        print(f"\n=== {asset}, global EW, vintages <= 2010 ({len(vs)} vintages) ===")

        print("Probe 1. beta_hat by criterion")
        print("  H | LS pool | U eq | U n-wt | U n>=5 | U pre | LS pre | LS n>=5")
        for h in HORIZONS:
            mm = int(12 * h)
            row = [
                argmin(lambda b: crit(vs, b, mm, "LS")),
                argmin(lambda b: crit(vs, b, mm, "U")),
                argmin(lambda b: crit(vs, b, mm, "U", wt="n")),
                argmin(lambda b: crit(vs, b, mm, "U", minn=5)),
                argmin(lambda b: crit(vp, b, mm, "U")),
                argmin(lambda b: crit(vp, b, mm, "LS")),
                argmin(lambda b: crit(vs, b, mm, "LS", minn=5)),
            ]
            print(f"{h:5.1f} | " + " | ".join(f"{x:6.2f}" for x in row))

        print("Probe 2. roots of the equal-weighted mean of A_{v,H} PV_v (alpha = 0)")
        print("  H | A(beta) | A fixed at beta=1 | same, n_v weights")
        for h in HORIZONS:
            mm = int(12 * h)
            r1 = roots(lambda b: pooled(vs, b, mm)[1].mean())
            r2 = roots(lambda b: pooled(vs, b, mm, fix_a=1.0)[1].mean())

            def nw(b):
                _, e, n = pooled(vs, b, mm, fix_a=1.0)
                return (n * e).sum() / n.sum()

            print(f"{h:5.1f} | {r1} | {r2} | {roots(nw)}")

        b0 = argmin(lambda b: crit(vs, b, 0, "LS"))
        print(f"Vintage pooled PV at the H=0 squared-loss beta ({b0:.2f}), with within-vintage s.e.:")
        for v, (rf, mk, cf) in vs.items():
            x = (cf / np.cumprod(1.0 + rf + b0 * mk)[None, :]).sum(1)
            n = len(x)
            se = x.std(ddof=1) / np.sqrt(n) if n > 1 else float("nan")
            print(f"  {v}  n={n:3d}  PV={x.mean():+.3f}  se={se:.3f}")

        print("Probe 3. joint (alpha per month, beta) grid minimisers")
        for h in [0, 15]:
            mm = int(12 * h)
            best_u = (np.inf, None)
            best_ls = (np.inf, None)
            for a in np.arange(-0.02, 0.0201, 0.001):
                for b in np.arange(-0.5, 3.01, 0.05):
                    u = crit(vs, b, mm, "U", alpha=a)
                    ls = crit(vs, b, mm, "LS", alpha=a)
                    if u < best_u[0]:
                        best_u = (u, (round(float(a), 3), round(float(b), 2)))
                    if ls < best_ls[0]:
                        best_ls = (ls, (round(float(a), 3), round(float(b), 2)))
            print(f"  H={h}: U argmin={best_u[1]} value={best_u[0]:.4f} | LS argmin={best_ls[1]} value={best_ls[0]:.4f}")


if __name__ == "__main__":
    main()
