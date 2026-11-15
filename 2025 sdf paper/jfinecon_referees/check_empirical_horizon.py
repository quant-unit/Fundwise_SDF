"""
Empirical verification and horizon-tilt decomposition for review_of_referee_reports.tex.

Verifies:
1. Local reproduction of Table 2 (BO) and Table 3 (VC) in sdf_paper.tex using
   stored North America VYP files (r_project/empirical/data_prepared_2026_b/)
   and the local monthly q5 factor file
   (2025 sdf paper/presentation/nav_returns/q5_factors_monthly_2022_fallback.csv).
2. Equivalence/difference between max.month == 0 (0 yr row) and max.month == 1
   (mathematical H=0 inception NPV, commented out as 0.1 yr in sdf_paper.tex).
3. Empirical decomposition of the falling beta(H) profile on North America BO & VC
   into:
   - Baseline squared loss: min_beta sum_v A_{v,H}(beta)^2 PV_v(beta)^2
   - Fixed-weight cross-vintage reweighting: min_beta sum_v A_{v,H}(1)^2 PV_v(beta)^2
   - Pure scale-Jacobian tilt: min_beta sum_v [A_{v,H}(beta)/A_{v,H}(1)]^2 PV_v(beta)^2
   alongside empirical mean S_bar_{v,H}(1) and cash-flow-weighted exit sensitivity S_exit(1).
4. Empirical decomposition on the stored Global Preqin fund-level dataset
   (preqin_cashflows_EW_NAV.csv, vintages <= 2010):
   - Reporter-mean VYP (current prep_preqin.R line 229)
   - Constant-composition VYP (sum(CF) / N_v)
   - Cross-fund product criterion U_{v,H}(beta)
"""

from pathlib import Path
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
FACTOR_CSV = ROOT / "2025 sdf paper/presentation/nav_returns/q5_factors_monthly_2022_fallback.csv"
DATA_DIR = ROOT / "r_project/empirical/data_prepared_2026_b"


def gss(f, a=-1.0, b=4.0, tol=1e-7):
    gr = (np.sqrt(5.0) + 1.0) / 2.0
    c = b - (b - a) / gr
    d = a + (b - a) / gr
    while abs(b - a) > tol:
        if f(c) < f(d):
            b = d
        else:
            a = c
        c = b - (b - a) / gr
        d = a + (b - a) / gr
    return 0.5 * (a + b)


def load_factors():
    q = pd.read_csv(FACTOR_CSV)
    q["Date"] = pd.date_range("1967-01-01", periods=len(q), freq="ME").strftime("%Y-%m-%d")
    q["RF"] = q["R_F"] / 100.0
    q["MKT"] = q["R_MKT"] / 100.0
    return q[["Date", "RF", "MKT"]]


def build_monthly_units(df_vyp, q, group_col="Fund.ID"):
    units = []
    for uid, grp in df_vyp.groupby(group_col):
        dmin = pd.to_datetime(grp["Date"].min())
        dmax = pd.to_datetime(grp["Date"].max())
        grid = pd.date_range(dmin, dmax, freq="ME").strftime("%Y-%m-%d")
        m = pd.DataFrame({"Date": grid}).merge(grp[["Date", "CF"]], on="Date", how="left")
        m["CF"] = m["CF"].fillna(0.0)
        m = m.merge(q, on="Date", how="left")
        units.append((uid, m))
    return units


def run_na_verification_and_decomposition(q):
    print("=== 1. North America VYP Reproduction & Scale-Jacobian Decomposition (V <= 2010) ===")
    horizons = [0, 1, 30, 60, 120, 150, 180, 210, 240, 300, 360]
    vyp_ew = pd.read_csv(DATA_DIR / "preqin_cashflows_EW_VYP_North America_NAV.csv")
    vyp_fw = pd.read_csv(DATA_DIR / "preqin_cashflows_FW_VYP_North America_NAV.csv")

    for asset in ["BO", "VC"]:
        units_ew = build_monthly_units(vyp_ew[(vyp_ew["type"] == asset) & (vyp_ew["Vintage"] <= 2010)], q)
        units_fw = build_monthly_units(vyp_fw[(vyp_fw["type"] == asset) & (vyp_fw["Vintage"] <= 2010)], q)

        print(f"\n--- {asset} (North America, V <= 2010) ---")
        print("max_m | H(yr) | EW Base | FW Base | EW Fixed-A(1) | EW Jac-Only | Mean S_bar(1) | Mean S_exit(1)")
        for mm in horizons:
            def solve_ls(units, mode="base"):
                def obj(beta):
                    sq = 0.0
                    for _, u in units:
                        ret = np.exp(np.cumsum(np.log(1.0 + u["RF"].values + beta * u["MKT"].values)))
                        dcf = np.sum(u["CF"].values / ret)
                        if mm == 0:
                            val = dcf
                        else:
                            imax = min(mm, len(ret))
                            if mode == "base":
                                val = dcf * np.mean(ret[:imax])
                            elif mode == "fixed_A":
                                ret1 = np.exp(np.cumsum(np.log(1.0 + u["RF"].values + 1.0 * u["MKT"].values)))
                                val = dcf * np.mean(ret1[:imax])
                            elif mode == "jac_only":
                                ret1 = np.exp(np.cumsum(np.log(1.0 + u["RF"].values + 1.0 * u["MKT"].values)))
                                val = dcf * (np.mean(ret[:imax]) / np.mean(ret1[:imax]))
                        sq += val ** 2
                    return sq
                return gss(obj, -1.0, 4.0)

            b_ew = solve_ls(units_ew, "base")
            b_fw = solve_ls(units_fw, "base")
            b_fix = solve_ls(units_ew, "fixed_A")
            b_jac = solve_ls(units_ew, "jac_only")

            sbars, sexits = [], []
            for _, u in units_ew:
                g1 = 1.0 + u["RF"].values + 1.0 * u["MKT"].values
                ret1 = np.exp(np.cumsum(np.log(g1)))
                s_t = np.cumsum(u["MKT"].values / g1)
                if mm == 0:
                    sbar = 0.0
                else:
                    imax = min(mm, len(ret1))
                    sbar = np.sum(ret1[:imax] * s_t[:imax]) / np.sum(ret1[:imax])
                pos_cf = np.maximum(u["CF"].values, 0.0) / ret1
                sexit = np.sum(pos_cf * s_t) / np.sum(pos_cf)
                sbars.append(sbar)
                sexits.append(sexit)

            h_label = 0.0 if mm == 0 else (0.083 if mm == 1 else mm / 12.0)
            print(
                f"{mm:5d} | {h_label:5.2f} | {b_ew:7.3f} | {b_fw:7.3f} | "
                f"{b_fix:13.3f} | {b_jac:11.3f} | {np.mean(sbars):13.3f} | {np.mean(sexits):14.3f}"
            )


def run_global_pool_and_product_check(q):
    print("\n=== 2. Global Fund-Level Sample (EW, V <= 2010): Reporter-Mean vs Constant Pool vs U_{v,H} ===")
    fl = pd.read_csv(DATA_DIR / "preqin_cashflows_EW_NAV.csv")
    horizons = [0, 30, 60, 120, 150, 180, 240, 300, 360]

    for asset in ["BO", "VC"]:
        sub = fl[(fl["type"] == asset) & (fl["Vintage"] <= 2010)].copy()
        rep_vyp = sub.groupby(["Vintage", "Date"], as_index=False)["CF"].mean()
        nv_map = sub.groupby("Vintage")["Fund.ID"].nunique().to_dict()
        const_vyp = sub.groupby(["Vintage", "Date"], as_index=False)["CF"].sum()
        const_vyp["CF"] = const_vyp["CF"] / const_vyp["Vintage"].map(nv_map)

        u_rep = build_monthly_units(rep_vyp, q, group_col="Vintage")
        u_const = build_monthly_units(const_vyp, q, group_col="Vintage")

        vint_funds = {}
        for v, grp in sub.groupby("Vintage"):
            fids = grp["Fund.ID"].unique()
            if len(fids) < 2:
                continue
            dmin = pd.to_datetime(grp["Date"].min())
            dmax = pd.to_datetime(grp["Date"].max())
            grid = pd.date_range(dmin, dmax, freq="ME").strftime("%Y-%m-%d")
            m_q = pd.DataFrame({"Date": grid}).merge(q, on="Date", how="left")
            cf_mat = np.zeros((len(fids), len(grid)))
            date_to_idx = {d: idx for idx, d in enumerate(grid)}
            for i, (_, fgrp) in enumerate(grp.groupby("Fund.ID")):
                for d, cf_val in zip(fgrp["Date"], fgrp["CF"]):
                    if d in date_to_idx:
                        cf_mat[i, date_to_idx[d]] += cf_val
            vint_funds[v] = (m_q["RF"].values, m_q["MKT"].values, cf_mat)

        print(f"\n--- Global {asset} (EW, V <= 2010, {len(vint_funds)} vintages with n_v >= 2) ---")
        print("H(yr) | Reporter-Mean VYP | Constant-Pool VYP | Cross-Fund Product U_{v,H}")
        for mm in horizons:
            def solve_vyp(units):
                def obj(beta):
                    sq = 0.0
                    for _, u in units:
                        ret = np.exp(np.cumsum(np.log(1.0 + u["RF"].values + beta * u["MKT"].values)))
                        dcf = np.sum(u["CF"].values / ret)
                        val = dcf if mm == 0 else dcf * np.mean(ret[:min(mm, len(ret))])
                        sq += val ** 2
                    return sq
                return gss(obj, -1.0, 4.0)

            def obj_u(beta):
                tot = 0.0
                for _, (rf_v, mkt_v, cf_mat) in vint_funds.items():
                    nv = cf_mat.shape[0]
                    ret = np.exp(np.cumsum(np.log(1.0 + rf_v + beta * mkt_v)))
                    dcfs = np.sum(cf_mat / ret[None, :], axis=1)
                    mult = 1.0 if mm == 0 else np.mean(ret[:min(mm, len(ret))])
                    e = dcfs * mult
                    tot += (np.sum(e) ** 2 - np.sum(e ** 2)) / (nv * (nv - 1))
                return tot

            bgrid = np.linspace(-0.5, 3.5, 801)
            b_u = bgrid[int(np.argmin([obj_u(b) for b in bgrid]))]
            print(f"{mm/12.0:5.1f} | {solve_vyp(u_rep):17.3f} | {solve_vyp(u_const):17.3f} | {b_u:26.3f}")


if __name__ == "__main__":
    factors = load_factors()
    run_na_verification_and_decomposition(factors)
    run_global_pool_and_product_check(factors)
