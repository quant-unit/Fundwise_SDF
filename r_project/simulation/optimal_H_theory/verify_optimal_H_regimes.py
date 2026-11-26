#!/usr/bin/env python3
"""
Exact analytical calculations and nonlinear cash-flow Monte Carlo verification
for the Optimal-H Theory across three dependence regimes:
  Regime 1 : i.i.d. fund/deal errors (c_v = 0, c_{v,w} = 0)
  Regime 2 : Independent static vintage shocks W_v (c_v > 0, rho_W = 0)
  Regime 3a: Persistent AR(1) cross-vintage shocks W_v (c_v > 0, rho_W in {0.6, 0.8})
  Regime 3b: Overlapping calendar-time latent PE shock e_u (i.i.d. & AR(1))

Also verifies:
  - Exact zero-bias horizon table H_0*(T_eff, g, sigma_r)
  - Effective duration T_eff across shock mechanisms (static W_v, calendar e_u, L2 idio)
  - Leave-band-out cross-vintage U-statistic U_cross_{0,m} spectrum & boundary failure

Outputs CSV data tables to ./data/ for TikZ/PGFPlots and prints summary tables.
"""

import os
import json
import numpy as np

OUT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(OUT_DIR, "data")
os.makedirs(DATA_DIR, exist_ok=True)


# =============================================================================
# 1. Exact Analytical Duration-Matching Functions (Deterministic & Stochastic)
# =============================================================================

def tau_bar_det(H: float, r: float) -> float:
    """Continuous-time return-weighted mean compounding date under constant r."""
    if H <= 0.0:
        return 0.0
    if abs(r) < 1e-10:
        return 0.5 * H
    x = r * H
    if x > 500.0:
        return H - 1.0 / r
    return H - 1.0 / r + H / (np.exp(x) - 1.0)


def solve_H_det(T_eff: float, g: float) -> float:
    """Solve tau_bar_det(H, ln(1+g)) = T_eff by bisection."""
    r = np.log(1.0 + g)
    lo, hi = 1e-6, max(4.0 * T_eff, 60.0)
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if tau_bar_det(mid, r) < T_eff:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def tau_bar_stoch(H: float, g: float, sigma_r: float, n_grid: int = 200) -> float:
    """
    Exact continuous-time stochastic lognormal compounding date solving
      E[A_H^2 * (S_bar_H - T_eff)] = 0
    where r_u ~ N(mu_r, sigma_r^2) with arithmetic return E[exp(r_1)] = 1 + g,
    so mu_r = ln(1+g) - 0.5 * sigma_r^2 and
      E[R_s R_t] = exp( ln(1+g)*(s+t) + sigma_r^2 * min(s,t) ).
    """
    if H <= 0.0:
        return 0.0
    r_arith = np.log(1.0 + g)
    if sigma_r <= 1e-10:
        return tau_bar_det(H, r_arith)
    # Gauss-Legendre / midpoint quadrature on [0, H]
    dt = H / n_grid
    t = (np.arange(n_grid) + 0.5) * dt
    S, T = np.meshgrid(t, t, indexing="ij")
    log_E = r_arith * (S + T) + (sigma_r ** 2) * np.minimum(S, T)
    log_E -= np.max(log_E)
    W = np.exp(log_E)
    return float(np.sum(W * T) / np.sum(W))


def solve_H_stoch(T_eff: float, g: float, sigma_r: float) -> float:
    """Solve tau_bar_stoch(H, g, sigma_r) = T_eff by bisection."""
    if sigma_r <= 1e-10:
        return solve_H_det(T_eff, g)
    lo, hi = 1e-4, max(3.5 * T_eff, 50.0)
    for _ in range(55):
        mid = 0.5 * (lo + hi)
        if tau_bar_stoch(mid, g, sigma_r, n_grid=140) < T_eff:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def compute_T_eff_mechanisms(I_max=5.0, h_min=0.25, h_max=11.5, g=0.133, sigma_m=0.154):
    """
    Compute exact covariance-weighted payoff duration T_eff and zero-bias H_0*
    under the core deal-timing calibration a ~ U[0, I_max], h ~ U[h_min, h_max].
    """
    a_c = np.linspace(0.0, I_max, 21)
    h_c = np.linspace(h_min, h_max, 25)
    Ac, Hc = np.meshgrid(a_c, h_c, indexing="ij")
    ac, hc, bc = Ac.ravel(), Hc.ravel(), (Ac + Hc).ravel()

    r_disc = np.log(1.0 + g)
    pv_w = np.exp(-r_disc * ac)  # D_{0,a}(theta_0)

    # 1. Static vintage shock W_v (rank-1 across all deals in vintage v)
    W_static = np.outer(pv_w, pv_w)
    T_static = float(np.sum(W_static * bc[None, :]) / np.sum(W_static))

    # 2. Idiosyncratic Brownian/multiplicative deal shocks (L2, diagonal, Var propto h)
    w_idio = (pv_w ** 2) * hc
    T_idio = float(np.sum(w_idio * bc) / np.sum(w_idio))

    # 3. Overlapping calendar-time i.i.d. PE shock e_u
    overlap = np.maximum(0.0, np.minimum(bc[:, None], bc[None, :]) - np.maximum(ac[:, None], ac[None, :]))
    W_cal_iid = np.outer(pv_w, pv_w) * overlap
    T_cal_iid = float(np.sum(W_cal_iid * bc[None, :]) / np.sum(W_cal_iid))

    # 4. Overlapping calendar-time AR(1) PE shock e_u (annual corr 0.5 and 0.8)
    dt = 0.25
    t_steps = np.arange(0.0, I_max + h_max + dt, dt)
    M = ((t_steps[None, :] > ac[:, None]) & (t_steps[None, :] <= bc[:, None])).astype(float) * dt
    T_cal_ar = {}
    for rho_yr in [0.5, 0.8]:
        lam = -np.log(rho_yr)
        K_time = np.exp(-lam * np.abs(t_steps[:, None] - t_steps[None, :]))
        Sigma_ar = np.outer(pv_w, pv_w) * (M @ K_time @ M.T)
        T_cal_ar[rho_yr] = float(np.sum(Sigma_ar * bc[None, :]) / np.sum(Sigma_ar))

    rows = [
        ("Static vintage shock ($W_v$)", T_static),
        ("Calendar-time shock ($e_u$, i.i.d.)", T_cal_iid),
        ("Calendar-time shock ($e_u$, AR(1) $\\rho_e=0.5$)", T_cal_ar[0.5]),
        ("Calendar-time shock ($e_u$, AR(1) $\\rho_e=0.8$)", T_cal_ar[0.8]),
        ("Idiosyncratic deal shocks (L2, $c_v=0$)", T_idio),
    ]
    out = []
    for name, Te in rows:
        out.append({
            "mechanism": name,
            "T_eff": Te,
            "H_flat": 2.0 * Te,
            "H_det": solve_H_det(Te, g),
            "H_stoch": solve_H_stoch(Te, g, sigma_m),
        })
    return out


# =============================================================================
# 2. Nonlinear Cash-Flow Monte Carlo Across Regimes 1, 2, 3a, 3b
# =============================================================================

def simulate_regimes_mc(
    V=25,
    n_funds=15,
    n_deals=15,
    n_reps=600,
    seed=20261126,
    sigma_W=0.35,
    sigma_idio_yr=0.45,
    g_ann=0.133,
    rf_ann=0.030,
    sigma_mkt_ann=0.154,
):
    """
    Quarterly nonlinear cash-flow panel simulation matching Preqin core deal geometry:
      - V vintages, annual origins y_v = 0, ..., V-1 (in quarters: t_v^0 = 4*v)
      - n_funds funds per vintage, n_deals deals per fund
      - Deal entry a ~ U{0, ..., 20} quarters (0..5 yr), hold h ~ U{1, ..., 46} quarters (0.25..11.5 yr)
      - True parameters: alpha_0 = 0, beta_0 = 1
      - Reciprocal linear SDF: G_u(alpha, beta) = 1 + rf_q + alpha_q + beta * F_q
      - Estimate alpha (fixing beta=1) AND beta (fixing alpha=0) on U_{v,H} (q_v=1) and L2 VYP
        across H in [0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 25, 30] yr.
      - Also evaluate leave-band-out cross-vintage U_cross_{0,m} (m=0, 5, 7) on exact same draws.
    """
    rng = np.random.default_rng(seed)

    H_grid = np.array([0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5, 20.0, 25.0, 30.0])
    H_q = np.round(4 * H_grid).astype(int)
    max_H_q = int(H_q.max())

    # Quarterly factor parameters (exact match to annual mean g_ann=13.3%, rf=3.0%, vol=15.4%)
    rf_q = (1.0 + rf_ann) ** 0.25 - 1.0
    mu_G_q = (1.0 + g_ann) ** 0.25
    mu_F_q = mu_G_q - (1.0 + rf_q)
    sig_F_q = sigma_mkt_ann / 2.0

    # Total calendar quarters needed: last vintage origin 4*(V-1) + max(max_H_q, 66)
    max_deal_q = 20 + 46  # 16.5 yr
    T_total = 4 * (V - 1) + max(max_H_q, max_deal_q) + 2

    # Fix a realistic market factor path (matching the fixed historical q5 path design in sdf_paper.tex!)
    # Draw one clean representative path with sample mean/vol matching US 1985-2021
    F_raw = rng.normal(0.0, 1.0, size=T_total)
    F_raw = (F_raw - F_raw.mean()) / F_raw.std()
    F_q = mu_F_q + sig_F_q * F_raw
    # Ensure gross return at beta=1 is strictly positive everywhere
    G0_q = 1.0 + rf_q + F_q
    assert np.all(G0_q > 0.25), "Unexpected extreme return draw"

    # Precompute cumulative log discount factors and A_{v,H}(theta) on fine candidate grids!
    # Using exact 1D grid search + parabolic refinement for blazing speed and zero local-optima bugs!
    # (a) Alpha estimation (fixing beta=1): annualized alpha in [-15%, +25%] (monthly [-1.25%, +2.08%])
    alpha_ann_grid = np.linspace(-0.15, 0.24, 391)  # 0.1% annual resolution
    alpha_q_grid = alpha_ann_grid / 4.0             # quarterly additive alpha
    # Check floor G_u >= 0.08
    G_alpha = G0_q[None, :] + alpha_q_grid[:, None]  # (nA, T_total)
    assert np.all(G_alpha > 0.08)
    # Cumulative gross return from t=0: R_cum_a[k, t] = prod_{u=0}^{t-1} G_u(alpha_k), R_cum_a[:, 0] = 1
    log_R_a = np.zeros((len(alpha_ann_grid), T_total + 1))
    log_R_a[:, 1:] = np.cumsum(np.log(G_alpha), axis=1)

    # (b) Beta estimation (fixing alpha=0): beta in [-1.0, 4.0]
    beta_grid = np.linspace(-1.0, 4.0, 501)  # 0.01 resolution
    G_beta = 1.0 + rf_q + beta_grid[:, None] * F_q[None, :]
    valid_beta = np.all(G_beta >= 0.08, axis=1)
    beta_grid = beta_grid[valid_beta]
    G_beta = G_beta[valid_beta]
    log_R_b = np.zeros((len(beta_grid), T_total + 1))
    log_R_b[:, 1:] = np.cumsum(np.log(G_beta), axis=1)

    # True log cumulative return (alpha=0, beta=1)
    log_R0 = np.zeros(T_total + 1)
    log_R0[1:] = np.cumsum(np.log(G0_q))

    # Precompute A_{v,H}(alpha_k)^2 and A_{v,H}(beta_k)^2 for all v=0..V-1 and H in H_grid
    # Shape: A2_alpha[h_idx, k, v] and A2_beta[h_idx, k, v]
    v_origins = 4 * np.arange(V)
    nA = len(alpha_ann_grid)
    nB = len(beta_grid)
    nH = len(H_grid)

    A2_alpha = np.zeros((nH, nA, V))
    A2_beta = np.zeros((nH, nB, V))
    for hi, hq in enumerate(H_q):
        taus = np.arange(0, hq + 1)
        for v in range(V):
            t0 = v_origins[v]
            # R_{t0, t0+tau} = exp(log_R[:, t0+tau] - log_R[:, t0])
            idx_t = t0 + taus
            Ra = np.exp(log_R_a[:, idx_t] - log_R_a[:, t0 : t0 + 1])
            A2_alpha[hi, :, v] = (Ra.mean(axis=1)) ** 2
            Rb = np.exp(log_R_b[:, idx_t] - log_R_b[:, t0 : t0 + 1])
            A2_beta[hi, :, v] = (Rb.mean(axis=1)) ** 2

    # Precompute fixed deal entry/exit template per vintage (or random per rep—fixed deal counts, random timing!)
    # To isolate shock regimes cleanly with common random numbers, draw deal timings (a, h) once per rep
    # Actually: for ultra-fast vectorized computation across 600 reps, precompute per-vintage discount ratios!
    # Look at PV_i(theta) for fund i in vintage v:
    #   PV_i(theta) = (1/J) sum_{j=1}^J [ - D_{t0, t0+a_j}(theta) + D_{t0, t0+b_j}(theta) * R_{t0+a_j, t0+b_j}(theta_0) * N_{i,j} ]
    # Write PV_i(theta) = -C_in(theta) + sum_t D_{t0, t0+t}(theta) * X_{i, t},
    # where X_{i, t} is the aggregated exit cash flow of fund i at offset t in {1, ..., 66}!
    # OMG THAT MAKES EVALUATION ACROSS ALL 391 ALPHAS AND 501 BETAS A SINGLE MATRIX MULTIPLY!!

    max_off = max_deal_q  # 66
    # Precompute D_in_a[v, k, a] for a in 0..20 and D_ex_a[v, k, b] for b in 0..66
    # Wait: even simpler! Store net cash flows on offset grid tau_cf = 0..66 (67 points per fund!)
    # D_a[v, k, off] = exp(-(log_R_a[k, t0+off] - log_R_a[k, t0])) -> shape (V, nA, 67)
    offs = np.arange(0, max_off + 1)
    D_a = np.zeros((V, nA, max_off + 1))
    D_b = np.zeros((V, nB, max_off + 1))
    for v in range(V):
        t0 = v_origins[v]
        D_a[v] = np.exp(-(log_R_a[:, t0 + offs] - log_R_a[:, t0 : t0 + 1]))
        D_b[v] = np.exp(-(log_R_b[:, t0 + offs] - log_R_b[:, t0 : t0 + 1]))

    # Calibrate calendar-time shock e_u variance so that its within-vintage cross-fund PV variance
    # matches Regime 2 static vintage shock sigma_W^2 at theta_0!
    # Let's compute average overlap E[|[a1,b1] cap [a2,b2]|] across deals in the same vintage (in quarters)
    a_tmp = rng.integers(0, 21, size=4000)
    h_tmp = rng.integers(1, 47, size=4000)
    b_tmp = a_tmp + h_tmp
    ov_iid_q = np.mean(np.maximum(0, np.minimum(b_tmp[:2000], b_tmp[2000:]) - np.maximum(a_tmp[:2000], a_tmp[2000:])))
    # For quarterly i.i.d. e_u with variance s_e^2: cross-deal log-payoff cov = s_e^2 * ov_iid_q
    # Match cross-fund PV variance ~ sigma_W^2: set s_e_iid = sigma_W / sqrt(ov_iid_q)
    sig_e_iid = sigma_W / np.sqrt(ov_iid_q)

    # For quarterly AR(1) e_u with quarterly phi_e = 0.85 (annual 0.85^4 = 0.522):
    phi_e = 0.85
    # Compute average double sum sum_{u in [a1,b1], u' in [a2,b2]} phi_e^{|u-u'|}
    k_mat = phi_e ** np.abs(np.arange(max_off + 1)[:, None] - np.arange(max_off + 1)[None, :])
    ov_ar_list = []
    for m_i in range(2000):
        s1, e1 = a_tmp[m_i], b_tmp[m_i]
        s2, e2 = a_tmp[2000 + m_i], b_tmp[2000 + m_i]
        ov_ar_list.append(k_mat[s1:e1, s2:e2].sum())
    ov_ar_q = float(np.mean(ov_ar_list))
    sig_e_ar = sigma_W / np.sqrt(ov_ar_q)

    # Helper to optimize 1D grid with parabolic sub-grid refinement
    def argmin_refined(obj_vals, grid):
        idx = int(np.argmin(obj_vals))
        if idx == 0 or idx == len(grid) - 1:
            return float(grid[idx]), True
        # Parabolic interpolation on (idx-1, idx, idx+1)
        y0, y1, y2 = obj_vals[idx - 1], obj_vals[idx], obj_vals[idx + 1]
        denom = (y0 - 2.0 * y1 + y2)
        if denom > 1e-14:
            step = 0.5 * (y0 - y2) / denom
            step = np.clip(step, -0.5, 0.5)
            return float(grid[idx] + step * (grid[1] - grid[0])), False
        return float(grid[idx]), False

    # Regimes to simulate with common idiosyncratic draws!
    regime_names = [
        "R1_iid",          # c_v = 0
        "R2_static_r00",   # static W_v, rho_W = 0.0
        "R3a_AR1_r06",     # AR(1) W_v, rho_W = 0.6
        "R3a_AR1_r08",     # AR(1) W_v, rho_W = 0.8
        "R3b_cal_iid",     # calendar-time e_u, i.i.d.
        "R3b_cal_AR1",     # calendar-time e_u, AR(1) phi=0.85/qtr
    ]

    # Storage: estimates[regime][est_type][param] of shape (n_reps, nH)
    results = {
        rname: {
            "U_alpha": np.zeros((n_reps, nH)),
            "L2_alpha": np.zeros((n_reps, nH)),
            "U_beta": np.zeros((n_reps, nH)),
            "L2_beta": np.zeros((n_reps, nH)),
            "U_alpha_bnd": np.zeros((n_reps, nH), dtype=bool),
            "U_beta_bnd": np.zeros((n_reps, nH), dtype=bool),
            "Ucross_m0_beta": np.zeros(n_reps),
            "Ucross_m0_bnd": np.zeros(n_reps, dtype=bool),
            "Ucross_m5_beta": np.zeros(n_reps),
            "Ucross_m5_bnd": np.zeros(n_reps, dtype=bool),
            "Ucross_m7_beta": np.zeros(n_reps),
            "Ucross_m7_bnd": np.zeros(n_reps, dtype=bool),
        }
        for rname in regime_names
    }

    # Leave-band-out cross-vintage masks at H=0
    v_idx = np.arange(V)
    masks = {
        0: (np.abs(v_idx[:, None] - v_idx[None, :]) > 0).astype(float),
        5: (np.abs(v_idx[:, None] - v_idx[None, :]) > 5).astype(float),
        7: (np.abs(v_idx[:, None] - v_idx[None, :]) > 7).astype(float),
    }

    # Quarterly deal idiosyncratic vol
    sig_idio_q = sigma_idio_yr / 2.0

    for rep in range(n_reps):
        # Draw deal entry a in 0..20, hold h in 1..46 for all (V, n_funds, n_deals)
        a_d = rng.integers(0, 21, size=(V, n_funds, n_deals))
        h_d = rng.integers(1, 47, size=(V, n_funds, n_deals))
        b_d = a_d + h_d

        # True cumulative return for each deal: R_{t0+a, t0+b}(theta_0)
        t0_arr = v_origins[:, None, None]
        R_deal0 = np.exp(log_R0[t0_arr + b_d] - log_R0[t0_arr + a_d])

        # Idiosyncratic mean-one multiplicative deal shock: lognormal with var = sig_idio_q^2 * h
        var_idio = (sig_idio_q ** 2) * h_d
        z_idio = rng.normal(0.0, 1.0, size=(V, n_funds, n_deals))
        N_idio = np.exp(np.sqrt(var_idio) * z_idio - 0.5 * var_idio)

        # Draw common vintage shocks (AR(1) with rho in {0, 0.6, 0.8}) using common innovations
        eps_W = rng.normal(0.0, 1.0, size=V)
        W_dict = {}
        for rho in [0.0, 0.6, 0.8]:
            W = np.zeros(V)
            W[0] = sigma_W * eps_W[0]
            s_inn = sigma_W * np.sqrt(1.0 - rho ** 2)
            for v in range(1, V):
                W[v] = rho * W[v - 1] + s_inn * eps_W[v]
            W_dict[rho] = np.exp(W - 0.5 * sigma_W ** 2)  # mean-one multiplicative

        # Draw calendar-time PE shocks e_u (i.i.d. and AR(1))
        eps_e = rng.normal(0.0, 1.0, size=T_total + 1)
        # (a) i.i.d. calendar shock
        e_iid = sig_e_iid * eps_e
        cum_e_iid = np.zeros(T_total + 1)
        cum_e_iid[1:] = np.cumsum(e_iid[1:])
        sum_e_iid = cum_e_iid[t0_arr + b_d] - cum_e_iid[t0_arr + a_d]
        var_e_iid = (sig_e_iid ** 2) * h_d
        N_cal_iid = np.exp(sum_e_iid - 0.5 * var_e_iid)

        # (b) AR(1) calendar shock (phi_e = 0.85/qtr)
        e_ar = np.zeros(T_total + 1)
        e_ar[0] = sig_e_ar * eps_e[0]
        s_e_inn = sig_e_ar * np.sqrt(1.0 - phi_e ** 2)
        for t in range(1, T_total + 1):
            e_ar[t] = phi_e * e_ar[t - 1] + s_e_inn * eps_e[t]
        cum_e_ar = np.zeros(T_total + 1)
        cum_e_ar[1:] = np.cumsum(e_ar[1:])
        sum_e_ar = cum_e_ar[t0_arr + b_d] - cum_e_ar[t0_arr + a_d]
        # Exact variance of sum of h consecutive stationary AR(1) terms:
        # Var(sum_{u=1}^h e_u) = sig_e_ar^2 * [ h + 2 * sum_{k=1}^{h-1} (h-k) phi^k ]
        var_h_lookup = np.zeros(47)
        for hh in range(1, 47):
            ks = np.arange(1, hh)
            var_h_lookup[hh] = (sig_e_ar ** 2) * (hh + 2.0 * np.sum((hh - ks) * (phi_e ** ks)))
        var_e_ar = var_h_lookup[h_d]
        N_cal_ar = np.exp(sum_e_ar - 0.5 * var_e_ar)

        # Shared multiplier per regime
        shared_map = {
            "R1_iid": np.ones((V, 1, 1)),
            "R2_static_r00": W_dict[0.0][:, None, None],
            "R3a_AR1_r06": W_dict[0.6][:, None, None],
            "R3a_AR1_r08": W_dict[0.8][:, None, None],
            "R3b_cal_iid": N_cal_iid,
            "R3b_cal_AR1": N_cal_ar,
        }

        # Pre-build investment cash-flow histogram (same across all regimes!)
        # CF_inv[v, i, off] = -(1/n_deals) * count(a_d == off)
        CF_inv = np.zeros((V, n_funds, max_off + 1))
        for v in range(V):
            for i in range(n_funds):
                np.add.at(CF_inv[v, i], a_d[v, i], -1.0 / n_deals)

        for rname in regime_names:
            payoffs = (R_deal0 * N_idio * shared_map[rname]) / n_deals
            CF = CF_inv.copy()
            for v in range(V):
                for i in range(n_funds):
                    np.add.at(CF[v, i], b_d[v, i], payoffs[v, i])

            # --- ALPHA (beta=1 fixed) ---
            # PV_a[v, i, k] = sum_off CF[v, i, off] * D_a[v, k, off]
            PV_a = np.einsum("vio,vko->vik", CF, D_a)  # (V, n_funds, nA)
            PV_v_a = PV_a.mean(axis=1)                 # (V, nA)
            sq_v_a = PV_v_a ** 2                       # (V, nA)
            diag_a = (PV_a ** 2).mean(axis=1) / n_funds
            U0_a = (sq_v_a - diag_a) / (1.0 - 1.0 / n_funds)  # (V, nA)

            # Objective across all H: obj[hi, k] = mean_v(A2_alpha[hi, k, v] * U0_a[v, k])
            obj_U_a = np.einsum("hkv,vk->hk", A2_alpha, U0_a) / V
            obj_L2_a = np.einsum("hkv,vk->hk", A2_alpha, sq_v_a) / V

            for hi in range(nH):
                est, bnd = argmin_refined(obj_U_a[hi], alpha_ann_grid)
                results[rname]["U_alpha"][rep, hi] = est
                results[rname]["U_alpha_bnd"][rep, hi] = bnd
                est_l2, _ = argmin_refined(obj_L2_a[hi], alpha_ann_grid)
                results[rname]["L2_alpha"][rep, hi] = est_l2

            # --- BETA (alpha=0 fixed) ---
            PV_b = np.einsum("vio,vko->vik", CF, D_b)  # (V, n_funds, nB)
            PV_v_b = PV_b.mean(axis=1)                 # (V, nB)
            sq_v_b = PV_v_b ** 2
            diag_b = (PV_b ** 2).mean(axis=1) / n_funds
            U0_b = (sq_v_b - diag_b) / (1.0 - 1.0 / n_funds)

            obj_U_b = np.einsum("hkv,vk->hk", A2_beta, U0_b) / V
            obj_L2_b = np.einsum("hkv,vk->hk", A2_beta, sq_v_b) / V

            for hi in range(nH):
                est, bnd = argmin_refined(obj_U_b[hi], beta_grid)
                results[rname]["U_beta"][rep, hi] = est
                results[rname]["U_beta_bnd"][rep, hi] = bnd
                est_l2, _ = argmin_refined(obj_L2_b[hi], beta_grid)
                results[rname]["L2_beta"][rep, hi] = est_l2

            # --- Leave-band-out cross-vintage U_cross at H=0 (beta) ---
            # obj(k) = PV_v_b[:, k]^T M_m PV_v_b[:, k] / sum(M_m)
            for m_val in [0, 5, 7]:
                M_m = masks[m_val]
                obj_cr = np.einsum("vk,vw,wk->k", PV_v_b, M_m, PV_v_b) / M_m.sum()
                est_c, bnd_c = argmin_refined(obj_cr, beta_grid)
                results[rname][f"Ucross_m{m_val}_beta"][rep] = est_c
                results[rname][f"Ucross_m{m_val}_bnd"][rep] = bnd_c

    return H_grid, results


if __name__ == "__main__":
    print("=" * 78)
    print("PART 1: Exact Analytical Zero-Bias Horizon Table H_0*(T_eff, g, sigma_r)")
    print("=" * 78)
    T_vals = [3.0, 5.0, 7.0, 8.25, 10.0]
    g_vals = [0.0, 0.08, 0.10, 0.12, 0.133, 0.15, 0.20]
    sig_m = 0.154

    table1_rows = []
    header = "  T_eff | " + " | ".join(f"g={g*100:4.1f}%" for g in g_vals) + " | stoch(13.3%,15.4%)"
    print(header)
    for T in T_vals:
        vals = [solve_H_det(T, g) for g in g_vals]
        hs = solve_H_stoch(T, 0.133, sig_m)
        row_str = f"  {T:5.2f} | " + " | ".join(f"{v:6.2f}" for v in vals) + f" | {hs:6.2f}"
        print(row_str)
        table1_rows.append({"T": T, "det": dict(zip(g_vals, vals)), "stoch_133": hs})

    # Export curve CSV for Panel B of Figure
    T_curve = np.linspace(1.0, 12.0, 45)
    with open(os.path.join(DATA_DIR, "panelB_Hstar_curves.csv"), "w") as f:
        f.write("T_eff,H_g0,H_g133_det,H_g133_stoch,H_g20_det\n")
        for T in T_curve:
            f.write(
                f"{T:.4f},{2*T:.4f},{solve_H_det(T, 0.133):.4f},"
                f"{solve_H_stoch(T, 0.133, sig_m):.4f},{solve_H_det(T, 0.20):.4f}\n"
            )

    print("\n" + "=" * 78)
    print("PART 2: Effective Payoff Duration T_eff Across Shock Mechanisms")
    print("=" * 78)
    mech_core = compute_T_eff_mechanisms(I_max=5.0, h_min=0.25, h_max=11.5, g=0.133, sigma_m=0.154)
    mech_short = compute_T_eff_mechanisms(I_max=2.0, h_min=0.25, h_max=2.0, g=0.133, sigma_m=0.154)
    print("[Core calibration: I=5 yr, h in [0.25, 11.5] yr, g=13.3%, sigma_m=15.4%]")
    for r in mech_core:
        print(f"  {r['mechanism']:45s}: T_eff={r['T_eff']:5.2f} | 2T={r['H_flat']:5.2f} | det={r['H_det']:5.2f} | stoch={r['H_stoch']:5.2f}")
    print("[Short-life M5 2/2: I=2 yr, h in [0.25, 2.0] yr, g=13.3%, sigma_m=15.4%]")
    for r in mech_short:
        print(f"  {r['mechanism']:45s}: T_eff={r['T_eff']:5.2f} | 2T={r['H_flat']:5.2f} | det={r['H_det']:5.2f} | stoch={r['H_stoch']:5.2f}")

    print("\n" + "=" * 78)
    print("PART 3: Nonlinear Cash-Flow Monte Carlo (V=25, n_v=15, J=15, R=600)")
    print("=" * 78)
    H_grid, res = simulate_regimes_mc(V=25, n_funds=15, n_deals=15, n_reps=600, sigma_W=0.35)

    # Print summary for alpha (annualized %/yr and monthly %/mo) and beta
    print("\n--- U_{v,H} ALPHA Bias (%/yr) and RMSE (%/yr) [true alpha=0, beta=1 fixed] ---")
    key_H_idx = [0, 2, 4, 5, 6, 8, 10]  # 0, 5, 10, 12.5, 15, 20, 30
    for rname in res:
        arr = res[rname]["U_alpha"] * 100.0  # %/yr
        biases = arr.mean(axis=0)
        sds = arr.std(axis=0)
        rmses = np.sqrt(biases**2 + sds**2)
        h_rmse_opt = H_grid[int(np.argmin(rmses))]
        h_zb_idx = int(np.argmin(np.abs(biases)))
        print(f"\n[{rname}]  min-RMSE H*={h_rmse_opt:4.1f} yr (RMSE={rmses.min():.2f}%/yr) | zero-bias H~{H_grid[h_zb_idx]:4.1f} yr")
        print("  H (yr) : " + " ".join(f"{H_grid[i]:6.1f}" for i in key_H_idx))
        print("  Bias   : " + " ".join(f"{biases[i]:+6.2f}" for i in key_H_idx))
        print("  SD     : " + " ".join(f"{sds[i]:6.2f}" for i in key_H_idx))
        print("  RMSE   : " + " ".join(f"{rmses[i]:6.2f}" for i in key_H_idx))

    print("\n--- L2 VYP ALPHA Bias (%/yr) and RMSE (%/yr) in R1_iid ---")
    arr_l2 = res["R1_iid"]["L2_alpha"] * 100.0
    b_l2 = arr_l2.mean(axis=0)
    sd_l2 = arr_l2.std(axis=0)
    rmse_l2 = np.sqrt(b_l2**2 + sd_l2**2)
    print("  H (yr) : " + " ".join(f"{H_grid[i]:6.1f}" for i in key_H_idx))
    print("  Bias   : " + " ".join(f"{b_l2[i]:+6.2f}" for i in key_H_idx))
    print("  RMSE   : " + " ".join(f"{rmse_l2[i]:6.2f}" for i in key_H_idx))

    print("\n--- U_{v,H} BETA Bias and RMSE [true beta=1, alpha=0 fixed] ---")
    for rname in res:
        arr = res[rname]["U_beta"] - 1.0
        biases = arr.mean(axis=0)
        sds = arr.std(axis=0)
        rmses = np.sqrt(biases**2 + sds**2)
        h_rmse_opt = H_grid[int(np.argmin(rmses))]
        print(f"[{rname:15s}] min-RMSE H*={h_rmse_opt:4.1f} yr | "
              + " | ".join(f"H={H_grid[i]:.0f}: b={biases[i]:+.3f}, rmse={rmses[i]:.3f}" for i in [0, 4, 6, 10]))

    print("\n--- Leave-Band-Out Cross-Vintage U_cross_{0,m} Failure (beta, true=1) ---")
    for rname in ["R1_iid", "R2_static_r00", "R3a_AR1_r06", "R3b_cal_iid"]:
        line = f"  {rname:15s}: "
        for m in [0, 5, 7]:
            est = res[rname][f"Ucross_m{m}_beta"]
            bnd = res[rname][f"Ucross_m{m}_bnd"].mean() * 100.0
            rmse = np.sqrt(np.mean((est - 1.0) ** 2))
            line += f"m={m}: bnd={bnd:5.1f}%, rmse={rmse:5.3f} | "
        print(line)

    # Export Panel C CSV for TikZ/PGFPlots (alpha bias & RMSE across H for key regimes!)
    with open(os.path.join(DATA_DIR, "panelC_mc_regimes.csv"), "w") as f:
        cols = [
            "H",
            "R1_U_bias", "R1_U_rmse", "R1_L2_bias", "R1_L2_rmse",
            "R2_U_bias", "R2_U_rmse",
            "R3a06_U_bias", "R3a06_U_rmse",
            "R3a08_U_bias", "R3a08_U_rmse",
            "R3b_iid_U_bias", "R3b_iid_U_rmse",
            "R3b_ar_U_bias", "R3b_ar_U_rmse",
        ]
        f.write(",".join(cols) + "\n")
        for hi, H in enumerate(H_grid):
            def br(rname, key="U_alpha"):
                a = res[rname][key][:, hi] * 100.0
                b = float(a.mean())
                r = float(np.sqrt(np.mean(a ** 2)))
                return b, r
            r1u = br("R1_iid", "U_alpha")
            r1l = br("R1_iid", "L2_alpha")
            r2u = br("R2_static_r00", "U_alpha")
            r3a6 = br("R3a_AR1_r06", "U_alpha")
            r3a8 = br("R3a_AR1_r08", "U_alpha")
            r3bi = br("R3b_cal_iid", "U_alpha")
            r3ba = br("R3b_cal_AR1", "U_alpha")
            vals = [H, *r1u, *r1l, *r2u, *r3a6, *r3a8, *r3bi, *r3ba]
            f.write(",".join(f"{v:.4f}" for v in vals) + "\n")

    # Also save JSON summary for exact table generation
    summary = {
        "table1": table1_rows,
        "mech_core": mech_core,
        "mech_short": mech_short,
        "H_grid": H_grid.tolist(),
        "mc": {},
    }
    for rname in res:
        a = res[rname]["U_alpha"] * 100.0
        b = res[rname]["U_beta"] - 1.0
        summary["mc"][rname] = {
            "U_alpha_bias": a.mean(axis=0).tolist(),
            "U_alpha_sd": a.std(axis=0).tolist(),
            "U_alpha_rmse": np.sqrt(np.mean(a ** 2, axis=0)).tolist(),
            "U_beta_bias": b.mean(axis=0).tolist(),
            "U_beta_sd": b.std(axis=0).tolist(),
            "U_beta_rmse": np.sqrt(np.mean(b ** 2, axis=0)).tolist(),
            "Ucross_m0_bnd_pct": float(res[rname]["Ucross_m0_bnd"].mean() * 100.0),
            "Ucross_m0_rmse": float(np.sqrt(np.mean((res[rname]["Ucross_m0_beta"] - 1.0) ** 2))),
            "Ucross_m5_bnd_pct": float(res[rname]["Ucross_m5_bnd"].mean() * 100.0),
            "Ucross_m5_rmse": float(np.sqrt(np.mean((res[rname]["Ucross_m5_beta"] - 1.0) ** 2))),
            "Ucross_m7_bnd_pct": float(res[rname]["Ucross_m7_bnd"].mean() * 100.0),
            "Ucross_m7_rmse": float(np.sqrt(np.mean((res[rname]["Ucross_m7_beta"] - 1.0) ** 2))),
        }
    with open(os.path.join(DATA_DIR, "verification_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print("\nSaved data to:", DATA_DIR)
