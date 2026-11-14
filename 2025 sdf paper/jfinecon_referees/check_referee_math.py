# Numerical checks for the mathematical audit of the referee reports (JFEC-2026-141).
# Stand-alone replica; NOT the paper's pipeline. Simulated i.i.d. market, one-factor benchmark.
# Each block prints the numbers quoted in review_of_referee_reports.tex (Section "Mathematical audit").
import numpy as np
from math import erf, sqrt, exp

rf, mu, sd, beta0 = 0.003, 0.006, 0.045, 1.0       # monthly risk-free, factor mean, factor sd, true beta
L, BASE = 181, 120                                   # fund offsets 0..180; BASE months of market before vintage 0
Phi = lambda x: 0.5 * (1 + erf(x / sqrt(2)))


def gen(V, nv, sig, seed):
    """Vintage portfolios as in the paper's baseline: 15 deals per fund, entry in years 1-5, hold 1-10 years."""
    rng = np.random.default_rng(seed)
    mkt = rng.normal(mu, sd, size=BASE + V * 12 + L + 12)
    CF = np.zeros((V, L)); k = np.arange(1, L)[None, :]
    for v in range(V):
        s = BASE + 12 * v
        G = 1 + rf + beta0 * mkt[s + 1: s + L]
        nd = nv * 15
        start = rng.integers(1, 61, size=nd); end = np.minimum(start + rng.integers(1, 121, size=nd), L - 1)
        inwin = (k > start[:, None]) & (k <= end[:, None])
        M = G[None, :] + rng.normal(0, sig, size=(nd, L - 1))
        mult = np.where(inwin, M, 1.0).prod(axis=1); mult[((M < 0) & inwin).any(axis=1)] = 0.0
        np.add.at(CF[v], start, -1.0); np.add.at(CF[v], end, mult)
    return CF / nv, mkt


def Dmat(mkt, V, beta, alpha=0.0, lo=-60):
    """D[v, j] = benchmark gross return from offset 0 to offset offs[j] (offs = lo..L-1; <1 for negative offsets)."""
    offs = np.arange(lo, L)
    idx = (BASE + 12 * np.arange(V))[:, None] + offs[None, :]
    c = np.cumsum(np.log(1 + alpha + rf + beta * mkt[idx]), axis=1)
    return np.exp(c - c[:, [-lo]])


def unit_dates(CF):
    nz = CF != 0
    return nz.argmax(axis=1), L - 1 - nz[:, ::-1].argmax(axis=1)


def eps_bar(CF, mkt, beta, H, alpha=0.0, lo=-60):
    """Horizon-averaged pricing error. H in years; H<0 averages over the |H| years BEFORE inception."""
    V = CF.shape[0]; D = Dmat(mkt, V, beta, alpha, lo); t0, tmax = unit_dates(CF); r = np.arange(V)
    pv0 = (CF / D[:, -lo:]).sum(axis=1)                              # value at offset 0
    cs = np.concatenate([np.zeros((V, 1)), np.cumsum(D, axis=1)], axis=1)
    a = t0 + min(12 * H, 0) - lo; b = np.minimum(t0 + max(12 * H, 0), tmax) - lo
    A = (cs[r, b + 1] - cs[r, a]) / (b - a + 1)
    return pv0 * A, pv0 * D[r, t0 - lo]                              # (eps_bar, inception PV)


def ls_beta(CF, mkt, H, grid=np.arange(0.0, 3.2001, 0.01)):
    obj = [np.sum(eps_bar(CF, mkt, b, H)[0] ** 2) for b in grid]
    return grid[int(np.argmin(obj))]


print("=" * 100)
print("C1  Identity eps_tau = PV * R^b(t0,tau), for a pre- and a post-inception date (brute force vs identity)")
CF1, mkt1 = gen(10000, 1, 0.1, 11)
t0s, _ = unit_dates(CF1)
for v, beta in ((7, 1.3), (123, 0.6)):
    s = BASE + 12 * v; t0 = t0s[v]
    Rb = lambda a, b: np.prod(1 + rf + beta * mkt1[s + a + 1: s + b + 1]) if b >= a else 1 / np.prod(1 + rf + beta * mkt1[s + b + 1: s + a + 1])
    pv = sum(CF1[v, k] / Rb(t0, k) for k in range(L) if CF1[v, k] != 0)
    for tau in (t0 - 24, t0 + 60):
        brute = sum(CF1[v, k] * Rb(k, tau) for k in range(L) if CF1[v, k] != 0)
        print(f"    unit {v}, beta={beta}, tau=t0{tau - t0:+d}: brute {brute:+.10f}  identity {pv * Rb(t0, tau):+.10f}  gap {abs(brute - pv * Rb(t0, tau)):.1e}")

print("=" * 100)
print("C2  Wedge under the DLP-type DGP: sample mean of eps_bar(theta0) over 10,000 vintages (n_v=1, sigma=0.1)")
for H in (-5, 0, 5, 15):
    e = eps_bar(CF1, mkt1, beta0, H)[0]
    print(f"    H={H:+3d}: mean {e.mean():+.4f}   s.e. {e.std() / sqrt(len(e)):.4f}   t = {e.mean() / (e.std() / sqrt(len(e))):+.2f}")

print("=" * 100)
print("C3  Wedge under a genuine SDF (lognormal market, exp-affine kernel that prices market and T-bill)")
rng = np.random.default_rng(5)
sig_m = 0.045; lrf = np.log(1 + rf); mu_m = 0.0055                     # log-return parameters, monthly
lam = (mu_m + sig_m ** 2 / 2 - lrf) / sig_m                             # price of risk so that E[M R_m] = 1
n, h = 400000, 60
z = rng.standard_normal((n, h))
logM = -lrf - lam * z - lam ** 2 / 2; logR = mu_m + sig_m * z
cM, cR = np.cumsum(logM, axis=1), np.cumsum(logR, axis=1)
delta = np.exp(cM[:, -1] + cR[:, -1]) - 1                               # deal = market held h months, Inv = -1
print(f"    lambda (monthly) {lam:.4f}; check E[M R_m]-1 = {np.exp(logM[:, 0] + logR[:, 0]).mean() - 1:+.5f}, E[M R_f]-1 = {np.exp(logM[:, 0] + lrf).mean() - 1:+.5f}")
print(f"    inception moment E[delta] = {delta.mean():+.4f} (s.e. {delta.std() / sqrt(n):.4f})")
for hp in (12, 36, 60):
    w = np.exp(-cM[:, hp - 1]) * delta
    cf = exp((lrf + sig_m * lam) * hp) - exp((lrf + lam ** 2) * hp)
    print(f"    valuation {hp:2d} months after investment: MC E[delta/Psi] = {w.mean():+.4f} (s.e. {w.std() / sqrt(n):.4f})   closed form {cf:+.4f}")
del z, logM, logR, cM, cR

print("=" * 100)
print("C4  Gradient of E[PV^2] at beta0 (single deals, 60-month hold, sigma=0.05) and 1/n scaling of the shift")
rng = np.random.default_rng(6); sig = 0.05; T = 60; hstep = 1e-3
grid = np.arange(0.90, 1.4001, 0.005)
for npool in (1, 4, 16):
    U = 200000 // npool * 2; chunks = 20; Uc = U // chunks
    acc = dict(cp=0.0, cm=0.0, Sv=0.0, NS=0.0, S2w=0.0, curv2=0.0); obj = np.zeros(len(grid))
    for _ in range(chunks):
        F = rng.normal(mu, sd, size=(Uc, T)); G = 1 + rf + beta0 * F
        Nbar = np.mean([np.prod(1 + rng.normal(0, sig, size=(Uc, T)) / G, axis=1) for _ in range(npool)], axis=0)
        S = (F / G).sum(axis=1); Q = (F ** 2 / G ** 2).sum(axis=1); v = np.prod(1 + sig ** 2 / G ** 2, axis=1) - 1
        LG0 = np.log(G).sum(axis=1)
        rho = lambda b: np.exp(LG0 - np.log(1 + rf + b * F).sum(axis=1))
        acc["cp"] += np.sum((Nbar * rho(beta0 + hstep) - 1) ** 2); acc["cm"] += np.sum((Nbar * rho(beta0 - hstep) - 1) ** 2)
        acc["Sv"] += np.sum(S * v); acc["NS"] += np.sum((Nbar - 1) * Nbar * S)
        acc["S2w"] += np.sum(S ** 2 * (1 + v / npool)); acc["curv2"] += np.sum((S ** 2 + Q) * v) / npool
        for i, b in enumerate(grid): obj[i] += np.sum((Nbar * rho(b) - 1) ** 2)
    fd = (acc["cp"] - acc["cm"]) / (2 * hstep) / U
    g_formula = -2 * acc["Sv"] / U / npool; g_sample = -2 * acc["NS"] / U
    newton = -g_formula / (2 * acc["S2w"] / U + 2 * acc["curv2"] / U)
    i = int(np.argmin(obj)); p = np.polyfit(grid[i - 2:i + 3], obj[i - 2:i + 3], 2); bstar = -p[1] / (2 * p[0])
    print(f"    n={npool:2d}: gradient finite-diff {fd:+.5f} | -2E[(N-1)N S] {g_sample:+.5f} | -2E[S Var(N|G)]/n {g_formula:+.5f}"
          f" || minimiser {bstar:.4f}  shift {bstar - 1:+.4f}  first-order prediction {newton:+.4f}  shift*n {npool * (bstar - 1):.3f}")

print("=" * 100)
print("C5  Pre-inception averaging (H<0) versus post-inception: least-squares plateau, n_v=1, sigma=0.1, V=10,000")
print("    " + "   ".join(f"H={H:+d}: {ls_beta(CF1, mkt1, H):.2f}" for H in (-5, 0, 5, 15)))

print("=" * 100)
print("C6  Noise of realised-SDF discounting versus benchmark deflation (lognormal, 5-year payoff)")
rng = np.random.default_rng(7); n = 400000; hh = 60; beta, sig_e = 1.2, 0.08
gam = lam / sig_m
rm = mu_m * hh + sig_m * sqrt(hh) * rng.standard_normal(n); e = sig_e * sqrt(hh) * rng.standard_normal(n)
log_sdf = -gam * rm + beta * rm + e                                     # log(M X) up to a constant
log_defl = -beta * rm + beta * rm + e                                   # log(X / R^b) up to a constant
print(f"    gamma = lambda/sigma_m = {gam:.2f}, beta = {beta}")
print(f"    Var log(M X):   MC {log_sdf.var():.4f}   formula (beta-gamma)^2 sigma_m^2 h + sigma_e^2 h = {(beta - gam) ** 2 * sig_m ** 2 * hh + sig_e ** 2 * hh:.4f}")
print(f"    Var log(X/R^b): MC {log_defl.var():.4f}   formula sigma_e^2 h = {sig_e ** 2 * hh:.4f}")

print("=" * 100)
print("C7  Realised-kernel sample means (Martin 2012 mechanism), normalised kernel K_h = M_h exp(rf h), E[K_h] = 1")
rng = np.random.default_rng(8)
for years in (5, 10, 15):
    hh = 12 * years; l2h = lam ** 2 * hh
    K = np.exp(-lam * sqrt(hh) * rng.standard_normal((20000, 30)) - l2h / 2)
    sm = K.mean(axis=1)
    print(f"    h={years:2d}y: lambda^2 h = {l2h:.2f} | median K {exp(-l2h / 2):.3f} | P(K<1) {Phi(sqrt(l2h) / 2):.3f} | Var K {exp(l2h) - 1:8.1f}"
          f" | 30-draw sample mean: median {np.median(sm):.2f}, P(mean<0.8) {np.mean(sm < 0.8):.2f}")
print("    Comparable variance terms exp(x)-1 for the deflator route, 10-year horizon:")
for label, x in (("beta gap 0.3: 0.3^2 sigma_m^2 h", 0.09 * sig_m ** 2 * 120), ("weight A_H at H=15: beta^2 sigma_m^2 180", sig_m ** 2 * 180),
                 ("one deal, idiosyncratic sigma=0.2/month", 0.04 * 120), ("300 pooled deals, same sigma (approx.)", None)):
    val = (exp(0.04 * 120) - 1) / 300 if x is None else exp(x) - 1
    print(f"      {label:45s} {val:8.3f}")

print("=" * 100)
print("C8  Alpha-beta separation: linearised information in samples of 20 vintages (n_v=20, sigma=0.2)")
CF2, mkt2 = gen(4000, 20, 0.2, 12)
pv = eps_bar(CF2, mkt2, beta0, 0)[1]
ga = (eps_bar(CF2, mkt2, beta0, 0, alpha=1e-5)[1] - eps_bar(CF2, mkt2, beta0, 0, alpha=-1e-5)[1]) / 2e-5
gb = (eps_bar(CF2, mkt2, beta0 + 1e-4, 0)[1] - eps_bar(CF2, mkt2, beta0 - 1e-4, 0)[1]) / 2e-4
B = 200; ga, gb, pvb = ga.reshape(B, 20), gb.reshape(B, 20), pv.reshape(B, 20)
rho = (ga * gb).sum(1) / np.sqrt((ga ** 2).sum(1) * (gb ** 2).sum(1)); vif = 1 / (1 - rho ** 2)
su = pv.std()
se_b = su / np.sqrt((gb ** 2).sum(1)); se_a = su / np.sqrt((ga ** 2).sum(1))
print(f"    sd of PV(theta0) per vintage {su:.3f}; cosine between alpha- and beta-gradients: median {np.median(rho):.3f} (5%-95%: {np.quantile(rho, .05):.3f} to {np.quantile(rho, .95):.3f})")
print(f"    variance inflation 1/(1-rho^2): median {np.median(vif):.1f} (5%-95%: {np.quantile(vif, .05):.1f} to {np.quantile(vif, .95):.1f})")
print(f"    linearised s.e. of beta: alone {np.median(se_b):.3f}  -> with alpha {np.median(se_b * np.sqrt(vif)):.3f}")
print(f"    linearised s.e. of alpha (% per month): beta known {100 * np.median(se_a):.3f} -> with beta {100 * np.median(se_a * np.sqrt(vif)):.3f}")

print("=" * 100)
print("C9  Terminal-value haircut and the one-factor beta (same sample; payoffs after year 10 scaled by c)")
late = np.arange(L)[None, :] > 120
for c in (1.0, 0.75, 0.5, 0.0):
    CFc = np.where(late & (CF2 > 0), c * CF2, CF2)
    print(f"    c = {c:.2f}: least-squares beta at H=0 = {ls_beta(CFc, mkt2, 0):.2f}")

print("=" * 100)
print("C4b Gradient formula with a larger sample (4 million single deals per seed)")
for seed in (1, 2):
    rng = np.random.default_rng(seed); a = b = 0.0; xs = []; n = 0; sig = 0.05; T = 60
    for _ in range(40):
        F = rng.normal(mu, sd, (100000, T)); G = 1 + rf + beta0 * F
        N = np.prod(1 + rng.normal(0, sig, (100000, T)) / G, axis=1)
        S = (F / G).sum(1); v = np.prod(1 + sig ** 2 / G ** 2, axis=1) - 1
        x = (N - 1) * N * S; xs.append(x); a += x.sum(); b += (S * v).sum(); n += 100000
    xs = np.concatenate(xs)
    print(f"    seed {seed}: -2E[(N-1)N S] = {-2 * a / n:.5f} (s.e. {2 * xs.std() / sqrt(n):.5f})   -2E[S Var(N|G)] = {-2 * b / n:.5f}")

print("=" * 100)
print("C10 Finite-sample distribution of the mean inception NPV at theta0 (cents per dollar invested)")
for nv, sig, seed in ((20, 0.2, 21), (20, 0.1, 22)):
    CF, mkt = gen(4000, nv, sig, seed)
    pvc = eps_bar(CF, mkt, beta0, 0)[1] / 15 * 100          # 15 units invested per fund
    for V in (20, 50, 200):
        m = pvc[: (4000 // V) * V].reshape(-1, V).mean(axis=1)
        print(f"    n_v={nv}, sigma={sig}, V={V:3d}: mean {m.mean():+6.2f}  median {np.median(m):+6.2f}  sd {m.std():6.2f}  share below zero {np.mean(m < 0):.2f}")
