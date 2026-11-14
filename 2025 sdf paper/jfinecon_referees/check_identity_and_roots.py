# Stress test of three claims: (1) identity eps_tau = PV * compounding factor,
# (2) H=0 least-squares plateau != beta0 for bounded n_v, (3) mean-moment estimator has no such tilt.
import numpy as np, sys
V, L = 10000, 181
rf, mu, sd, beta0 = 0.003, 0.006, 0.045, 1.0
betas = np.arange(0.0, 3.2001, 0.01)

def run(seed, sig, nv):
    rng = np.random.default_rng(seed)
    mkt = rng.normal(mu, sd, size=V * 12 + L + 12)
    idx = (12 * np.arange(V))[:, None] + np.arange(1, L)[None, :]
    CF = np.zeros((V, L))
    for v in range(V):
        G = 1 + rf + beta0 * mkt[idx[v]]
        nd = nv * 15
        start = rng.integers(1, 61, size=nd); end = np.minimum(start + rng.integers(1, 121, size=nd), L - 1)
        k = np.arange(1, L)[None, :]
        inwin = (k > start[:, None]) & (k <= end[:, None])
        M = G[None, :] + rng.normal(0, sig, size=(nd, L - 1))
        mult = np.where(inwin, M, 1.0).prod(axis=1); mult[((M < 0) & inwin).any(axis=1)] = 0.0
        np.add.at(CF[v], start, -1.0); np.add.at(CF[v], end, mult)
    CF /= nv
    nz = CF != 0; t0 = nz.argmax(axis=1); tmax = L - 1 - nz[:, ::-1].argmax(axis=1); r = np.arange(V)
    ls0, ls15, mm = [], [], []
    for b in betas:
        D = np.ones((V, L)); D[:, 1:] = np.cumprod(1 + rf + b * mkt[idx], axis=1)
        pv = (CF / D).sum(axis=1) * D[r, t0]                 # inception NPV
        cs = np.concatenate([np.zeros((V, 1)), np.cumsum(D, axis=1)], axis=1)
        hi = np.minimum(t0 + 180, tmax)
        A = (cs[r, hi + 1] - cs[r, t0]) / (hi - t0 + 1) / D[r, t0]
        if abs(b - 1.3) < 1e-9:                              # identity check: brute-force average of date-tau errors
            v = 7; taus = np.arange(t0[v], hi[v] + 1)
            brute = np.mean([(CF[v] * D[v, tau] / D[v]).sum() for tau in taus])
            ident = abs(brute - pv[v] * A[v])
        ls0.append((pv**2).sum()); ls15.append(((pv * A)**2).sum()); mm.append(pv.mean())
    mm = np.array(mm); roots = betas[:-1][np.sign(mm[:-1]) != np.sign(mm[1:])]
    print(f"sig={sig} n_v={nv} seed={seed}: LS H=0 {betas[np.argmin(ls0)]:.2f} | LS H=15 {betas[np.argmin(ls15)]:.2f} | mean-moment roots at H=0 {np.round(roots,2)} | identity gap {ident:.1e}", flush=True)

for seed in (3, 4, 5):
    run(seed, 0.1, 1)
run(3, 0.1, 5); run(4, 0.1, 5)
run(3, 0.2, 1); run(4, 0.2, 1)
