# Numerical checks for Appendix A/B of review_of_referee_reports.tex (horizon theory).
# Stand-alone replica; NOT the paper's pipeline. Simulated i.i.d. market, one-factor benchmark, alpha = 0.
#   D1  closed-form tilt of the squared-loss criterion at beta0 for every H  vs  numerical derivative
#   D2  least-squares plateau by horizon (large V)                            -> crosses beta0 where the tilt changes sign
#   D3  cross-fund product criterion (no squared idiosyncratic noise)         -> beta0 at every H, bounded funds per vintage
import numpy as np, sys

rf, mu, sd, beta0 = 0.003, 0.006, 0.045, 1.0
L = 181
HS = (0, 5, 10, 15, 20)
GRID = np.arange(0.0, 3.2001, 0.01)


def gen(V, nv, sig, seed):
    """Fund-level cash flows CF[v, i, k]; funds of a vintage share the factor path. 15 deals per fund."""
    rng = np.random.default_rng(seed)
    mkt = rng.normal(mu, sd, size=V * 12 + L + 12)
    idx = (12 * np.arange(V))[:, None] + np.arange(1, L)[None, :]          # months 1..180 after vintage date
    CF = np.zeros((V, nv, L)); k = np.arange(1, L)[None, :]
    W = np.zeros(V); WS = np.zeros(V); T0 = np.zeros(V, dtype=int)
    for v in range(V):
        F = mkt[idx[v]]; G = 1 + rf + beta0 * F
        nd = nv * 15
        start = rng.integers(1, 61, size=nd); end = np.minimum(start + rng.integers(1, 121, size=nd), L - 1)
        inwin = (k > start[:, None]) & (k <= end[:, None])
        M = G[None, :] + rng.normal(0, sig, size=(nd, L - 1))
        mult = np.where(inwin, M, 1.0).prod(axis=1); mult[((M < 0) & inwin).any(axis=1)] = 0.0
        fund = np.repeat(np.arange(nv), 15)
        np.add.at(CF[v], (fund, start), -1.0); np.add.at(CF[v], (fund, end), mult)
        # ingredients of the closed-form tilt (deal level, conditional on the path)
        t0 = start.min(); T0[v] = t0
        cumlogG = np.concatenate([[0.0], np.cumsum(np.log(G))])            # log R^b(0, k)
        cumS = np.concatenate([[0.0], np.cumsum(F / G)])                   # S(0, k)
        cumV = np.concatenate([[0.0], np.cumsum(np.log1p(sig**2 / G**2))])
        c = np.exp(-(cumlogG[start] - cumlogG[t0]))                        # 1 / R^b(t0, start_j)
        vj = np.exp(cumV[end] - cumV[start]) - 1.0                         # Var(N_j | path)
        Sb = cumS[end] - cumS[t0]                                          # S(t0, exit_j)
        W[v] = (c**2 * vj).sum() / nv**2; WS[v] = (c**2 * vj * Sb).sum() / nv**2
    return CF, mkt, idx, W, WS, T0


def pieces(CF, mkt, idx, beta):
    """Fund-level inception NPVs (inception of the vintage pool) and the horizon weights A_H."""
    V, nv, _ = CF.shape; r = np.arange(V)
    D = np.ones((V, L)); D[:, 1:] = np.cumprod(1 + rf + beta * mkt[idx], axis=1)
    pool = CF.sum(axis=1); nz = pool != 0
    t0 = nz.argmax(axis=1); tmax = L - 1 - nz[:, ::-1].argmax(axis=1)
    pv = (CF / D[:, None, :]).sum(axis=2) * D[r, t0][:, None]              # (V, nv)
    cs = np.concatenate([np.zeros((V, 1)), np.cumsum(D, axis=1)], axis=1)
    A = {}
    for H in HS:
        hi = np.minimum(t0 + 12 * H, tmax)
        A[H] = (cs[r, hi + 1] - cs[r, t0]) / (hi - t0 + 1) / D[r, t0]
    return pv, A, t0, tmax


def crit(pv, A):
    """Squared-loss criterion and cross-fund product criterion, summed over vintages."""
    nv = pv.shape[1]; e = pv * A[:, None]
    m = e.mean(axis=1); ls = (m**2).sum()
    u = ((e.sum(axis=1)**2 - (e**2).sum(axis=1)) / (nv * (nv - 1))).sum() if nv > 1 else np.nan
    return ls, u


def plateau(V, nv, sig, seed):
    CF, mkt, idx, *_ = gen(V, nv, sig, seed)
    ls = {H: [] for H in HS}; us = {H: [] for H in HS}
    for b in GRID:
        pv, A, _, _ = pieces(CF, mkt, idx, b)
        for H in HS:
            l, u = crit(pv, A[H]); ls[H].append(l); us[H].append(u)
    out_ls = [GRID[int(np.argmin(ls[H]))] for H in HS]
    out_u = [GRID[int(np.argmin(us[H]))] if nv > 1 else np.nan for H in HS]
    return out_ls, out_u


def small_sample(Vbig, nv, sig, seed, block=20):
    """Mean and s.d. of the two estimators over blocks of 20 consecutive vintages."""
    CF, mkt, idx, *_ = gen(Vbig, nv, sig, seed)
    nb = Vbig // block
    ls = {H: np.zeros((nb, len(GRID))) for H in HS}; us = {H: np.zeros((nb, len(GRID))) for H in HS}
    for g, b in enumerate(GRID):
        pv, A, _, _ = pieces(CF, mkt, idx, b)
        for H in HS:
            e = pv * A[H][:, None]; m = e.mean(axis=1)
            l = m**2; u = (e.sum(axis=1)**2 - (e**2).sum(axis=1)) / (nv * (nv - 1))
            ls[H][:, g] = l[:nb * block].reshape(nb, block).sum(axis=1)
            us[H][:, g] = u[:nb * block].reshape(nb, block).sum(axis=1)
    res = {}
    for H in HS:
        bl = GRID[ls[H].argmin(axis=1)]; bu = GRID[us[H].argmin(axis=1)]
        res[H] = (bl.mean(), bl.std(), np.median(bl), bu.mean(), bu.std(), np.median(bu),
                  np.mean((bu <= GRID[0] + 1e-9) | (bu >= GRID[-1] - 1e-9)))
    return res


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"

    if what in ("D1", "all"):
        print("=" * 100)
        print("D1  Tilt of the squared-loss criterion at beta0: closed form vs numerical derivative (n_v=1, sigma=0.1)")
        for seed in (21, 22):
            V = 40000
            CF, mkt, idx, W, WS, T0 = gen(V, 1, 0.1, seed)
            h = 0.005
            pvp, Ap, _, _ = pieces(CF, mkt, idx, beta0 + h); pvm, Am, _, _ = pieces(CF, mkt, idx, beta0 - h)
            pv0, A0, t0, tmax = pieces(CF, mkt, idx, beta0)
            F = mkt[idx]; G = 1 + rf + beta0 * F
            D = np.ones((V, L)); D[:, 1:] = np.cumprod(G, axis=1)
            cumS = np.concatenate([np.zeros((V, 1)), np.cumsum(F / G, axis=1)], axis=1)
            r = np.arange(V)
            for H in HS:
                dv = ((pvp[:, 0] * Ap[H])**2 - (pvm[:, 0] * Am[H])**2) / (2 * h); num = dv.mean(); nse = dv.std() / np.sqrt(V)
                hi = np.minimum(t0 + 12 * H, tmax)
                mask = (np.arange(L)[None, :] >= t0[:, None]) & (np.arange(L)[None, :] <= hi[:, None])
                Rb = D / D[r, t0][:, None]; St = cumS - cumS[r, t0][:, None]
                Sbar = (Rb * St * mask).sum(axis=1) / (Rb * mask).sum(axis=1)
                terms = 2 * A0[H]**2 * (W * Sbar - WS)
                print(f"    seed {seed}  H={H:2d}: closed form {terms.mean():+.4f} (s.e. {terms.std() / np.sqrt(V):.4f})   numerical {num:+.4f} (s.e. {nse:.4f})")

    if what in ("D2", "all"):
        print("=" * 100)
        print("D2/D3  Large-V plateaus (10,000 vintages). LS = squared loss, U = cross-fund product criterion. H =", HS)
        for sig, nv, seed in ((0.1, 1, 3), (0.1, 2, 3), (0.1, 5, 3), (0.1, 5, 4), (0.2, 5, 3), (0.2, 5, 4), (0.2, 20, 3)):
            l, u = plateau(10000, nv, sig, seed)
            print(f"    sigma={sig} n_v={nv:2d} seed={seed}:  LS {np.round(l, 2)}   U {np.round(u, 2)}", flush=True)

    if what in ("D4", "all"):
        print("=" * 100)
        print("D4  Samples of 20 vintages (500 blocks): mean / s.d. / median of beta-hat; share of U estimates on the grid boundary")
        for sig, nv, seed in ((0.1, 5, 7), (0.2, 5, 7), (0.2, 20, 7)):
            res = small_sample(10000, nv, sig, seed)
            for H in HS:
                a = res[H]
                print(f"    sigma={sig} n_v={nv:2d} H={H:2d}:  LS {a[0]:.2f} / {a[1]:.2f} / {a[2]:.2f}    U {a[3]:.2f} / {a[4]:.2f} / {a[5]:.2f}   boundary {a[6]:.2f}", flush=True)


def pieces_ab(CF, mkt, idx, alpha, beta):
    V, nv, _ = CF.shape; r = np.arange(V)
    D = np.ones((V, L)); D[:, 1:] = np.cumprod(1 + alpha + rf + beta * mkt[idx], axis=1)
    pool = CF.sum(axis=1); nz = pool != 0
    t0 = nz.argmax(axis=1); tmax = L - 1 - nz[:, ::-1].argmax(axis=1)
    pv = (CF / D[:, None, :]).sum(axis=2) * D[r, t0][:, None]
    cs = np.concatenate([np.zeros((V, 1)), np.cumsum(D, axis=1)], axis=1)
    A = {}
    for H in (0, 15):
        hi = np.minimum(t0 + 12 * H, tmax)
        A[H] = (cs[r, hi + 1] - cs[r, t0]) / (hi - t0 + 1) / D[r, t0]
    return pv, A


if __name__ == "__main__" and (len(sys.argv) > 1 and sys.argv[1] == "D5"):
    print("=" * 100)
    print("D5  Two parameters (monthly alpha, beta), truth (0, 1). Large-V plateaus on a grid; LS vs cross-fund product criterion")
    AG = np.arange(-0.010, 0.0301, 0.001); BG = np.arange(0.0, 3.0001, 0.05)
    for sig, nv, seed, V in ((0.1, 5, 3, 4000), (0.2, 5, 3, 4000)):
        CF, mkt, idx, *_ = gen(V, nv, sig, seed)
        res = {(k, H): np.zeros((len(AG), len(BG))) for k in ("LS", "U") for H in (0, 15)}
        for ia, a in enumerate(AG):
            for ib, b in enumerate(BG):
                pv, A = pieces_ab(CF, mkt, idx, a, b)
                for H in (0, 15):
                    l, u = crit(pv, A[H]); res[("LS", H)][ia, ib] = l; res[("U", H)][ia, ib] = u
        for H in (0, 15):
            out = []
            for k in ("LS", "U"):
                ia, ib = np.unravel_index(res[(k, H)].argmin(), res[(k, H)].shape)
                out.append(f"{k}: alpha {100 * AG[ia]:+.1f}% p.m., beta {BG[ib]:.2f}")
            print(f"    sigma={sig} n_v={nv} H={H:2d}:  " + "   ".join(out), flush=True)
