# Check of Referee 2 comment 3.1 / Referee 1 comment 2:
# does the H=0 (and H>0) quadratic criterion target beta0 when n_v is bounded and V grows?
import numpy as np, sys, time
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
V = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
L = 181            # offsets 0..180 months (5y invest + 10y hold)
rf, mu, sd, beta0, sig = 0.003, 0.006, 0.045, 1.0, (float(sys.argv[3]) if len(sys.argv) > 3 else 0.2)
mkt = rng.normal(mu, sd, size=V * 12 + L + 12)
Hs = [0, 5, 10, 15, 20, 30]
betas = np.arange(0.0, 3.0001, 0.01)

def sim_vyp(v, nv):
    s = 12 * v
    G = 1 + rf + beta0 * mkt[s + 1: s + L]          # gross benchmark return, offsets 1..180
    nd = nv * 15
    start = rng.integers(1, 61, size=nd)            # invest offset 1..60
    hold = rng.integers(1, 121, size=nd)
    end = np.minimum(start + hold, L - 1)
    k = np.arange(1, L)[None, :]
    inwin = (k > start[:, None]) & (k <= end[:, None])
    M = G[None, :] + rng.normal(0, sig, size=(nd, L - 1))
    dflt = ((M < 0) & inwin).any(axis=1)
    mult = np.where(inwin, M, 1.0).prod(axis=1)
    mult[dflt] = 0.0
    cf = np.zeros(L)
    np.add.at(cf, start, -1.0)
    np.add.at(cf, end, mult)
    return cf / nv

def eps_bar(CF, beta):
    # CF: [V, L]; returns dict H -> eps_bar vector [V]
    idx = (12 * np.arange(V))[:, None] + np.arange(1, L)[None, :]
    D = np.ones((V, L)); D[:, 1:] = np.cumprod(1 + rf + beta * mkt[idx], axis=1)
    nz = CF != 0
    t0 = nz.argmax(axis=1); tmax = L - 1 - nz[:, ::-1].argmax(axis=1)
    pv0 = (CF / D).sum(axis=1)                      # value at offset 0
    cs = np.concatenate([np.zeros((V, 1)), np.cumsum(D, axis=1)], axis=1)
    out = {}
    r = np.arange(V)
    for H in Hs:
        hi = np.minimum(t0 + 12 * H, tmax)
        A = (cs[r, hi + 1] - cs[r, t0]) / (hi - t0 + 1)   # mean compounding factor over tau in T_{i,H}
        out[H] = pv0 * A
    return out

for nv in [1, 5, 20, 80]:
    t = time.time()
    CF = np.vstack([sim_vyp(v, nv) for v in range(V)])
    E = {H: np.empty((len(betas), V)) for H in Hs}
    for b, beta in enumerate(betas):
        e = eps_bar(CF, beta)
        for H in Hs: E[H][b] = e[H]
    line = f"n_v={nv:3d} | plateau (V={V}): "
    for H in Hs:
        line += f"H={H}: {betas[(E[H]**2).sum(axis=1).argmin()]:.2f}  "
    print(line)
    nb = V // 20
    line = f"        | mean of beta-hat over {nb} blocks of V=20: "
    for H in Hs:
        sq = (E[H]**2)[:, :nb * 20].reshape(len(betas), nb, 20).sum(axis=2)
        bh = betas[sq.argmin(axis=0)]
        line += f"H={H}: {bh.mean():.2f} (sd {bh.std():.2f})  "
    print(line, f"[{time.time()-t:.0f}s]", flush=True)
