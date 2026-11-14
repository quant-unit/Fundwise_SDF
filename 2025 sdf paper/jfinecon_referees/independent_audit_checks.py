"""Independent finite-support checks for the referee-review audit, 14 November 2026.

Run with Python and numpy. No source data are modified. Numerical examples are
counterexamples/algebra checks, not a reproduction of the paper's estimates.
"""
from itertools import product
import json
import numpy as np


def gradient(fn, theta, step=1e-6):
    theta = np.asarray(theta, dtype=float)
    return np.array([
        (fn(theta + step * e) - fn(theta - step * e)) / (2 * step)
        for e in np.eye(len(theta))
    ])


def exact_tilt():
    """Enumerate independent two-point deal shocks; check both alpha and beta."""
    factors = np.array([.05, -.02, .03, .06, -.01])
    theta0 = np.array([.002, 1.2])
    rf, nv = .01, 2
    starts, ends = np.array([0, 1, 1, 2]), np.array([3, 4, 5, 5])
    amplitudes = np.array([.1, .2, .3, .15])
    noise = 1 + np.array(list(product([-1, 1], repeat=4))) * amplitudes
    g0 = 1 + rf + theta0[0] + theta0[1] * factors
    d0 = np.r_[1., np.cumprod(g0)]
    proceeds = d0[ends] / d0[starts]
    cf = np.zeros((len(noise), 6))
    for j, (a, b) in enumerate(zip(starts, ends)):
        cf[:, a] -= 1 / nv
        cf[:, b] += proceeds[j] * noise[:, j] / nv
    s = np.vstack([np.zeros(2), np.cumsum(
        np.column_stack([np.ones(5), factors]) / g0[:, None], axis=0)])
    out = []
    for horizon in [0, 2, 5]:
        def criterion(theta):
            d = np.r_[1., np.cumprod(1 + rf + theta[0] + theta[1] * factors)]
            return np.mean((np.sum(cf / d, axis=1) * d[:horizon + 1].mean()) ** 2)
        a = d0[:horizon + 1].mean()
        sbar = (d0[:horizon + 1, None] * s[:horizon + 1]).sum(0) / d0[:horizon + 1].sum()
        analytic = 2 * a**2 / nv**2 * np.sum(
            (amplitudes**2 / d0[starts]**2)[:, None] * (sbar - s[ends]), axis=0)
        numeric = gradient(criterion, theta0)
        assert np.max(np.abs(analytic - numeric)) < 1e-7
        out.append(dict(horizon=horizon, analytic=analytic.tolist(), numeric=numeric.tolist()))
    return out


def sign_counterexample():
    """Positive E[S] does not imply positive E[S Var(N|G)]."""
    f = np.array([.2, -.1])
    g = 1 + f
    sigma = .2
    s, nu = f / g, sigma**2 / g**2
    score = -2 * np.mean(s * nu)
    assert np.mean(s) > 0 and score > 0
    return dict(mean_S=float(s.mean()), mean_S_nu=float(np.mean(s*nu)),
                beta_score_at_truth=float(score), growth_optimal_beta=2.5)


def heterogeneous_pair():
    """Two noise-free, correctly benchmarked funds with different schedules."""
    rows = []
    for beta in [.5, 1., 1.5]:
        # Fund 1: -1 at 0, +1.2 at 1. Fund 2: -1 at 1, +.8 at 2.
        g1, g2 = 1 + .2 * beta, 1 - .2 * beta
        p1, p2 = -1 + 1.2 / g1, (-1 + .8 / g2) / g1
        rows.append(dict(beta=beta, pv1=p1, pv2=p2, pair_product=p1*p2))
    assert rows[0]['pair_product'] < 0 and rows[2]['pair_product'] < 0
    assert abs(rows[1]['pair_product']) < 1e-20
    return rows


def positive_product_identity():
    """The proposed correction works for conditional iid errors."""
    rows = []
    for mu in [-.3, 0., .4]:
        samples = mu + np.array(list(product([-.2, .2], repeat=3)))
        pair = (samples.sum(1)**2 - (samples**2).sum(1)) / 6
        corrected = samples.mean(1)**2 - samples.var(1, ddof=1)/3
        assert np.max(abs(pair-corrected)) < 1e-15
        assert abs(pair.mean() - mu**2) < 1e-15
        rows.append(dict(mu=mu, expected_pair=float(pair.mean()), squared_mean=mu**2))
    return rows


def horizon_counterexamples():
    # At H=2, for an exit at 1 and true gross return 1,
    # Q(r)=[(r+1+1/r)^2/9]*[(r-1)^2+variance], globally minimized at r=1.
    grid = np.linspace(.5, 1.5, 10001)
    variance = .04
    q = ((grid + 1 + 1/grid)**2/9)*((grid-1)**2 + variance)
    assert grid[q.argmin()] == 1
    # A negative factor realization can make fitted beta RISE with H.
    beta = np.linspace(.5, 1.5, 100001)
    g0, g = .9, 1 - .1*beta
    loss0 = (g0/g - 1)**2 + .01*(g0/g)**2
    loss1 = loss0 * ((1+g)/2)**2
    b0, b1 = beta[loss0.argmin()], beta[loss1.argmin()]
    assert b1 > b0
    return dict(positive_noise_consistent_horizon_gross_return=float(grid[q.argmin()]),
                negative_factor_beta_h0=float(b0), negative_factor_beta_h1=float(b1))


def simulation_timing():
    gross = np.array([1.1, 1.2, 1.3])
    cumulative = gross.cumprod()
    # Investment at array position 0 and exit at position 2.
    inclusive_payoff = gross.prod()
    consistent_payoff = gross[1:].prod()
    def inception_npv(payoff):
        return -1 + payoff * cumulative[0] / cumulative[2]
    assert abs(inception_npv(inclusive_payoff) - .1) < 1e-14
    assert abs(inception_npv(consistent_payoff)) < 1e-14
    return dict(inclusive_start_payoff=float(inclusive_payoff),
                inception_npv_inclusive=float(inception_npv(inclusive_payoff)),
                inception_npv_start_plus_one=float(inception_npv(consistent_payoff)))


def clt_centering_counterexample():
    # Independent bounded errors, deterministic drift a_i=i^(-1/4).
    # The limiting criterion is theta^2+1, but sqrt(n)*E[theta_hat] diverges.
    out = []
    for n in [100, 10000, 1000000]:
        drift = np.arange(1, n+1, dtype=float)**(-.25)
        out.append(dict(n=n, population_center=float(drift.mean()),
                        root_n_center=float(np.sqrt(n)*drift.mean())))
    return out


def sdf_pricing_and_recovery():
    """Check Appendix A's restricted SDF and bounded-pool recovery examples."""
    payoff = np.array([1.2, .8])
    benchmark = 1 / payoff
    common_sdf = np.ones(2)
    assert abs(np.mean(benchmark * payoff) - 1) < 1e-14
    assert abs(benchmark.mean() - 25/24) < 1e-14
    assert abs(np.mean((common_sdf - benchmark) * payoff)) < 1e-14
    assert abs(np.mean(common_sdf * payoff) - 1) < 1e-14
    assert common_sdf.mean() == 1

    # Exhaust all independent two-point noise realizations, not just a formula.
    r0, amplitude = 1.1, .3
    pools = []
    for n in [1, 2, 4, 8]:
        noise = 1 + amplitude * np.array(list(product([-1, 1], repeat=n)))
        pooled_payoffs = r0 * noise.mean(axis=1)
        variance = amplitude**2 / n
        r_star = r0 * (1 + variance)
        def loss(r):
            return np.mean((pooled_payoffs / r - 1)**2)
        step = 1e-5
        slope = (loss(r_star+step) - loss(r_star-step)) / (2*step)
        mean_error = np.mean(pooled_payoffs/r_star - 1)
        assert abs(slope) < 1e-8
        assert loss(r_star) < loss(r0)
        assert abs(mean_error + variance/(1+variance)) < 1e-14
        pools.append(dict(funds=n, fitted_gross_return=r_star,
                          mean_pricing_error=float(mean_error)))
    return dict(pe_price=float(np.mean(benchmark*payoff)),
                benchmark_bond_price=float(benchmark.mean()),
                common_sdf_bond_price=float(common_sdf.mean()), pools=pools)


if __name__ == '__main__':
    print(json.dumps(dict(
        exact_tilt=exact_tilt(),
        conditional_iid_product=positive_product_identity(),
        sign_counterexample=sign_counterexample(),
        heterogeneous_pair=heterogeneous_pair(),
        horizon_counterexamples=horizon_counterexamples(),
        simulation_timing=simulation_timing(),
        clt_centering_counterexample=clt_centering_counterexample(),
        sdf_pricing_and_recovery=sdf_pricing_and_recovery(),
    ), indent=2))
