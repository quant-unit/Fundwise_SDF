"""Independent, deterministic checks for revision_plan_audit.md.

Run from any directory with Python, numpy and pandas. Reads the existing
global-data probe module; does not modify data, the plan, or estimation code.
"""

import importlib.util
import sys
from math import erf, exp, pi, sqrt
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True


def main():
    # T5: an infinite-alpha limit can beat truth but lose to a finite choice.
    def loss(gross):
        return 1 - 2 / gross + 10 / gross**2

    assert loss(1) == 9
    assert np.isclose(loss(10), 0.9)
    assert loss(10) < 1 < loss(1)
    print("T5: Q(truth)=9, Q(finite gross=10)=0.9, Q(infinity)=1.")

    # Fixed-V infill retains realized common shocks.
    shocks = np.array([0.2, 0.4])
    realized_argmin = 1 - shocks.mean()
    assert np.isclose(realized_argmin, 0.7)
    print("Fixed-V infill: realized minimizer=0.7, expected-loss minimizer=1.")

    # Leaving out v does not remove its covariance with the training sample.
    v, count, rho = 12, 25, 0.8
    covariance = sum(rho ** abs(v - w) for w in range(count) if w != v) / (count - 1)
    assert covariance > 0
    print(f"AR(1) leave-one-out training/residual covariance={covariance:.6f}.")

    # Negative market return reverses the asserted beta monotonicity.
    def horizon_weight(beta):
        return (1 + (1 - 0.1 * beta)) / 2

    derivative = (horizon_weight(1.0001) - horizon_weight(0.9999)) / 0.0002
    assert np.isclose(derivative, -0.05)
    print("One-period negative factor: dA/dbeta=-0.05.")

    # Exact one-period Gaussian truncation expectation, not Monte Carlo.
    gross, sigma = 1.0, 2.0
    x = gross / sigma
    truncated_mean = gross * (1 + erf(x / sqrt(2))) / 2 + sigma * exp(-x*x / 2) / sqrt(2*pi)
    assert truncated_mean > gross
    print(f"Gaussian default truncation: E[max(1+shock,0)]={truncated_mean:.6f}>1.")

    source = Path(__file__).with_name("check_revision_probes.py")
    spec = importlib.util.spec_from_file_location("revision_probes", source)
    probes = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(probes)
    q, flows = probes.load()

    # The manuscript's VC LS columns use n>=2, unlike the cited script's default.
    sample = probes.build(flows, q, "VC", "sample")
    predetermined = probes.build(flows, q, "VC", "pre")
    expected_sample = {0: 1.93, 5: 2.20, 10: 1.78, 15: 1.32, 20: 1.24, 30: 0.90}
    expected_pre = {0: 1.97, 5: 2.16, 10: 1.80, 15: 1.34, 20: 1.23, 30: -0.48}
    print("VC squared loss: H | script all-vintage | n>=2 | plan | pre all | pre n>=2 | plan")
    for h in expected_sample:
        values = [
            probes.argmin(lambda b: probes.crit(panel, b, 12*h, "LS", minn=minimum))
            for panel, minimum in [(sample, 1), (sample, 2), (predetermined, 1), (predetermined, 2)]
        ]
        assert np.isclose(values[1], expected_sample[h])
        assert np.isclose(values[3], expected_pre[h])
        print(f"{h:2d} | {values[0]:.2f} | {values[1]:.2f} | {expected_sample[h]:.2f}"
              f" | {values[2]:.2f} | {values[3]:.2f} | {expected_pre[h]:.2f}")

    buyout = probes.build(flows, q, "BO", "sample")
    print("Buyout residual/within-vintage SE at beta=1.21:")
    for year in [1996, 2006, 2009, 2010]:
        rf, market, cf = buyout[year]
        errors = (cf / np.cumprod(1 + rf + 1.21*market)[None, :]).sum(1)
        se = errors.std(ddof=1) / sqrt(len(errors))
        ratio = errors.mean() / se
        assert abs(ratio) < 2
        print(f"{year}: PV={errors.mean():+.6f}, SE={se:.6f}, ratio={ratio:+.6f}")

    for asset in ["BO", "VC"]:
        panel = probes.build(flows, q, asset, "pre")
        zero_origins = sum(np.all(cf[:, 0] == 0) for _, _, cf in panel.values())
        print(f"{asset}: zero initial flow at predetermined origin in {zero_origins}/{len(panel)} vintages.")

    print("All deterministic assertions passed.")


if __name__ == "__main__":
    main()
