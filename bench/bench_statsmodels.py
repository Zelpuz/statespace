"""statsmodels side of the timing benchmark (see example/bench.f90).

Same models and parameters; each case is the minimum over repeats, in ms.
`loglike` uses statsmodels' memory-conserving mode, like our `loglike`.
Run from the repo root:

    .venv/bin/python bench/bench_statsmodels.py
"""

import time
import warnings

import numpy as np
import statsmodels.api as sm
from statsmodels.tsa.statespace.kalman_filter import (
    FILTER_CONVENTIONAL,
    FILTER_UNIVARIATE,
    MEMORY_CONSERVE,
)
from statsmodels.tsa.statespace.kalman_smoother import KalmanSmoother

warnings.simplefilter("ignore")
rng = np.random.default_rng(1)


def rw_noise(n):
    return np.cumsum(np.sqrt(0.1) * rng.standard_normal(n)) + rng.standard_normal(n)


def best(f, reps):
    b = np.inf
    for _ in range(reps):
        t0 = time.perf_counter()
        f()
        b = min(b, time.perf_counter() - t0)
    return 1000 * b


def nreps(n):
    return max(3, min(20, 2_000_000 // n))


def time_ssm(name, ks, n):
    """ks: a KalmanSmoother with data and matrices set."""
    r = nreps(n)

    def ll():
        ks.set_conserve_memory(MEMORY_CONSERVE)
        ks.loglike()
        ks.set_conserve_memory(0)

    t = [best(ll, r), best(ks.filter, r), best(ks.smooth, r)]
    print(f"{name:>8}{n:>10}" + "".join(f"{x:14.3f}" for x in t))


def ll_case(n):
    ks = KalmanSmoother(1, 1, 1)
    ks.bind(rw_noise(n)[None, :])
    ks["design"] = [[1.0]]
    ks["transition"] = [[1.0]]
    ks["selection"] = [[1.0]]
    ks["obs_cov"] = [[1.0]]
    ks["state_cov"] = [[0.1]]
    ks.initialize_known(np.zeros(1), np.array([[1e6]]))
    time_ssm("ll", ks, n)


def model_case(name, mod, params, n):
    mod.update(params)
    time_ssm(name, mod.ssm, n)


def bsm_case(n):
    mod = sm.tsa.UnobservedComponents(
        rw_noise(n),
        level="llevel",
        freq_seasonal=[{"period": 12}],
        use_exact_diffuse=True,
    )
    model_case("bsm", mod, [1.0, 0.1, 0.01], n)


def mv_case(n):
    tau = np.array([3, 6, 12, 24, 36, 60, 84, 120.0])
    lam = 0.06
    z = np.exp(-lam * tau)
    Z = np.column_stack([np.ones(8), (1 - z) / (lam * tau), (1 - z) / (lam * tau) - z])
    ks = KalmanSmoother(8, 3, 3)
    ks.bind(np.asfortranarray(rng.standard_normal((8, n))))
    ks["design"] = Z
    ks["obs_cov"] = 0.01 * np.eye(8)
    ks["transition"] = np.diag([0.99, 0.95, 0.9])
    ks["selection"] = np.eye(3)
    ks["state_cov"] = 0.1 * np.eye(3)
    ks.initialize_known(np.zeros(3), np.eye(3))
    time_ssm("mv", ks, n)
    ks.filter_method = FILTER_UNIVARIATE
    time_ssm("mv-uv", ks, n)
    ks.filter_method = FILTER_CONVENTIONAL
    y = rng.standard_normal((8, n))
    y[rng.random((8, n)) < 0.1] = np.nan
    ks.bind(np.asfortranarray(y))
    time_ssm("mv-miss", ks, n)


def arma_case(n):
    mod = sm.tsa.SARIMAX(rng.standard_normal(n), order=(2, 0, 1))
    model_case("arma", mod, [0.5, 0.2, 0.3, 1.0], n)


def many_case(nseries, n):
    ys = [rw_noise(n) for _ in range(nseries)]
    t0 = time.perf_counter()
    for y in ys:
        sm.tsa.UnobservedComponents(y, level="llevel", use_exact_diffuse=True).fit(
            disp=False
        )
    print(
        f"\nfit {nseries} local level series of n = {n}: "
        f"{1000 * (time.perf_counter() - t0):10.1f} ms"
    )


if __name__ == "__main__":
    print(f"{'case':>8}{'n':>10}{'loglike':>14}{'filter':>14}{'filter+smooth':>14}")
    for n in (10_000, 100_000, 1_000_000):
        ll_case(n)
    for n in (10_000, 100_000):
        bsm_case(n)
    for n in (10_000, 100_000):
        mv_case(n)
    for n in (10_000, 100_000):
        arma_case(n)
    many_case(1000, 100)
