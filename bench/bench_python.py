"""ssfortran against statsmodels, both from Python (times in ms, minimum
over repeats). Run from the repo root after building the library:

    .venv/bin/python bench/bench_python.py
"""

import sys
import time
import warnings
from pathlib import Path

import numpy as np
import statsmodels.api as sm

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))
import ssfortran as ss  # noqa: E402

warnings.simplefilter("ignore")
rng = np.random.default_rng(1)


def best(f, reps=5):
    b = np.inf
    for _ in range(reps):
        t0 = time.perf_counter()
        f()
        b = min(b, time.perf_counter() - t0)
    return 1000 * b


def rw_noise(n):
    return np.cumsum(np.sqrt(0.1) * rng.standard_normal(n)) + rng.standard_normal(n)


def row(name, ours, theirs):
    print(f"{name:<44}{ours:12.2f}{theirs:14.2f}{theirs / ours:10.1f}x")


class LocalLevel(ss.MLEModel):
    def __init__(self, y):
        super().__init__(y, k_states=1, k_params=2)
        self["design"] = [[1.0]]
        self["transition"] = [[1.0]]
        self["selection"] = [[1.0]]
        self.ssm.initialize_diffuse()
        self._v = np.var(y)

    @property
    def start_params(self):
        return np.array([self._v / 2] * 2)

    def transform_params(self, x):
        return np.exp(2 * np.asarray(x))

    def untransform_params(self, p):
        return 0.5 * np.log(p)

    def update(self, params):
        self["obs_cov"] = [[params[0]]]
        self["state_cov"] = [[params[1]]]


if __name__ == "__main__":
    print(f"{'':<44}{'ssfortran':>12}{'statsmodels':>14}{'ratio':>11}")
    y = rw_noise(100_000)
    ours = ss.StructuralModel(y, [ss.Irregular(), ss.Level()])
    theirs = sm.tsa.UnobservedComponents(y, "llevel", use_exact_diffuse=True)
    p = [1.0, 0.1]
    row("loglike, local level, n = 1e5", best(lambda: ours.loglike(p)),
        best(lambda: theirs.loglike(p)))
    row("smooth, local level, n = 1e5", best(lambda: ours.smooth(p)),
        best(lambda: theirs.smooth(p)))

    y = rw_noise(1000)
    ours = ss.StructuralModel(y, [ss.Irregular(), ss.Level(), ss.Seasonal(12, "trig")])
    theirs = sm.tsa.UnobservedComponents(y, "llevel", freq_seasonal=[{"period": 12}],
                                         use_exact_diffuse=True)
    row("fit, level + trig seasonal, n = 1000", best(lambda: ours.fit(), 3),
        best(lambda: theirs.fit(disp=False), 3))

    ys = [rw_noise(100) for _ in range(1000)]
    t_ours = best(lambda: ss.fit_many([ss.StructuralModel(v, [ss.Irregular(), ss.Level()])
                                       for v in ys]), 3)
    t_serial = best(lambda: [ss.StructuralModel(v, [ss.Irregular(), ss.Level()]).fit()
                             for v in ys], 1)
    t_sm = best(lambda: [sm.tsa.UnobservedComponents(v, "llevel", use_exact_diffuse=True)
                         .fit(disp=False) for v in ys], 1)
    row("1000 local level fits, n = 100 (fit_many)", t_ours, t_sm)
    row("1000 local level fits, n = 100 (serial)", t_serial, t_sm)

    y = rw_noise(1000)
    row("fit, local level as a Python MLEModel", best(lambda: LocalLevel(y).fit(), 3),
        best(lambda: sm.tsa.UnobservedComponents(y, "llevel", use_exact_diffuse=True)
             .fit(disp=False), 3))
