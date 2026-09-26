"""Generate reference results from statsmodels for the Fortran test suite.

Run from the repo root:  .venv/bin/python test/fixtures/make_fixtures.py

Each fixture is a plain-text file of named arrays. Every array is written as a
header line ``name ndim d1 [d2 [d3]]`` followed by its values in Fortran
(column-major) order, one per line.
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.tsa.statespace.kalman_smoother import KalmanSmoother

HERE = "test/fixtures"


def write_fixture(name, arrays):
    with open(f"{HERE}/{name}.txt", "w") as f:
        for key, val in arrays.items():
            val = np.atleast_1d(np.asarray(val, dtype=float))
            f.write(f"{key} {val.ndim} {' '.join(map(str, val.shape))}\n")
            for x in val.ravel(order="F"):
                f.write(f"{float(x)!r}\n")


def results_arrays(res):
    return {
        "a": res.predicted_state,
        "yhat": res.forecasts,
        "P": res.predicted_state_cov,
        "v": res.forecasts_error,
        "F": res.forecasts_error_cov,
        "att": res.filtered_state,
        "Ptt": res.filtered_state_cov,
        "K": res.kalman_gain,
        "alphahat": res.smoothed_state,
        "V": res.smoothed_state_cov,
        "epshat": res.smoothed_measurement_disturbance,
        "epsvar": res.smoothed_measurement_disturbance_cov,
        "etahat": res.smoothed_state_disturbance,
        "etavar": res.smoothed_state_disturbance_cov,
        "llf_obs": res.llf_obs,
    }


# Initialization kinds, matching INIT_* in src/statespace_rep.f90
INIT_KNOWN, INIT_APPROX_DIFFUSE, INIT_STATIONARY, INIT_DIFFUSE = 1, 2, 3, 4


def diffuse_arrays(res):
    """Exact diffuse output. In diffuse periods statsmodels reports v, F, Finf
    and the eps variances element by element in LDL-transformed coordinates."""
    return {
        "nobs_diffuse": res.nobs_diffuse,
        "Pinf": res.predicted_diffuse_state_cov,
        "Finf": res.forecasts_error_diffuse_cov,
    }


def ssm_inputs(ssm):
    """System matrices of a statsmodels representation, time-varying ones
    with a trailing time axis."""
    return dict(Z=ssm["design"], H=ssm["obs_cov"], T=ssm["transition"],
                R=ssm["selection"], Q=ssm["state_cov"], c=ssm["state_intercept"],
                d=ssm["obs_intercept"])


def generic(name, y, Z, H, T, R, Q, c, d, a1=None, P1=None):
    """Low-level representation; known initialization if a1 and P1 are given,
    otherwise stationary."""
    p, m, r = Z.shape[0], T.shape[0], Q.shape[0]
    ks = KalmanSmoother(k_endog=p, k_states=m, k_posdef=r)
    ks.bind(np.asfortranarray(y))
    ks["design"], ks["obs_cov"], ks["transition"] = Z, H, T
    ks["selection"], ks["state_cov"] = R, Q
    ks["state_intercept"], ks["obs_intercept"] = c, d
    if a1 is None:
        ks.initialize_stationary()
        inputs = dict(y=y, Z=Z, H=H, T=T, R=R, Q=Q, c=c, d=d)
    else:
        ks.initialize_known(a1, P1)
        inputs = dict(y=y, Z=Z, H=H, T=T, R=R, Q=Q, c=c, d=d, a1=a1, P1=P1)
    res = ks.smooth()
    write_fixture(name, {**inputs, **results_arrays(res), "llf": res.llf_obs.sum()})
    print(name, "llf =", res.llf_obs.sum())


def spd(rng, k, scale=1.0):
    A = rng.standard_normal((k, k))
    return scale * (A @ A.T + k * np.eye(k))


def nile():
    y = pd.read_csv("data/nile.csv")["volume"].to_numpy(dtype=float)
    params = [15099.0, 1469.1]  # sigma2_eps, sigma2_eta (DK section 2.2.4)
    one = np.ones((1, 1))
    generic("nile_llevel_known", y[None, :], one, params[0] * one, one, one,
            params[1] * one, np.zeros(1), np.zeros(1), np.zeros(1), 1e7 * one)

    # DK 2.7: observations 21-40 and 61-80 treated as missing
    ymiss = y.copy()
    ymiss[20:40] = np.nan
    ymiss[60:80] = np.nan
    generic("nile_llevel_missing", ymiss[None, :], one, params[0] * one, one, one,
            params[1] * one, np.zeros(1), np.zeros(1), np.zeros(1), 1e7 * one)

    # MLE with statsmodels' default approximate diffuse init (kappa = 1e6,
    # loglikelihood_burn = 1), and with exact diffuse init as in DK.
    for tag, exact in [("approx", False), ("exact", True)]:
        mod = sm.tsa.UnobservedComponents(y, "llevel", use_exact_diffuse=exact)
        fit = mod.fit(disp=False, cov_type="approx")
        # The likelihood is flat near the optimum, so statsmodels' default
        # tolerances stop early. Also record a tightly converged optimum.
        tight = mod.fit(disp=False, pgtol=1e-12, factr=10, maxiter=5000)
        tight = mod.filter(tight.params, cov_type="approx")
        write_fixture(f"nile_llevel_mle_{tag}", {
            "start_params": mod.start_params, "params": fit.params, "llf": fit.llf,
            "params_tight": tight.params, "llf_tight": tight.llf,
            "bse_tight": tight.bse, "cov_params_tight": tight.cov_params(),
            "aic_tight": tight.aic, "bic_tight": tight.bic})
        print(f"nile_llevel_mle_{tag}", fit.params, fit.llf)


def multivariate(time_varying):
    rng = np.random.default_rng(20260925)
    p, m, r, n = 2, 3, 2, 50
    Z = rng.standard_normal((p, m))
    H = spd(rng, p, 0.5)
    T = 0.3 * rng.standard_normal((m, m)) + 0.5 * np.eye(m)
    R = rng.standard_normal((m, r))
    Q = spd(rng, r, 0.2)
    c = 0.1 * rng.standard_normal(m)
    d = rng.standard_normal(p)
    if time_varying:
        Z = Z[:, :, None] + 0.2 * rng.standard_normal((p, m, n))
        T = T[:, :, None] + 0.05 * rng.standard_normal((m, m, n))
        H = np.stack([spd(rng, p, 0.5) for _ in range(n)], axis=-1)
        d = d[:, None] + 0.1 * rng.standard_normal((p, n))
    a1 = rng.standard_normal(m)
    P1 = spd(rng, m, 2.0)
    y = rng.standard_normal((p, n)) * 3
    name = "mv_timevarying" if time_varying else "mv_invariant"
    generic(name, y, Z, H, T, R, Q, c, d, a1, P1)


def multivariate_missing():
    """Correlated H, partially and fully missing y_t, and 8 trailing missing
    periods whose one-step forecasts equal the out-of-sample forecasts."""
    rng = np.random.default_rng(7)
    p, m, r, n = 3, 2, 2, 40
    Z = rng.standard_normal((p, m))
    H = spd(rng, p, 0.5)
    T = np.array([[0.8, 0.1], [-0.2, 0.6]])
    R = np.eye(m)
    Q = spd(rng, r, 0.3)
    c = np.array([0.1, -0.2])
    d = rng.standard_normal(p)
    y = rng.standard_normal((p, n)) * 2
    y[0, 3] = np.nan
    y[[0, 2], 10] = np.nan
    y[:, 15] = np.nan
    y[1, 20:25] = np.nan
    y[:, n - 8:] = np.nan
    generic("mv_missing", y, Z, H, T, R, Q, c, d, rng.standard_normal(m), spd(rng, m))


def diffuse_models():
    y = pd.read_csv("data/nile.csv")["volume"].to_numpy(dtype=float)
    one = np.ones((1, 1))

    # Local level and local linear trend on Nile, exact diffuse (DK 5.2)
    for name, Z, T, Q in [
            ("nile_llevel_exact", one, one, 1469.1 * one),
            ("nile_lltrend_exact", np.array([[1.0, 0.0]]), np.array([[1.0, 1.0], [0.0, 1.0]]),
             np.diag([1500.0, 10.0]))]:
        m = T.shape[0]
        ks = KalmanSmoother(k_endog=1, k_states=m, k_posdef=m)
        ks.bind(y[None, :].copy(order="F"))
        ks["design"], ks["obs_cov"], ks["transition"] = Z, 15099.0 * one, T
        ks["selection"], ks["state_cov"] = np.eye(m), Q
        ks.initialize_diffuse()
        res = ks.smooth()
        inputs = dict(y=y[None, :], Z=Z, H=15099.0 * one, T=T, R=np.eye(m), Q=Q,
                      c=np.zeros(m), d=np.zeros(1), init_blocks=[[1, m, INIT_DIFFUSE]])
        write_fixture(name, {**inputs, **results_arrays(res), **diffuse_arrays(res),
                             "llf": res.llf_obs.sum()})
        print(name, "nobs_diffuse", res.nobs_diffuse, "llf", res.llf_obs.sum())

    # Correlated H with several elements per period (see test notes on coordinates)
    rng = np.random.default_rng(20260925)
    p, m, r, n = 2, 3, 2, 50
    Z = rng.standard_normal((p, m))
    H = spd(rng, p, 0.5)
    T = 0.3 * rng.standard_normal((m, m)) + 0.5 * np.eye(m)
    R = rng.standard_normal((m, r))
    Q = spd(rng, r, 0.2)
    c = 0.1 * rng.standard_normal(m)
    d = rng.standard_normal(p)
    ymv = rng.standard_normal((p, n)) * 3
    ymv[1, 0] = np.nan          # partly missing inside the diffuse period
    # ... and the same with Z, T, H and d time-varying
    Ztv = Z[:, :, None] + 0.2 * rng.standard_normal((p, m, n))
    Ttv = T[:, :, None] + 0.05 * rng.standard_normal((m, m, n))
    Htv = np.stack([spd(rng, p, 0.5) for _ in range(n)], axis=-1)
    dtv = d[:, None] + 0.1 * rng.standard_normal((p, n))
    for name, (Zx, Hx, Tx, dx) in [("mv_diffuse", (Z, H, T, d)),
                                   ("mv_diffuse_timevarying", (Ztv, Htv, Ttv, dtv))]:
        ks = KalmanSmoother(k_endog=p, k_states=m, k_posdef=r)
        ks.bind(np.asfortranarray(ymv))
        ks["design"], ks["obs_cov"], ks["transition"] = Zx, Hx, Tx
        ks["selection"], ks["state_cov"] = R, Q
        ks["state_intercept"], ks["obs_intercept"] = c, dx
        ks.initialize_diffuse()
        res = ks.smooth()
        inputs = dict(y=ymv, Z=Zx, H=Hx, T=Tx, R=R, Q=Q, c=c, d=dx,
                      init_blocks=[[1, m, INIT_DIFFUSE]])
        write_fixture(name, {**inputs, **results_arrays(res), **diffuse_arrays(res),
                             "llf": res.llf_obs.sum()})
        print(name, "nobs_diffuse", res.nobs_diffuse, "llf", res.llf_obs.sum())

    # Trend + seasonal with missing data, and level + stationary AR(1) (mixed init)
    ymiss = y.copy()
    ymiss[[2, 3, 30, 31, 32, 60]] = np.nan
    for name, kwargs, params, blocks in [
            ("uc_trend_seasonal_exact", dict(level="lltrend", seasonal=4),
             [15000.0, 1500.0, 10.0, 200.0], [[1, 5, INIT_DIFFUSE]]),
            ("uc_level_ar1_mixed", dict(level="llevel", autoregressive=1),
             [15000.0, 1000.0, 500.0, 0.5], [[1, 1, INIT_DIFFUSE], [2, 2, INIT_STATIONARY]])]:
        mod = sm.tsa.UnobservedComponents(ymiss, use_exact_diffuse=True, **kwargs)
        mod.update(params)
        res = mod.ssm.smooth()
        inputs = dict(y=ymiss[None, :], **ssm_inputs(mod.ssm), init_blocks=blocks)
        write_fixture(name, {**inputs, **results_arrays(res), **diffuse_arrays(res),
                             "llf": res.llf_obs.sum()})
        print(name, "nobs_diffuse", res.nobs_diffuse, "llf", res.llf_obs.sum())


def multivariate_stationary():
    rng = np.random.default_rng(11)
    p, m, r, n = 2, 3, 2, 60
    Z = rng.standard_normal((p, m))
    H = spd(rng, p, 0.5)
    T = np.array([[0.5, 0.3, 0.0], [0.2, 0.4, 0.1], [0.0, -0.3, 0.7]])
    R = rng.standard_normal((m, r))
    Q = spd(rng, r, 0.2)
    c = np.array([0.5, -0.3, 0.2])
    d = rng.standard_normal(p)
    y = rng.standard_normal((p, n)) * 2
    generic("mv_stationary", y, Z, H, T, R, Q, c, d)


def ar2():
    """AR(2) with SARIMAX: loglike at fixed params, a tight MLE, and the
    stationarity transform."""
    rng = np.random.default_rng(3)
    n, phi, sigma2 = 200, np.array([0.5, 0.3]), 1.2
    e = rng.standard_normal(n + 100) * np.sqrt(sigma2)
    y = np.zeros(n + 100)
    for t in range(2, n + 100):
        y[t] = phi[0] * y[t - 1] + phi[1] * y[t - 2] + e[t]
    y = y[100:]
    mod = sm.tsa.SARIMAX(y, order=(2, 0, 0), trend="n")
    params = np.r_[phi, sigma2]
    tight = mod.fit(disp=False, pgtol=1e-12, factr=10, maxiter=5000)
    tight = mod.filter(tight.params, cov_type="approx")
    x = np.array([0.7, -1.3, 2.0])
    from statsmodels.tsa.statespace.tools import constrain_stationary_univariate
    write_fixture("ar2", {
        "y": y[None, :], "params": params, "llf": mod.loglike(params),
        "llf_obs": mod.loglikeobs(params),
        "params_tight": tight.params, "llf_tight": tight.llf, "bse_tight": tight.bse,
        "transform_in": x, "transform_out": constrain_stationary_univariate(x)})
    print("ar2 llf =", mod.loglike(params), "tight", tight.params, tight.llf)


def read_fixture(name):
    arrays, lines, i = {}, open(f"{HERE}/{name}.txt").read().split("\n"), 0
    while i < len(lines) and lines[i]:
        head = lines[i].split()
        ndim = int(head[1])
        shape = tuple(int(x) for x in head[2:2 + ndim])
        size = int(np.prod(shape))
        arrays[head[0]] = np.array([float(x) for x in lines[i + 1:i + 1 + size]]).reshape(
            shape, order="F")
        i += 1 + size
    return arrays


def build_from_fixture(cls, f):
    """A statsmodels representation of class `cls` for a fixture's model."""
    from statsmodels.tsa.statespace.initialization import Initialization
    kinds = {INIT_KNOWN: "known", INIT_STATIONARY: "stationary", INIT_DIFFUSE: "diffuse"}
    p, n = f["y"].shape
    m, r = f["T"].shape[0], f["Q"].shape[0]
    ss = cls(k_endog=p, k_states=m, k_posdef=r)
    ss.bind(np.asfortranarray(f["y"]))
    for key, arr in [("design", "Z"), ("obs_cov", "H"), ("transition", "T"),
                     ("selection", "R"), ("state_cov", "Q")]:
        ss[key] = f[arr]
    ss["state_intercept"] = f["c"].reshape(m, -1)
    ss["obs_intercept"] = f["d"].reshape(p, -1)
    if "init_blocks" in f:
        init = Initialization(m)
        for first, last, kind in np.atleast_2d(f["init_blocks"]).astype(int):
            init.set((first - 1, last), kinds[kind])
        ss.initialization = init
    elif "a1" in f:
        ss.initialize_known(f["a1"], f["P1"])
    else:
        ss.initialize_stationary()
    return ss


def simulation_smoother_fixtures():
    """DK 4.9 simulation smoother on existing fixture models, with fixed N(0, 1)
    variates so the draws are deterministic."""
    from statsmodels.tsa.statespace.simulation_smoother import SimulationSmoother
    rng = np.random.default_rng(99)
    for base in ["mv_invariant", "mv_timevarying", "mv_missing", "mv_stationary",
                 "mv_diffuse", "uc_trend_seasonal_exact"]:
        # (Not uc_level_ar1_mixed: with a singular P_star, statsmodels scales the
        # initial variate by the variance instead of its square root.)
        f = read_fixture(base)
        p, n = f["y"].shape
        m, r = f["T"].shape[0], f["Q"].shape[0]
        ss = build_from_fixture(SimulationSmoother, f)
        u_eps = rng.standard_normal((p, n))
        u_eta = rng.standard_normal((r, n))
        u_init = rng.standard_normal(m)
        sim = ss.simulation_smoother()
        sim.simulate(measurement_disturbance_variates=u_eps.ravel(order="F"),
                     state_disturbance_variates=u_eta.ravel(order="F"),
                     initial_state_variates=u_init)
        write_fixture(f"sim_{base}", {
            "u_eps": u_eps, "u_eta": u_eta, "u_init": u_init,
            "generated_obs": sim.generated_obs, "generated_state": sim.generated_state,
            "state": sim.simulated_state, "eps": sim.simulated_measurement_disturbance,
            "eta": sim.simulated_state_disturbance})
        print("sim", base)


def smoothing_extras_fixtures():
    """DK 4.7 cross-time covariances of the smoothed state and DK 4.8 weights."""
    from statsmodels.tsa.statespace.tools import compute_smoothed_state_weights
    for base in ["mv_invariant", "mv_timevarying", "mv_missing", "mv_stationary",
                 "mv_diffuse", "uc_trend_seasonal_exact"]:
        f = read_fixture(base)
        n = f["y"].shape[1]
        # compute_smoothed_state_weights needs model results, so wrap the
        # representation in a parameterless MLEModel.
        ks = build_from_fixture(KalmanSmoother, f)
        mod = sm.tsa.statespace.MLEModel(f["y"].T, k_states=ks.k_states, k_posdef=ks.k_posdef)
        mod.ssm = ks
        res = mod.smooth([])
        weights, c_weights, prior_weights = compute_smoothed_state_weights(res)
        write_fixture(f"extras_{base}", {
            # Cov(alpha_t+1, alpha_t), t = 1..n-1
            "autocov": res.smoothed_state_autocov[:, :, :n - 1],
            # Cov(alpha_t, alpha_t+3), t = 1..n-3
            "cov_shift3": np.transpose(res.smoother_results._smoothed_state_autocovariance(3, 0, n - 3),
                                       (1, 2, 0)),
            # Weight of y_j,i in alphahat_t,k as (k, i, t, j), and likewise
            # c_j,l as (k, l, t, j) and a1_l as (k, l, t)
            "weights": np.transpose(weights, (2, 3, 0, 1)),
            "c_weights": np.transpose(c_weights, (2, 3, 0, 1)),
            "prior_weights": np.transpose(prior_weights, (1, 2, 0))})
        print("extras", base)


def concentrated_fixtures():
    """Scale concentrated out of the likelihood (filter_concentrated)."""
    y = pd.read_csv("data/nile.csv")["volume"].to_numpy(dtype=float)
    ymiss = y.copy()
    ymiss[[5, 6, 40]] = np.nan
    q = 1469.1 / 15099.0
    out = {"q": q}
    for init in ["diffuse", "known"]:
        for data, tag in [(y, "full"), (ymiss, "missing")]:
            ks = KalmanSmoother(1, 1, 1)
            ks.bind(data[None, :].copy())
            ks["design"], ks["obs_cov"], ks["transition"] = [[1.0]], [[1.0]], [[1.0]]
            ks["selection"], ks["state_cov"] = [[1.0]], [[q]]
            if init == "diffuse":
                ks.initialize_diffuse()
            else:
                ks.initialize_known([1000.0], [[100.0]])
            ks.filter_concentrated = True
            res = ks.filter()
            out[f"{init}_{tag}_llf"] = res.llf
            out[f"{init}_{tag}_scale"] = res.scale
    out["y_missing"] = ymiss[None, :]
    # p = 2, correlated H, diffuse, with a missing element in the diffuse period
    f = read_fixture("mv_diffuse")
    ks = build_from_fixture(KalmanSmoother, f)
    ks.filter_concentrated = True
    res = ks.filter()
    out["mv_diffuse_llf"] = res.llf
    out["mv_diffuse_scale"] = res.scale
    # Concentrated AR(2) MLE
    ya = read_fixture("ar2")["y"][0]
    fit = sm.tsa.SARIMAX(ya, order=(2, 0, 0), trend="n", concentrate_scale=True).fit(
        disp=False, pgtol=1e-12, factr=10)
    out.update({"ar2_params": fit.params, "ar2_scale": fit.scale, "ar2_llf": fit.llf,
                "ar2_aic": fit.aic, "ar2_bic": fit.bic})
    write_fixture("concentrated", out)
    print("concentrated", {k: v for k, v in out.items() if k.endswith(("llf", "scale"))})


def diagnostics_fixtures():
    """DK 7.5 diagnostics: standardized one-step errors and the statsmodels
    tests on them."""
    from scipy import stats
    y = pd.read_csv("data/nile.csv")["volume"].to_numpy(dtype=float)
    mod = sm.tsa.UnobservedComponents(y, "llevel", use_exact_diffuse=True)
    res = mod.smooth([15099.0, 1469.1])
    lb = res.test_serial_correlation("ljungbox", lags=10)[0]
    jb = res.test_normality("jarquebera")[0]
    het = res.test_heteroskedasticity("breakvar")[0]
    out = {"std_err": res.filter_results.standardized_forecasts_error,
           "lb_stat": lb[0], "lb_pvalue": lb[1], "jb": jb, "het": het,
           # special functions for the p-values
           "chi2_sf_in": [3.7, 2.0, 25.0, 17.5], "chi2_sf_out": [
               stats.chi2.sf(3.7, 2.0), stats.chi2.sf(25.0, 17.5)],
           "f_cdf_in": [0.8, 33.0, 32.0], "f_cdf_out": stats.f.cdf(0.8, 33.0, 32.0)}
    # Multivariate standardization with missing values
    f = read_fixture("mv_missing")
    ks = build_from_fixture(KalmanSmoother, f)
    r = ks.filter()
    out["mv_std_err"] = r.standardized_forecasts_error
    t = 2   # check the convention: lower Cholesky factor of F
    L = np.linalg.cholesky(r.forecasts_error_cov[:, :, t])
    print("mv std err convention (lower Cholesky) max diff:",
          np.abs(np.linalg.solve(L, r.forecasts_error[:, t]) - r.standardized_forecasts_error[:, t]).max())
    write_fixture("diagnostics", out)
    print("diagnostics", lb, jb, het)


def structural_fixtures():
    """Structural components (DK 3.2) against statsmodels' UnobservedComponents
    with fixed parameters and exact diffuse initialization. Recorded are the
    representation-invariant outputs: llf_obs and each component's smoothed
    signal."""
    from statsmodels.tsa.statespace.initialization import Initialization
    rng = np.random.default_rng(2024)

    def bsm_series(n, s):
        mu = np.cumsum(np.cumsum(0.05 * rng.standard_normal(n)) + 0.3 * rng.standard_normal(n))
        gam = np.tile(rng.standard_normal(s), n // s + 1)[:n]
        return 10 + mu + gam + 0.5 * rng.standard_normal(n)

    def record(name, mod, params, signals, init=None, extra=None):
        if init is not None:
            mod.ssm.initialization = init
        res = mod.smooth(params)
        out = {"y": mod.endog.T, "params": params, "llf": res.llf, "llf_obs": res.llf_obs}
        for key, comp in signals.items():
            out[key] = comp(res)["smoothed"]
        out.update(extra or {})
        write_fixture(name, out)
        print(name, "llf", res.llf)

    yq = bsm_series(80, 4)
    yq[[5, 17, 18, 40]] = np.nan
    record("uc_bsm_dummy",
           sm.tsa.UnobservedComponents(yq, "lltrend", seasonal=4, use_exact_diffuse=True),
           [0.25, 0.09, 0.0025, 0.01],
           {"level": lambda r: r.level, "trend": lambda r: r.trend,
            "seasonal": lambda r: r.seasonal})
    ym = bsm_series(96, 12)
    record("uc_llevel_trig12",
           sm.tsa.UnobservedComponents(ym, "llevel", freq_seasonal=[{"period": 12}],
                                       use_exact_diffuse=True),
           [0.25, 0.09, 0.004],
           {"level": lambda r: r.level, "seasonal": lambda r: r.freq_seasonal[0]})
    record("uc_smooth_trend_dummy",
           sm.tsa.UnobservedComponents(yq, "smooth trend", seasonal=4, use_exact_diffuse=True),
           [0.25, 0.0025, 0.01],
           {"level": lambda r: r.level, "trend": lambda r: r.trend,
            "seasonal": lambda r: r.seasonal})
    # Damped cycle, initialized as stationary (DK 3.2.4)
    yc = np.cumsum(0.3 * rng.standard_normal(100)) + 2 * np.sin(np.arange(100) * 2 * np.pi / 20) \
        + 0.5 * rng.standard_normal(100)
    mod = sm.tsa.UnobservedComponents(yc, "llevel", cycle=True, damped_cycle=True,
                                      stochastic_cycle=True, use_exact_diffuse=True)
    init = Initialization(mod.k_states)
    init.set((0, 1), "diffuse")
    init.set((1, 3), "stationary")
    record("uc_llevel_cycle", mod, [0.25, 0.09, 0.3, 2 * np.pi / 20, 0.9],
           {"level": lambda r: r.level, "cycle": lambda r: r.cycle}, init=init)
    # Regression and a step intervention with fixed coefficients (DK 3.2.5)
    n = 80
    x = rng.standard_normal(n)
    step = (np.arange(1, n + 1) >= 50).astype(float)
    yr = np.cumsum(0.3 * rng.standard_normal(n)) + 2 * x + 3 * step + 0.5 * rng.standard_normal(n)
    mod = sm.tsa.UnobservedComponents(yr, "llevel", exog=np.c_[x, step], mle_regression=False,
                                      use_exact_diffuse=True)
    res = mod.smooth([0.25, 0.09])
    write_fixture("uc_llevel_regression", {
        "y": yr[None, :], "x": np.c_[x, step], "params": [0.25, 0.09], "llf": res.llf,
        "llf_obs": res.llf_obs, "level": res.level["smoothed"],
        "beta": res.smoothed_state[1:, -1]})
    print("uc_llevel_regression llf", res.llf)
    # Random-walk coefficient (DK 3.15), from explicit matrices
    ks = KalmanSmoother(k_endog=1, k_states=2, k_posdef=2)
    ks.bind(yr[None, :].copy())
    Zt = np.zeros((1, 2, n))
    Zt[0, 0, :] = 1.0
    Zt[0, 1, :] = x
    ks["design"], ks["obs_cov"], ks["transition"] = Zt, [[0.25]], np.eye(2)
    ks["selection"], ks["state_cov"] = np.eye(2), np.diag([0.09, 0.01])
    ks.initialize_diffuse()
    res = ks.smooth()
    write_fixture("uc_llevel_rw_regression", {
        "y": yr[None, :], "x": x, "llf": res.llf_obs.sum(), "llf_obs": res.llf_obs,
        "alphahat": res.smoothed_state})
    print("uc_llevel_rw_regression llf", res.llf_obs.sum())


def arima_fixtures():
    """ARMA/ARIMA (DK 3.4, 5.6) against SARIMAX with the differenced states
    exact diffuse (SARIMAX's default is approximate diffuse)."""
    from statsmodels.tsa.statespace.initialization import Initialization
    rng = np.random.default_rng(77)
    n = 120
    e = rng.standard_normal(n + 50)
    x = np.zeros(n + 50)
    for t in range(2, n + 50):
        x[t] = 0.5 * x[t - 1] + 0.2 * x[t - 2] + e[t] + 0.4 * e[t - 1]
    x = x[50:]
    series = {"stationary": x, "i1": np.cumsum(x), "i2": np.cumsum(np.cumsum(x))}
    seas = np.cumsum(x) + np.tile([2.0, -1.0, 0.5, -1.5], n // 4)
    seas[[10, 11, 60]] = np.nan
    series["seasonal"] = seas

    def exact(mod, n_diffuse):
        if n_diffuse > 0:
            init = Initialization(mod.k_states)
            init.set((0, n_diffuse), "diffuse")
            init.set((n_diffuse, mod.k_states), "stationary")
            mod.ssm.initialization = init
        return mod

    for name, key, order, sorder, params in [
            ("arima_201", "stationary", (2, 0, 1), (0, 0, 0, 0), [0.5, 0.2, 0.4, 1.0]),
            ("arima_211", "i1", (2, 1, 1), (0, 0, 0, 0), [0.5, 0.2, 0.4, 1.0]),
            ("arima_221", "i2", (2, 2, 1), (0, 0, 0, 0), [0.5, 0.2, 0.4, 1.0]),
            ("arima_111_011_4", "seasonal", (1, 1, 1), (0, 1, 1, 4), [0.5, 0.3, -0.4, 1.0])]:
        y = series[key]
        mod = sm.tsa.SARIMAX(y, order=order, seasonal_order=sorder, simple_differencing=False)
        n_diffuse = order[1] + sorder[1] * sorder[3]
        exact(mod, n_diffuse)
        res = mod.smooth(params)
        out = {"y": y[None, :], "order": list(order) + list(sorder), "params": params,
               "llf": res.llf, "llf_obs": res.llf_obs, "alphahat": res.smoothed_state}
        if name == "arima_211":
            fit = exact(sm.tsa.SARIMAX(y, order=order, simple_differencing=False),
                        n_diffuse).fit(disp=False, pgtol=1e-12, factr=10, maxiter=5000)
            # SARIMAX burns the first d observations in `llf`; the exact
            # diffuse log likelihood is the sum over all of them.
            out.update({"params_tight": fit.params, "llf_tight": fit.llf_obs.sum()})
        write_fixture(name, out)
        print(name, "llf", res.llf)

    # Regression with ARMA(1,1) errors and a constant (DK 3.6.2, 5.6.4)
    z = rng.standard_normal(n)
    y = 3.0 + 1.5 * z + x
    mod = sm.tsa.SARIMAX(y, exog=np.c_[np.ones(n), z], order=(1, 0, 1), mle_regression=False)
    # SARIMAX orders the states ARMA first, then the regression coefficients.
    init = Initialization(mod.k_states)
    init.set((0, 2), "stationary")
    init.set((2, 4), "diffuse")
    mod.ssm.initialization = init
    res = mod.smooth([0.6, 0.3, 1.0])
    write_fixture("arima_regression", {"y": y[None, :], "z": z, "params": [0.6, 0.3, 1.0],
                                       "llf": res.llf, "llf_obs": res.llf_obs,
                                       "beta": res.smoothed_state[2:, -1]})
    print("arima_regression llf", res.llf)


def spline_fixtures():
    """Cubic smoothing spline at irregular points (DK 3.9.2), from scipy."""
    from scipy.interpolate import make_smoothing_spline
    rng = np.random.default_rng(5)
    x = np.sort(rng.uniform(0, 10, 60))
    y = np.sin(x) + 0.3 * rng.standard_normal(60)
    lam = 0.5
    spl = make_smoothing_spline(x, y, lam=lam)
    write_fixture("smoothing_spline", {"x": x, "y": y[None, :], "lam": lam, "fitted": spl(x)})
    print("smoothing_spline")


if __name__ == "__main__":
    nile()
    multivariate(False)
    multivariate(True)
    multivariate_missing()
    multivariate_stationary()
    ar2()
    diffuse_models()
    simulation_smoother_fixtures()
    smoothing_extras_fixtures()
    concentrated_fixtures()
    diagnostics_fixtures()
    structural_fixtures()
    arima_fixtures()
    spline_fixtures()
