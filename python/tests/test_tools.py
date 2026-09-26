"""Forecasting, simulation and diagnostics through the Python layer."""

import numpy as np

import ssfortran as ss
from _fixtures import close, read_fixture, rep_from_fixture


def test_diagnostics_match_statsmodels():
    rep = rep_from_fixture(read_fixture("nile_llevel_exact"))
    dx = read_fixture("diagnostics")
    f = rep.filter()
    d = f.diagnostic_start
    assert d == 1 and np.isnan(f.standardized_forecasts_error[0, 0])
    e = f.standardized_forecasts_error[0, d:]
    close(e, dx["std_err"][0, d:])
    stat, pvalue = ss.diagnostics.ljung_box(e, 10)
    close(stat, dx["lb_stat"])
    close(pvalue, dx["lb_pvalue"])
    close(ss.diagnostics.jarque_bera(e), dx["jb"])
    close(ss.diagnostics.breakvar(e), dx["het"])


def test_simulation_matches_statsmodels():
    # statsmodels' simulated series are comparable given the same variates
    # for the fixture with missing data (as in test/test_simsmooth.f90)
    rep = rep_from_fixture(read_fixture("mv_missing"))
    sx = read_fixture("sim_mv_missing")
    y, alpha, eps, eta = rep.simulate(variates=(sx["u_init"], sx["u_eps"], sx["u_eta"]))
    close(y, sx["generated_obs"].reshape(y.shape, order="F"))
    close(alpha, sx["generated_state"][:, :y.shape[1]])

    rep = rep_from_fixture(read_fixture("mv_invariant"))
    sx = read_fixture("sim_mv_invariant")
    v = (sx["u_init"], sx["u_eps"], sx["u_eta"])
    state, eps, eta = rep.simulation_smoother(variates=v)
    close(state, sx["state"].reshape(state.shape, order="F"))
    close(eps, sx["eps"].reshape(eps.shape, order="F"))
    close(eta, sx["eta"].reshape(eta.shape, order="F"))
    a1, _, _ = rep.simulation_smoother(rng=1)
    a2, _, _ = rep.simulation_smoother(rng=1)
    assert np.array_equal(a1, a2)


def test_forecast_equals_missing_observations():
    fx = read_fixture("mv_invariant")
    rep = rep_from_fixture(fx)
    mean, cov = rep.forecast(5)
    y = np.hstack([fx["y"], np.full((fx["y"].shape[0], 5), np.nan)])
    fx_long = dict(fx, y=y)
    ext = rep_from_fixture(fx_long).filter()
    n = fx["y"].shape[1]
    close(mean, ext.forecasts[:, n:])
    close(cov, ext.forecasts_error_cov[:, :, n:])


def test_steady_state_local_level():
    rep = rep_from_fixture(read_fixture("nile_llevel_known"))
    rep["obs_cov"] = [[15099.0]]
    rep["state_cov"] = [[1469.1]]
    P, F = rep.steady_state()
    q = 1469.1 / 15099.0
    x = (q + np.sqrt(q * q + 4 * q)) / 2
    close([P[0, 0], F[0, 0]], [15099.0 * x, 15099.0 * (1 + x)])
