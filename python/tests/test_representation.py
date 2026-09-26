"""The Python layer against the statsmodels fixtures used by the Fortran tests
(test/fixtures, written by test/fixtures/make_fixtures.py)."""

import numpy as np
import pytest

import ssfortran as ss

from _fixtures import close, read_fixture, rep_from_fixture


# statsmodels' smoother is wrong in the diffuse period of the time-varying
# diffuse fixture (docs/statsmodels_differences.md), so it is left out.
CASES = ["nile_llevel_known", "nile_llevel_exact", "nile_llevel_missing", "nile_lltrend_exact",
         "mv_invariant", "mv_timevarying", "mv_missing", "mv_stationary", "mv_diffuse",
         "uc_trend_seasonal_exact", "uc_level_ar1_mixed"]


@pytest.mark.parametrize("name", CASES)
def test_filter_and_smoother(name):
    fx = read_fixture(name)
    rep = rep_from_fixture(fx)
    res = rep.smooth()
    f = res.filter
    n = rep.nobs
    close(f.llf_obs, fx["llf_obs"])
    close([rep.loglike()], [f.llf])
    close(f.predicted_state, fx["a"])
    close(f.predicted_state_cov, fx["P"])
    close(res.smoothed_state, fx["alphahat"])
    close(res.smoothed_state_cov, fx["V"])
    # v and F: in diffuse periods statsmodels reports them in transformed
    # coordinates for p > 1, so compare those after the diffuse period only.
    late = np.arange(n) >= (f.nobs_diffuse if rep.k_endog > 1 else 0)
    close(f.forecasts_error, fx["v"], np.broadcast_to(late, fx["v"].shape))
    close(f.forecasts_error_cov, fx["F"], np.broadcast_to(late, fx["F"].shape))


def test_options_and_errors():
    fx = read_fixture("mv_invariant")
    rep = rep_from_fixture(fx)
    full = rep_from_fixture(fx)
    full.tol_steady = -1.0
    assert full.filter().t_steady == 0
    close([rep.loglike()], [full.loglike()])
    rep.filter_method = ss.FILTER_UNIVARIATE
    close([rep.loglike()], [full.loglike()])

    with pytest.raises(ValueError):
        rep["design"] = np.ones((5, 5))
    bad = ss.Representation(np.ones(10), k_states=1)
    bad["design"] = [[1.0]]
    bad["transition"] = [[1.0]]
    bad["selection"] = [[1.0]]
    bad["state_cov"] = [[1.0]]
    bad["obs_cov"] = [[1.0]]
    with pytest.raises(ss.StateSpaceError) as e:
        bad.initialize_stationary()   # a unit root has no stationary distribution
        bad.loglike()
    assert e.value.code == 6


def test_version():
    assert ss.__version__ == "0.1.0"
