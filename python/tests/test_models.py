"""Built-in models and estimation through the Python layer."""

import numpy as np
import pytest
import ssfortran as ss
from _fixtures import close, read_fixture


def test_bsm_dummy_matches_statsmodels():
    fx = read_fixture("uc_bsm_dummy")
    mod = ss.StructuralModel(fx["y"].T, [ss.Irregular(), ss.Trend(), ss.Seasonal(4)])
    assert mod.param_names == [
        "sigma2.irregular",
        "sigma2.level",
        "sigma2.slope",
        "sigma2.seasonal",
    ]
    res = mod.smooth(fx["params"])
    close(res.filter.llf_obs, fx["llf_obs"])
    close(res.smoothed_state[0], fx["level"])
    close(res.smoothed_state[1], fx["trend"])


def test_trig_seasonal():
    fx = read_fixture("uc_llevel_trig12")
    mod = ss.StructuralModel(
        fx["y"].T, [ss.Irregular(), ss.Level(), ss.Seasonal(12, form="trig")]
    )
    res = mod.smooth(fx["params"])
    close(res.filter.llf_obs, fx["llf_obs"])
    assert mod.k_states == 12


def test_arima_matches_sarimax():
    fx = read_fixture("arima_211")
    mod = ss.StructuralModel(fx["y"].T, [ss.ARIMA(order=(2, 1, 1))])
    res = mod.smooth(fx["params"])
    close(res.filter.llf_obs, fx["llf_obs"])
    close(
        res.smoothed_state, fx["alphahat"].reshape(res.smoothed_state.shape, order="F")
    )


def test_nile_mle():
    fx = read_fixture("nile_llevel_mle_exact")
    y = read_fixture("nile_llevel_exact")["y"].T
    mod = ss.StructuralModel(y, [ss.Irregular(), ss.Level()])
    res = mod.fit(factr=10.0, pgtol=1e-9)
    assert res.converged and res.analytic_gradient
    assert abs(res.llf - fx["llf_tight"][0]) < 1e-6
    assert np.allclose(res.params, fx["params_tight"], rtol=1e-4)
    assert np.allclose(res.bse, fx["bse_tight"], rtol=1e-3)
    assert "sigma2.level" in res.summary()


def test_representation_at_params_and_options():
    y = read_fixture("nile_llevel_exact")["y"].T
    mod = ss.StructuralModel(y, [ss.Irregular(), ss.Level()])
    params = [15099.0, 1469.1]
    rep = mod.representation(params)
    assert rep["obs_cov"][0, 0] == 15099.0
    assert abs(rep.loglike() - mod.loglike(params)) < 1e-9
    mod.ssm.filter_method = ss.FILTER_UNIVARIATE
    assert abs(mod.loglike(params) - rep.loglike()) < 1e-8
    x = mod.untransform_params(params)
    assert np.allclose(mod.transform_params(x), params)


def test_fit_many_matches_fit():
    rng = np.random.default_rng(3)
    ys = [
        np.cumsum(0.3 * rng.standard_normal(80)) + rng.standard_normal(80)
        for _ in range(6)
    ]
    comps = [ss.Irregular(), ss.Level()]
    many = ss.fit_many([ss.StructuralModel(y, comps) for y in ys])
    for y, r in zip(ys, many):
        one = ss.StructuralModel(y, comps).fit()
        assert np.array_equal(r.params, one.params) and r.llf == one.llf


def test_errors():
    with pytest.raises(ValueError):
        ss.StructuralModel(np.ones(10), [ss.Level(cov="bogus")])
    with pytest.raises(TypeError):
        ss.StructuralModel(np.ones(10), [object()])


def test_cycle_component():
    fx = read_fixture("uc_llevel_cycle")
    mod = ss.StructuralModel(fx["y"].T, [ss.Irregular(), ss.Level(), ss.Cycle()])
    assert mod.param_names[-2:] == ["frequency.cycle", "damping.cycle"]
    res = mod.smooth(fx["params"])
    close(res.filter.llf_obs, fx["llf_obs"])
    close(res.smoothed_state[1], fx["cycle"])  # Z picks the first cycle state


def test_regression_components():
    fx = read_fixture("uc_llevel_regression")
    x = fx["x"].copy()
    x[:, 1] = (np.arange(x.shape[0]) >= 49).astype(float)  # step at t = 50
    mod = ss.StructuralModel(fx["y"].T, [ss.Irregular(), ss.Level(), ss.Regression(x)])
    res = mod.smooth(fx["params"])
    close(res.filter.llf_obs, fx["llf_obs"])
    close(res.smoothed_state[0], fx["level"])
    close(res.smoothed_state[1:3, -1], fx["beta"])

    fx = read_fixture("uc_llevel_rw_regression")
    mod = ss.StructuralModel(
        fx["y"].T,
        [ss.Irregular(), ss.Level(), ss.Regression(fx["x"], random_walk=True)],
    )
    assert mod.param_names[-1] == "sigma2.beta.1"
    res = mod.smooth([0.25, 0.09, 0.01])
    close(res.filter.llf_obs, fx["llf_obs"])
    close(
        res.smoothed_state, fx["alphahat"].reshape(res.smoothed_state.shape, order="F")
    )


def test_continuous_spline():
    fx = read_fixture("smoothing_spline")
    lam = fx["lam"].sum()
    mod = ss.StructuralModel(fx["y"].T, [ss.Irregular(), ss.ContinuousTrend(fx["x"])])
    res = mod.smooth([1.0, 1.0 / lam])
    close(res.smoothed_state[0], fx["fitted"])


def test_common_levels_loadings():
    y = read_fixture("mv_invariant")["y"].T  # (n, 2)
    a2 = 0.7
    mod = ss.StructuralModel(
        y,
        [
            ss.Irregular(),
            ss.Level(),
            ss.Regression(np.ones(y.shape[0]), series=1, at_observations=True),
        ],
        loading=[[1.0], [0.0]],
        loading_free=[[False], [True]],
    )
    assert mod.k_params == 4
    rep = ss.Representation(y, k_states=2, k_posdef=1)
    rep["design"] = [[1.0, 0.0], [a2, 1.0]]
    rep["transition"] = np.eye(2)
    rep["selection"] = [[1.0], [0.0]]
    rep["state_cov"] = [[0.25]]
    rep["obs_cov"] = np.diag([1.0, 2.0])
    rep.initialize_diffuse()
    close([mod.loglike([1.0, 2.0, 0.25, a2])], [rep.loglike()])


def local_level_q(y):
    """Local level with H = 1 and Q = q, for a concentrated scale."""
    rep = ss.Representation(y, k_states=1)
    rep["design"] = [[1.0]]
    rep["transition"] = [[1.0]]
    rep["selection"] = [[1.0]]
    rep["obs_cov"] = [[1.0]]
    rep.initialize_diffuse()
    return (
        ss.MappedModel(rep, 1, ["q"], start_params=[0.1])
        .map(0, "state_cov", 0, 0)
        .constrain(0, "positive")
    )


def test_concentrated_scale():
    cx = read_fixture("concentrated")
    y = read_fixture("nile_llevel_known")["y"].ravel()
    mod = local_level_q(y)
    mod.concentrate_scale = True
    close([mod.loglike(cx["q"])], cx["diffuse_full_llf"])
    close([mod.scale], cx["diffuse_full_scale"])

    # AR(2) with the variance concentrated out, as test_concentrated_ar2_fit
    y2 = read_fixture("ar2")["y"].ravel()
    rep = ss.Representation(y2, k_states=2, k_posdef=1)
    rep["design"] = [[1.0, 0.0]]
    rep["transition"] = [[0.0, 0.0], [1.0, 0.0]]
    rep["selection"] = [[1.0], [0.0]]
    rep["state_cov"] = [[1.0]]
    rep.initialize_stationary()
    ar2 = (
        ss.MappedModel(rep, 2, start_params=[0.0, 0.0])
        .map(0, "transition", 0, 0)
        .map(1, "transition", 0, 1)
        .constrain([0, 1], "stationary")
    )
    ar2.concentrate_scale = True
    res = ar2.fit(factr=10.0, pgtol=1e-9)
    assert abs(res.llf - cx["ar2_llf"][0]) < 1e-8
    assert np.allclose(res.params, cx["ar2_params"], rtol=1e-4)
    assert np.isclose(res.scale, cx["ar2_scale"][0], rtol=1e-4)
    assert np.allclose(
        [res.aic, res.bic], [cx["ar2_aic"][0], cx["ar2_bic"][0]], rtol=1e-6
    )
