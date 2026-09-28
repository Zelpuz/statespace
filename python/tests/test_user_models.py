"""Models defined from Python: MappedModel (declarative) and MLEModel
(callbacks), checked against the built-in models."""

import numpy as np
import pytest
import ssfortran as ss
from _fixtures import read_fixture

NILE = read_fixture("nile_llevel_exact")["y"].ravel()


def local_level(y, **kwargs):
    """The local level model as a MappedModel, set up directly."""
    mod = ss.MappedModel(
        y, k_states=1, k_params=2, start_params=[y.var() / 2] * 2, **kwargs
    )
    mod["design"] = mod["transition"] = mod["selection"] = [[1.0]]
    mod.initialize_diffuse()
    return (
        mod.map(0, "obs_cov", 0, 0)
        .map(1, "state_cov", 0, 0)
        .constrain([0, 1], "positive")
    )


def local_level_rep(y):
    rep = ss.Representation(y, k_states=1)
    rep["design"] = [[1.0]]
    rep["transition"] = [[1.0]]
    rep["selection"] = [[1.0]]
    rep.initialize_diffuse()
    return rep


def test_mapped_local_level_matches_structural():
    mod = local_level(NILE, param_names=["sigma2.irregular", "sigma2.level"])
    ref = ss.StructuralModel(NILE, [ss.Irregular(), ss.Level()])
    assert mod.loglike([15099.0, 1469.1]) == pytest.approx(
        ref.loglike([15099.0, 1469.1]), rel=1e-12
    )
    r1 = mod.fit(factr=10.0, pgtol=1e-9)
    r2 = ref.fit(factr=10.0, pgtol=1e-9)
    assert r1.param_names == r2.param_names
    assert np.allclose(r1.params, r2.params, rtol=1e-5) and r1.llf == pytest.approx(
        r2.llf
    )
    assert r1.analytic_gradient


def test_mapped_from_representation():
    """A Representation in place of the data: by keyword, or with the 0.1
    positional arguments, which warn."""
    ref = local_level(NILE).loglike([15099.0, 1469.1])
    kw = ss.MappedModel(local_level_rep(NILE), k_params=2, param_names=["a", "b"])
    with pytest.warns(DeprecationWarning):
        pos = ss.MappedModel(local_level_rep(NILE), 2, ["a", "b"], [1.0, 2.0])
    assert pos.param_names == ["a", "b"] and np.allclose(pos.start_params, [1, 2])
    for m in (kw, pos):
        m.map(0, "obs_cov", 0, 0).map(1, "state_cov", 0, 0)
        assert m.loglike([15099.0, 1469.1]) == pytest.approx(ref, rel=1e-12)


def test_mapped_setup_errors():
    with pytest.raises(TypeError, match="k_states"):
        ss.MappedModel(NILE, k_params=2)
    with pytest.raises(TypeError, match="k_params"):
        ss.MappedModel(NILE, k_states=1)
    mod = local_level(NILE)
    assert mod["transition"].shape == (1, 1)
    with pytest.raises(ss.StateSpaceError):  # period 5 of a time-invariant matrix
        mod.map(0, "design", 0, 0, t=5)


def test_mapped_covariance_blocks():
    fx = read_fixture("mv_invariant")
    y = fx["y"].T[:, :2]
    mod = ss.MappedModel(y, k_states=2, k_params=6)
    mod["design"] = mod["transition"] = mod["selection"] = np.eye(2)
    mod.initialize_diffuse()
    mod.cov(0, "obs_cov", 0, 2).cov(3, "state_cov", 0, 2)
    ref = ss.StructuralModel(y, [ss.Irregular(cov="full"), ss.Level(cov="full")])
    params = [1.0, 0.3, 2.0, 0.5, 0.1, 0.4]
    assert mod.loglike(params) == pytest.approx(ref.loglike(params), rel=1e-12)
    x = mod.untransform_params(params)
    assert np.allclose(x, ref.untransform_params(params))
    assert np.allclose(mod.transform_params(x), params)


def test_mapped_arma_matches_arima():
    y = read_fixture("arima_201")["y"].ravel()
    mod = ss.MappedModel(
        y, k_states=2, k_posdef=1, k_params=3, start_params=[0.0, 0.0, y.var()]
    )
    mod["design"] = [[1.0, 0.0]]
    mod["transition"] = [[0.0, 1.0], [0.0, 0.0]]
    mod["selection"] = [[1.0], [0.0]]
    mod.initialize_stationary()
    mod.map(0, "transition", 0, 0).map(1, "selection", 1, 0).map(2, "state_cov", 0, 0)
    mod.constrain(0, "stationary").constrain(1, "invertible").constrain(2, "positive")
    ref = ss.StructuralModel(y, [ss.ARIMA(order=(1, 0, 1))])
    params = [0.5, 0.3, 1.2]
    assert mod.loglike(params) == pytest.approx(ref.loglike(params), rel=1e-12)
    assert np.allclose(mod.untransform_params(params), ref.untransform_params(params))


class LocalLevel(ss.MLEModel):
    def __init__(self, y):
        super().__init__(y, k_states=1, k_params=2)
        self["design"] = [[1.0]]
        self["transition"] = [[1.0]]
        self["selection"] = [[1.0]]
        self.initialize_diffuse()
        self.calls = 0

    @property
    def param_names(self):
        return ["sigma2.irregular", "sigma2.level"]

    @property
    def start_params(self):
        return np.array([NILE.var() / 2] * 2)

    def transform_params(self, x):
        return np.exp(2 * np.asarray(x))

    def untransform_params(self, p):
        return 0.5 * np.log(p)

    def update(self, params):
        self.calls += 1
        if params[0] <= 0:
            raise ValueError("negative variance")
        self["obs_cov"] = [[params[0]]]
        self["state_cov"] = [[params[1]]]


def test_mlemodel_matches_structural():
    mod = LocalLevel(NILE)
    ref = ss.StructuralModel(NILE, [ss.Irregular(), ss.Level()])
    assert mod.loglike([15099.0, 1469.1]) == pytest.approx(
        ref.loglike([15099.0, 1469.1]), rel=1e-12
    )
    r1 = mod.fit(factr=10.0, pgtol=1e-9)
    r2 = ref.fit(factr=10.0, pgtol=1e-9)
    assert np.allclose(r1.params, r2.params, rtol=1e-5)
    assert r1.param_names == ["sigma2.irregular", "sigma2.level"]
    assert mod.calls > 0
    assert np.allclose(
        mod.transform_params(mod.untransform_params([2.0, 3.0])), [2.0, 3.0]
    )
    sm = r1.smooth()
    assert sm.smoothed_state.shape == (1, NILE.size)


def test_mlemodel_errors_propagate():
    mod = LocalLevel(NILE)
    with pytest.raises(ValueError, match="negative variance"):
        mod.loglike([-1.0, 1.0])
    assert mod.loglike([15099.0, 1469.1]) < 0  # usable afterwards


def test_fit_many_models():
    mapped = [local_level(NILE[i * 10 :]) for i in range(3)]
    res = ss.fit_many(mapped)
    assert all(r is not None and r.converged for r in res)
    with pytest.raises(TypeError):
        ss.fit_many([LocalLevel(NILE)])


class FixedH(ss.MLEModel):
    """H fixed in __init__; update sets only Q."""

    def __init__(self, y):
        super().__init__(y, k_states=1, k_params=1)
        self["design"] = [[1.0]]
        self["transition"] = [[1.0]]
        self["selection"] = [[1.0]]
        self["obs_cov"] = [[15099.0]]
        self.ssm.initialize_diffuse()

    def update(self, params):
        if params[0] < 0:
            raise ValueError("negative variance")
        self["state_cov"] = [[params[0]]]


def test_callback_failure_does_not_poison_the_model():
    mod = FixedH(NILE)
    good = mod.loglike([1469.1])
    with pytest.raises(ValueError):
        mod.loglike([-1.0])
    assert mod.loglike([1469.1]) == pytest.approx(good, rel=1e-14)
    assert mod["obs_cov"][0, 0] == 15099.0


def test_fit_many_rejects_callbacks_in_the_library():
    import ctypes

    mod = FixedH(NILE)
    hs = (ctypes.c_void_p * 1)(mod._h.value)
    fhs = (ctypes.c_void_p * 1)()
    infos = (ctypes.c_int * 1)()
    from ssfortran._lib import lib

    code = lib.ss_fit_many(1, hs, 100, 10, 1e7, 1e-5, 0, 0, fhs, infos)
    assert code == 4  # SS_ERR_UNSUPPORTED
